"""Post hoc sensitivity reanalysis of the frozen record frame. No network calls.
Model SE/CI values are not converted to individual SDs in raw-summary analyses.
"""
from pathlib import Path
import argparse, collections, gzip, hashlib, json, math, re
import numpy as np
import pandas as pd
import families as F

EXPECTED_SHA='6c2d72d295407261d769399f43ca02c3288f2b8317af944dc6395287c7bd826f'

def num(v):
    try:
        f=float(v)
        return f if math.isfinite(f) else None
    except (TypeError,ValueError): return None

def norm(x): return re.sub(r'\s+',' ',str(x or '').lower()).strip()

def days(text):
    t=str(text or '')
    if re.search(r'\b(male|female|subgroup|age|region|responder|baseline)\b',t,re.I): return None
    mult=lambda u: 7 if u.lower().startswith('week') else 30.4375 if u.lower().startswith('month') else 365.25 if u.lower().startswith('year') else 1
    v=[float(n)*mult(u) for n,u in re.findall(r'(\d+(?:\.\d+)?)\s*(days?|weeks?|months?|years?)\b',t,re.I)]
    v += [float(n)*mult(u) for u,n in re.findall(r'\b(days?|weeks?|months?|years?)\s*(\d+(?:\.\d+)?)\b',t,re.I)]
    return max(v) if v else None

def pick_class(o):
    cs=o.get('classes') or []
    if len(cs)==1: return 0,cs[0],'single_class'
    ts=[days(c.get('title')) for c in cs]
    if not ts or any(t is None for t in ts) or len(set(ts))!=len(ts): return None,None,'ambiguous_classes'
    j=int(np.argmax(ts)); return j,cs[j],'latest_explicit_time'

def arm(o,c,gid):
    cats=c.get('categories') or []
    if len(cats)!=1: return None
    ms=[m for m in cats[0].get('measurements') or [] if m.get('groupId')==gid]
    if len(ms)!=1: return None
    m=ms[0]; mean=num(m.get('value')); n=None; ns='outcome_participant_denominator'
    for d in o.get('denoms') or []:
        if str(d.get('units') or '').lower().startswith('participant'):
            nn=[num(x.get('value')) for x in d.get('counts') or [] if x.get('groupId')==gid]
            if nn: n=nn[0]; break
    gs=o.get('groups') or []
    mm=re.search(r'\bn\s*=\s*(\d+(?:\s*,\s*\d+)*)',c.get('title') or '',re.I)
    if mm and len(mm.group(1).split(','))==len(gs):
        ns='class_label_ordered_groups'; n=float(mm.group(1).split(',')[[g.get('id') for g in gs].index(gid)].strip())
    elif len(o.get('classes') or [])>1: ns='time_specific_n_unverified'
    if n is None or n<3 or mean is None: return None
    spread=num(m.get('spread'))
    sd=spread if o.get('paramType')=='MEAN' and o.get('disp')=='Standard Deviation' and spread is not None and spread>0 else None
    return dict(n=n,n_source=ns,mean=mean,spread=spread,lower=num(m.get('lowerLimit')),upper=num(m.get('upperLimit')),raw_sd=sd)

def proto(s,g):
    a=[a for a in s.get('arms') or [] if norm(a.get('label'))==norm(g.get('title'))]
    return a[0].get('type') if len(a)==1 else 'UNMATCHED'

def parse_p(v):
    if v is None: return None,'missing'
    t=str(v).strip().replace(' ','').replace('p','').replace('P','')
    m=re.fullmatch(r'(<=|>=|[<>=≤≥]?)(\d*\.?\d+(?:[eE]-?\d+)?)',t)
    if not m: return None,'unparsed'
    op,x=m.group(1),float(m.group(2))
    if not 0<=x<=1: return None,'out_of_range'
    if op in ('>','>=','≥'): sig=False if x>=.05 else None
    elif op=='<': sig=True if x<=.05 else None
    elif op in ('<=','≤'): sig=True if x<.05 else None
    else: sig=x<.05
    return sig,'determinate' if sig is not None else 'interval_crosses_threshold'

def p_audit(o,pg,active):
    if len(active)!=1: return 'multiple_or_missing_experimental_arms',None
    aid=active[0][0]['id']; pid=pg['id']
    aa=[a for a in o.get('analyses') or [] if len(a.get('g') or [])==2 and set(a.get('g') or [])=={aid,pid}]
    if not aa: return 'no_exact_two_group_analysis',None
    if len(aa)!=1: return 'multiple_exact_group_analyses',None
    if len(o.get('classes') or [])!=1: return 'analysis_timepoint_not_identified',None
    txt=' '.join(str(aa[0].get(k) or '') for k in ('desc','mdesc','sm'))
    if re.search(r'non.?infer|equivalence|one.?sided|adjust|multiplic',txt,re.I): return 'requires_analysis_adjudication',None
    sig,status=parse_p(aa[0].get('pval')); return 'single_exact_contrast_'+status,sig

def decorate(raw,records):
    wanted=collections.defaultdict(list)
    for r in records.to_dict('records'): wanted[r['nct']].append(r)
    out=[]; pairs=[]
    for line in gzip.open(raw,'rt',encoding='utf8'):
        s=json.loads(line)
        for base in wanted.get(s['nct'],[]):
            r=dict(base); hits=[(i,o) for i,o in enumerate(s.get('om') or []) if o.get('title')==r['title']]
            if not hits: raise RuntimeError('Frozen record outcome missing from verified snapshot')
            oi,o=hits[0]; gs=o.get('groups') or []; pgs=[g for g in gs if g.get('title')==r['group_title']]
            ci,c,cr=pick_class(o); vv=arm(o,c,pgs[0]['id']) if c is not None and len(pgs)==1 else None
            pp=[g for g in gs if F.PLACEBO.search(g.get('title') or '')]
            alternatives=[oo for oo in s.get('om') or [] if F.classify_title(oo.get('title') or '')==r['family'] and F.CHANGE_WORD.search(oo.get('title') or '') and oo.get('type')==r['omtype'] and oo.get('paramType') in F.PARAM_OK and oo.get('disp') in ('Standard Deviation','Standard Error','95% Confidence Interval') and any(F.PLACEBO.search(g.get('title') or '') for g in oo.get('groups') or [])]
            r.update(outcome_index=oi,outcome_title_occurrences=len(hits),param_type=o.get('paramType'),dispersion_type=o.get('disp'),outcome_population=o.get('pop'),outcome_description=o.get('desc'),class_index=ci,class_rule=cr,class_title=c.get('title') if c else None,number_classes=len(o.get('classes') or []),category_count=len(c.get('categories') or []) if c else 0,number_placebo_groups=len(pp),same_priority_outcomes=len(alternatives))
            why=[]
            if len(hits)!=1: why.append('duplicate_outcome_title')
            if len(pgs)!=1: why.append('ambiguous_placebo_title')
            if len(pp)!=1: why.append('multiple_placebo_groups')
            if c is None: why.append(cr)
            if vv is None: why.append('measurement_or_n_ambiguous')
            if vv is not None:
                r.update({('extracted_'+k if k not in ('raw_sd','n_source') else k):v for k,v in vv.items()}); r['placebo_group_id']=pgs[0]['id']
                if vv['raw_sd'] is None: why.append('not_arithmetic_mean_with_reported_SD')
                if vv['n_source']=='time_specific_n_unverified': why.append('time_specific_n_unverified')
            r['raw_eligible']=not why; r['raw_exclusion_reason']=';'.join(why)
            sign=r['resp']/r['mean'] if r['mean']!=0 else (-F.FAMILIES[r['family']][1] if r['sign_rule']=='change' else 1)
            r['improvement_sign']=sign
            if not why:
                y=sign*vv['mean']/vv['raw_sd']; r.update(raw_y=y,raw_v=1/vv['n']+y*y/(2*vv['n']))
            active=[]; allraw=[]; metadata=[]; experimental_groups=[]
            if c is not None:
                for g in gs:
                    t=g.get('title') or ''
                    if F.PLACEBO.search(t) or re.search(r'\bthen\b|->|→|followed\s+by|/\s*placebo|placebo\s*/|\btotal\b|\boverall\b|all\s+(participants|subjects|patients)',t,re.I): continue
                    v=arm(o,c,g['id']); typ=proto(s,g)
                    metadata.append(dict(group_id=g['id'],title=t,protocol_type=typ,raw_statistics_usable=bool(v and v['raw_sd'])))
                    if typ=='EXPERIMENTAL': experimental_groups.append((g,v))
                    if not v or not v['raw_sd'] or v['n_source']=='time_specific_n_unverified': continue
                    allraw.append((g,v))
                    if typ=='EXPERIMENTAL': active.append((g,v))
            r.update(nonplacebo_arm_metadata=json.dumps(metadata,ensure_ascii=False),experimental_raw_arms=len(active),all_nonplacebo_raw_arms=len(allraw))
            ps,sig=p_audit(o,pgs[0],experimental_groups) if len(pgs)==1 else ('ambiguous_placebo',None)
            r.update(p_audit_status=ps,reported_contrast_p_below_005=sig)
            if not why:
                for label,aa in [('experimental_only',active),('all_nonplacebo',allraw)]:
                    if not aa: continue
                    na=sum(v['n'] for _,v in aa); ma=sum(v['n']*v['mean'] for _,v in aa)/na
                    sa=math.sqrt(sum((v['n']-1)*v['raw_sd']**2+v['n']*(v['mean']-ma)**2 for _,v in aa)/(na-1))
                    n0,m0,s0=vv['n'],vv['mean'],vv['raw_sd']; df=na+n0-2
                    sp=math.sqrt(((na-1)*sa**2+(n0-1)*s0**2)/df); d=sign*(ma-m0)/sp; J=1-3/(4*df-1)
                    pairs.append(dict(nct=r['nct'],family=r['family'],cat=r['cat'],year=r['year'],omtype=r['omtype'],pairing=label,active_groups=json.dumps([g['id'] for g,_ in aa]),n_active=na,n_placebo=n0,mean_active=ma,mean_placebo=m0,sd_active=sa,sd_placebo=s0,y=J*d,v=J*J*((na+n0)/(na*n0)+d*d/(2*df)),active_count=len(aa),same_priority_outcomes=len(alternatives),outcome_index=oi,class_index=ci,title=r['title'],pure_placebo=not r['combo']))
            out.append(r)
    return pd.DataFrame(out),pd.DataFrame(pairs)

def fit(d,yc='y',vc='v',B=0,seed=2026102401):
    d=d.loc[np.isfinite(d[yc])&np.isfinite(d[vc])&(d[vc]>0)].reset_index(drop=True); ng=d.nct.nunique()
    if len(d)<10 or ng<10 or d.year.nunique()<3: return dict(n_records=len(d),n_trials=ng,status='too_few_for_trend')
    X=np.column_stack([pd.get_dummies(d.family,dtype=float).to_numpy(),(d.year.to_numpy(float)-2012)/10]); y=d[yc].to_numpy(float); v=d[vc].to_numpy(float); rank=np.linalg.matrix_rank(X); n=len(y)
    def est(xx,yy,vv):
        w=1/vv; ai=np.linalg.pinv((xx.T*w)@xx); b=ai@((xx.T*w)@yy); e=yy-xx@b
        ptr=w.sum()-np.trace(ai@((xx.T*w*w)@xx)); tau=max(0,((w*e*e).sum()-(len(yy)-np.linalg.matrix_rank(xx)))/ptr) if ptr>0 else 0
        w=1/(vv+tau); ai=np.linalg.pinv((xx.T*w)@xx); b=ai@((xx.T*w)@yy)
        return b,tau,ai,w
    b,tau,ai,w=est(X,y,v); ids,inv=np.unique(d.nct.to_numpy(),return_inverse=True)
    us=np.zeros((ng,X.shape[1])); np.add.at(us,inv,X*(w*(y-X@b))[:,None]); V=ai@(us.T@us)@ai*ng/(ng-1)*(n-1)/max(n-rank,1)
    se=math.sqrt(max(0,V[-1,-1])); beta=float(b[-1]); r=dict(n_records=n,n_trials=ng,n_families=int(d.family.nunique()),year_min=int(d.year.min()),year_max=int(d.year.max()),beta_per_decade=beta,cluster_se=se,ci95_cluster=[beta-1.96*se,beta+1.96*se],tau2=float(tau),status='estimated')
    if B:
        rng=np.random.default_rng(seed); pos=[np.flatnonzero(inv==i) for i in range(ng)]; vals=[]
        for _ in range(B):
            ii=np.concatenate([pos[i] for i in rng.integers(0,ng,size=ng)]); vals.append(est(X[ii],y[ii],v[ii])[0][-1])
        r.update(ci95_bootstrap=np.quantile(vals,[.025,.975]).tolist(),bootstrap_B=B)
    return r

def default(x):
    if isinstance(x,np.integer): return int(x)
    if isinstance(x,np.floating): return float(x)
    if isinstance(x,np.bool_): return bool(x)
    raise TypeError(type(x).__name__)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--raw',required=True); ap.add_argument('--records',required=True); ap.add_argument('--out',required=True); ap.add_argument('--bootstrap',type=int,default=500); a=ap.parse_args()
    raw=Path(a.raw); out=Path(a.out); out.mkdir(parents=True,exist_ok=True); sha=hashlib.sha256(raw.read_bytes()).hexdigest()
    if sha!=EXPECTED_SHA: raise RuntimeError('Raw snapshot SHA256 mismatch')
    orig=pd.read_parquet(a.records); frame,pairs=decorate(raw,orig); print('provenance reconstructed',len(frame),flush=True)
    frame.to_csv(out/'record_provenance.csv',index=False); pairs.to_csv(out/'raw_paired_contrasts.csv',index=False)
    pure=frame[~frame.combo].copy(); pure['y']=pure.resp/pure.sd; pure['v']=1/pure.n+pure.y**2/(2*pure.n); rr=pure[pure.raw_eligible].copy()
    r=dict(snapshot_sha256=sha,snapshot_bytes=raw.stat().st_size,record_table_sha256=hashlib.sha256(Path(a.records).read_bytes()).hexdigest(),analysis_status='post_hoc_after_inspection_of_original_results',legacy_all=dict(records=len(frame),trials=frame.nct.nunique()),legacy_pure=dict(records=len(pure),trials=pure.nct.nunique()),raw_pure=dict(records=len(rr),trials=rr.nct.nunique()),param_dispersion_counts=pure.groupby(['param_type','dispersion_type']).size().reset_index(name='n').to_dict('records'),raw_exclusions=pure.raw_exclusion_reason.value_counts().to_dict(),trends={},paired_trends={},pvalue_audit_counts=pure[pure.omtype=='PRIMARY'].p_audit_status.value_counts().to_dict())
    for name,cats in {'subjective':['CR','PR'],'laboratory':['OB'],'blood_pressure_weight':['OBB']}.items():
        for lab,dd,yc,vc in [('legacy_mixed_reporting',pure,'y','v'),('raw_mean_sd',rr,'raw_y','raw_v'),('raw_primary_outcomes',rr[rr.omtype=='PRIMARY'],'raw_y','raw_v'),('raw_unique_outcome',rr[rr.same_priority_outcomes==1],'raw_y','raw_v'),('raw_2006_2020',rr[rr.year.between(2006,2020)],'raw_y','raw_v')]:
            r['trends'][name+'__'+lab]=fit(dd[dd.cat.isin(cats)],yc,vc,a.bootstrap if lab in ('legacy_mixed_reporting','raw_mean_sd') else 0)
        for lab in ('experimental_only','all_nonplacebo'):
            pp=pairs[(pairs.pairing==lab)&pairs.pure_placebo&pairs.cat.isin(cats)]; r['paired_trends'][name+'__'+lab]=fit(pp,B=a.bootstrap)
    r['family_raw_trends']={f:fit(d,'raw_y','raw_v') for f,d in rr.groupby('family')}
    r['raw_by_category']=rr.groupby('cat').agg(records=('nct','size'),trials=('nct','nunique')).reset_index().to_dict('records')
    pure['period']=pd.cut(pure.year,[1999,2009,2015,2023],labels=['2000-2009','2010-2015','2016-2023']); mix=pure.groupby(['period','cat','param_type','dispersion_type'],observed=True).size().reset_index(name='records'); mix.to_csv(out/'reporting_mix_by_period.csv',index=False)
    r['reporting_mix_by_period']=mix.to_dict('records'); (out/'results.json').write_text(json.dumps(r,indent=2,ensure_ascii=False,default=default),encoding='utf8')
    print(json.dumps({k:r[k] for k in ('legacy_pure','raw_pure','raw_by_category','trends','paired_trends')},indent=2,default=default),flush=True)

if __name__=='__main__': main()
