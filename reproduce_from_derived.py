"""Recompute the manuscript's models from published numerical derivatives.
This does not repeat extraction from the live registry or claim that source text
was externally validated. The extraction code is supplied separately.
"""
from pathlib import Path
import argparse,json
from types import SimpleNamespace
import numpy as np,pandas as pd
from reanalyse_placebo import fit,default,parse_p,pick_class

def main():
    base=Path(__file__).resolve().parent
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=base/'results');args=ap.parse_args()
    out=args.out;out.mkdir(parents=True,exist_ok=True)
    a=SimpleNamespace(bootstrap=500)
    source_meta=json.loads((base/'source_identity.json').read_text(encoding='utf-8'))
    frame=pd.read_csv(base/'inputs/record_numeric.csv')
    frame['raw_exclusion_reason']=frame.raw_exclusion_reason.fillna('')
    pairs=pd.read_csv(base/'inputs/contrast_numeric.csv')
    # The first operation uses only the base contrast fields; enriched fields
    # remain available for the extended checks below.
    pure=frame[~frame.combo].copy(); pure['y']=pure.resp/pure.sd; pure['v']=1/pure.n+pure.y**2/(2*pure.n); rr=pure[pure.raw_eligible].copy()
    r=dict(snapshot_sha256=source_meta['snapshot_sha256'],snapshot_bytes=source_meta['snapshot_bytes'],record_table_sha256=source_meta['record_table_sha256'],analysis_status='post_hoc_after_inspection_of_original_results',legacy_all=dict(records=len(frame),trials=frame.nct.nunique()),legacy_pure=dict(records=len(pure),trials=pure.nct.nunique()),raw_pure=dict(records=len(rr),trials=rr.nct.nunique()),param_dispersion_counts=pure.groupby(['param_type','dispersion_type']).size().reset_index(name='n').to_dict('records'),raw_exclusions=pure.raw_exclusion_reason.value_counts().to_dict(),trends={},paired_trends={},pvalue_audit_counts=pure[pure.omtype=='PRIMARY'].p_audit_status.value_counts().to_dict())
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

    rec=frame.copy();pp=pairs.copy()
    rec.to_csv(out/'record_provenance.csv',index=False)
    pp.to_csv(out/'raw_paired_contrasts_enriched.csv',index=False)
    fields=rec[['nct','family','outcome_index','strict_field_summary_frozen']].rename(columns={'strict_field_summary_frozen':'strict_field_summary'})
    strict=pp[(pp.pairing=='experimental_only')&pp.pure_placebo&pp.pharmacological].copy();s=strict[strict.cat.isin(['CR','PR'])].copy()
    r={'pharmacological_by_category':strict.groupby('cat').agg(records=('nct','size'),trials=('nct','nunique')).reset_index().to_dict('records'),'subjective_sensitivities':{},'leave_one_family_out':{},'within_intervention':{}}
    settings={'all':s,'primary_outcomes':s[s.omtype=='PRIMARY'],'single_active_arm':s[s.active_count==1],'unique_outcome':s[s.same_priority_outcomes==1],'2006_2020':s[s.year.between(2006,2020)],'phase3':s[s.phase=='P3'],'abs_g_at_most_3':s[s.y.abs()<=3],'strict_fields_single_active':s[s.strict_field_summary&(s.active_count==1)],'strict_fields_single_active_abs_g_at_most_3':s[s.strict_field_summary&(s.active_count==1)&(s.y.abs()<=3)]}
    for lab,d in settings.items():r['subjective_sensitivities'][lab]=fit(d,B=1000 if lab in ('all','strict_fields_single_active') else 0)
    genuine_single=s[(s.active_count==1)&(s.experimental_group_count==1)]
    r['subjective_sensitivities']['single_experimental_group']=fit(genuine_single)
    r['subjective_sensitivities']['strict_fields_single_experimental_group']=fit(genuine_single[genuine_single.strict_field_summary],B=1000)
    r['subjective_sensitivities']['strict_fields_single_experimental_group_abs_g_at_most_3']=fit(genuine_single[genuine_single.strict_field_summary&(genuine_single.y.abs()<=3)])
    for fam in sorted(s.family.unique()):r['leave_one_family_out'][fam]=fit(s[s.family!=fam])
    r['pharmacological_trends']={lab:fit(strict[strict.cat.isin(cats)],B=1000) for lab,cats in {'subjective':['CR','PR'],'laboratory':['OB'],'blood_pressure_weight':['OBB']}.items()}
    w=s.copy();w['stratum']=w.intervention_key+' | '+w.family
    gg=w.groupby('stratum').agg(trials=('nct','nunique'),span=('year',lambda z:z.max()-z.min()));keep=gg[(gg.trials>=2)&(gg.span>=3)].index;w=w[w.stratum.isin(keep)].copy()
    r['within_intervention'].update(eligible_strata=len(keep),eligible_trials=w.nct.nunique(),eligible_records=len(w))
    if len(w):w['family']=w.stratum;r['within_intervention']['trend']=fit(w,B=1000);w.to_csv(out/'within_exact_intervention.csv',index=False)
    rr=rec[rec.raw_eligible&~rec.combo].merge(fields[['nct','family','outcome_index','strict_field_summary']],on=['nct','family','outcome_index'],validate='one_to_one');rr=rr[rr.strict_field_summary]
    r['strict_field_placebo_trends']={lab:fit(rr[rr.cat.isin(cats)],'raw_y','raw_v',B=1000) for lab,cats in {'subjective':['CR','PR'],'laboratory':['OB'],'blood_pressure_weight':['OBB']}.items()}
    r['strict_field_retention']={'records':len(rr),'trials':rr.nct.nunique()}
    r['same_paired_cohort_arm_changes']={}
    for lab,mcol,sdcol,ncol in [('placebo','mean_placebo','sd_placebo','n_placebo'),('active','mean_active','sd_active','n_active')]:
        x=s.copy();x['arm_y']=x.improvement_sign*x[mcol]/x[sdcol];x['arm_v']=1/x[ncol]+x.arm_y**2/(2*x[ncol]);r['same_paired_cohort_arm_changes'][lab]=fit(x,'arm_y','arm_v',B=1000)
    expected={'<0.05':True,'<=0.05':None,'>=0.01':None,'>0.05':False,'0.05':False,'0.049':True,'<0.1':None,'≥0.05':False}
    assert all(parse_p(k)[0] is v for k,v in expected.items())
    assert pick_class({'classes':[{'title':'Week 12'},{'title':'Week 4'}]})[0]==0
    assert pick_class({'classes':[{'title':'Male'},{'title':'Female'}]})[1] is None
    raw=rec[rec.raw_eligible];assert raw.param_type.eq('MEAN').all() and raw.dispersion_type.eq('Standard Deviation').all()
    assert (raw.extracted_spread==raw.raw_sd).all();assert not raw.n_source.eq('time_specific_n_unverified').any()
    assert len(raw)==len(raw.drop_duplicates(['nct','family']))
    df=pp.n_active+pp.n_placebo-2;sp=np.sqrt(((pp.n_active-1)*pp.sd_active**2+(pp.n_placebo-1)*pp.sd_placebo**2)/df)
    independently_recalculated=(1-3/(4*df-1))*pp.improvement_sign*(pp.mean_active-pp.mean_placebo)/sp
    assert np.allclose(pp.y,independently_recalculated,rtol=1e-12,atol=1e-12)
    r['boundary_checks']={'pvalue_interval_cases':len(expected),'nonchronological_and_subgroup_class_cases':2,'raw_summary_type_and_n_invariants':'PASS','one_record_per_trial_family':'PASS','all_paired_g_recalculated':'PASS'}
    (out/'extended_results.json').write_text(json.dumps(r,indent=2,ensure_ascii=False,default=default),encoding='utf8');print(json.dumps(r,indent=2,default=default),flush=True)

if __name__=='__main__':main()
