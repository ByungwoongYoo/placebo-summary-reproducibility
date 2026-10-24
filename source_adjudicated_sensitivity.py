"""Post hoc source-informed sensitivity, preserving all earlier results.

Exclude only the frozen, confirmed adjusted-mean record. Unresolved records
and registry numeric values are not modified. No network operations occur.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json
import numpy as np
import pandas as pd
from reanalyse_placebo import fit, default

B=Path(__file__).resolve().parent
R=B/'results'
V=B/'validation'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
manifest=json.loads((V/'source_adjudication_manifest.json').read_text(encoding='utf-8'))
assert sha(V/'source_adjudication_12.csv')==manifest['adjudication_sha256']
assert sha(B/'validation_sample_frozen.csv')==manifest['original_sample_sha256']
ad=pd.read_csv(V/'source_adjudication_12.csv')
assert len(ad)==12 and ad.nct.nunique()==12
exc=ad[ad.confirmed_ineligible_for_arithmetic_sensitivity]
keys=set(zip(exc.nct,exc.family,exc.outcome_index))
assert keys=={('NCT01145755','MADRS',0)}
protected=[R/'extended_results.json',R/'raw_paired_contrasts_enriched.csv',R/'record_provenance.csv']
original_hashes={p.name:sha(p) for p in protected}
pp=pd.read_csv(R/'raw_paired_contrasts_enriched.csv')
rec=pd.read_csv(R/'record_provenance.csv')
old=json.loads((R/'extended_results.json').read_text(encoding='utf-8'))
s=pp[(pp.pairing=='experimental_only')&pp.pure_placebo&pp.pharmacological&pp.cat.isin(['CR','PR'])].copy()
single=s[(s.active_count==1)&(s.experimental_group_count==1)&s.strict_field_summary].copy()
joint=single[single.y.abs()<=3].copy()
strict=rec.strict_field_summary_frozen.astype(bool)
placebo=rec[rec.raw_eligible&~rec.combo&strict&rec.cat.isin(['CR','PR'])].copy()
datasets={
 'subjective_pharmacological_all':(s,'y','v',old['subjective_sensitivities']['all']),
 'strict_single_experimental_group':(single,'y','v',old['subjective_sensitivities']['strict_fields_single_experimental_group']),
 'strict_single_experimental_group_abs_g_at_most_3':(joint,'y','v',old['subjective_sensitivities']['strict_fields_single_experimental_group_abs_g_at_most_3']),
 'subjective_placebo_strict_fields':(placebo,'raw_y','raw_v',old['strict_field_placebo_trends']['subjective'])}
results={}
retention=[]
for name,(d,yc,vc,expected) in datasets.items():
    before=fit(d,yc,vc)
    assert before['n_records']==expected['n_records'] and before['n_trials']==expected['n_trials']
    assert np.isclose(before['beta_per_decade'],expected['beta_per_decade'],atol=1e-12,rtol=0)
    assert np.allclose(before['ci95_cluster'],expected['ci95_cluster'],atol=1e-12,rtol=0)
    flagged=np.array([(n,f,int(o)) in keys for n,f,o in zip(d.nct,d.family,d.outcome_index)])
    assert flagged.sum()==1, (name,flagged.sum())
    after=d.loc[~flagged].copy()
    result=fit(after,yc,vc,B=1000)
    results[name]={'original_reproduced':before,'source_informed':result,'excluded_records':d.loc[flagged,['nct','family','outcome_index',yc,vc]].to_dict('records')}
    for row,ex in zip(d.to_dict('records'),flagged):
        retention.append({'analysis':name,'nct':row['nct'],'family':row['family'],'outcome_index':row['outcome_index'],'excluded_confirmed_adjusted_mean':bool(ex),'value':row[yc],'variance':row[vc]})
assert all(sha(p)==original_hashes[p.name] for p in protected)
out=dict(analysis_generated_at_utc=datetime.now(timezone.utc).isoformat(),post_hoc=True,adjudication_frozen_before_reestimation=manifest['frozen_at_utc'],adjudication_sha256=manifest['adjudication_sha256'],original_inputs_sha256=original_hashes,excluded_keys=manifest['confirmed_excluded_keys'],source_sample_denominator_trials=12,source_sample_not_a_corpus_error_rate=True,numeric_replacement_performed=False,bootstrap_seed=2026102401,results=results)
(R/'source_adjudicated_sensitivity.json').write_text(json.dumps(out,indent=2,default=default),encoding='utf-8')
pd.DataFrame(retention).to_csv(R/'source_adjudicated_retention.csv',index=False)
print(json.dumps(out,indent=2,default=default))
