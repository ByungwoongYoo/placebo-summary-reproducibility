"""Build readable review tables directly from saved numerical results."""
from pathlib import Path
import json,csv,hashlib
BASE=Path(__file__).resolve().parent
r=json.loads((BASE/'results/results.json').read_text(encoding='utf8'))
e=json.loads((BASE/'results/extended_results.json').read_text(encoding='utf8'))
s=json.loads((BASE/'results/source_adjudicated_sensitivity.json').read_text(encoding='utf8'))
manifest=json.loads((BASE/'validation/source_adjudication_manifest.json').read_text(encoding='utf8'))
adjudication_path=BASE/'validation/source_adjudication_12.csv'
adjudication_sha=hashlib.sha256(adjudication_path.read_bytes()).hexdigest()
assert adjudication_sha==manifest['adjudication_sha256']==s['adjudication_sha256']
with adjudication_path.open(newline='',encoding='utf8') as f:
    adjudication=list(csv.DictReader(f))
assert len(adjudication)==12 and len({d['nct'] for d in adjudication})==12
assert [(d['nct'],d['family'],int(d['outcome_index'])) for d in adjudication if d['confirmed_ineligible_for_arithmetic_sensitivity'].lower()=='true']==[('NCT01145755','MADRS',0)]
def cell(x): return str(x).replace('|','/')
def table(title,heads,rows,notes):
    return title+'\n\n'+' | '.join(heads)+'\n'+' | '.join(['---']*len(heads))+'\n'+'\n'.join(' | '.join(cell(x) for x in row) for row in rows)+'\n\n'+notes+'\n\n'
def trendrow(lab,d):
    ci=d.get('ci95_bootstrap',d.get('ci95_cluster'))
    return [lab,d['n_records'],d['n_trials'],f"{d['beta_per_decade']:.3f}",f"{ci[0]:.3f} to {ci[1]:.3f}",f"{d['bootstrap_B']}-replicate trial bootstrap" if 'ci95_bootstrap' in d else 'Trial-cluster normal interval']
labels={'subjective':'Subjective','laboratory':'Biochemical/spirometric','blood_pressure_weight':'Blood pressure/weight','legacy_mixed_reporting':'Mixed-summary reference','raw_mean_sd':'Reported mean/SD','all':'All eligible pharmacological contrasts','primary_outcomes':'Primary outcomes','single_active_arm':'One usable experimental summary','unique_outcome':'One eligible outcome per family','2006_2020':'Trials starting in 2006–2020','phase3':'Phase 3 trials','abs_g_at_most_3':'Excluding contrasts with absolute g >3','strict_fields_single_active':'Strict fields and one usable experimental summary','strict_fields_single_active_abs_g_at_most_3':'Strict fields, one usable summary and absolute g ≤3','single_experimental_group':'One identified experimental group','strict_fields_single_experimental_group':'Strict fields and one experimental group','strict_fields_single_experimental_group_abs_g_at_most_3':'Strict fields, one experimental group and absolute g ≤3'}
txt=table('Table 1. Summary types and record retention',['Reference-frame summary type','Dispersion label','Records'],[[{'MEAN':'Mean','LEAST_SQUARES_MEAN':'Least-squares mean'}[x['param_type']],x['dispersion_type'],x['n']] for x in r['param_dispersion_counts']], 'Denominator: 3,022 outcome-family records from 2,084 trials. The reported mean/SD and mapping rules retained 731 records from 538 trials; strict field conditions retained 636 records from 465 trials. Restrictions overlap; exclusions are listed jointly in Supplementary Table S1. A trial may contribute more than one family.')
rows=[]
for cat in ['subjective','laboratory','blood_pressure_weight']:
    for lab in ['legacy_mixed_reporting','raw_mean_sd']:
        rows.append(trendrow(labels[cat]+' / '+labels[lab],r['trends'][cat+'__'+lab]))
    rows.append(trendrow(labels[cat]+' / Strict fields',e['strict_field_placebo_trends'][cat]))
txt+=table('Table 2. Within-placebo temporal slopes per decade',['Analysis','Records','Trials','Slope','95% interval','Interval method'],rows,'Outcome-family-adjusted meta-regressions. Slopes are standardized change units per decade. Laboratory includes biochemical and spirometric outcomes. Cohort composition changes across restrictions, so differences between rows are not same-sample correction effects.')
duplicate_keys=['single_active_arm','strict_fields_single_active','strict_fields_single_active_abs_g_at_most_3']
rows=[trendrow(labels[k],v) for k,v in e['subjective_sensitivities'].items() if k not in duplicate_keys]
joint=s['results']['strict_single_experimental_group_abs_g_at_most_3']['source_informed']
rows.append(['Joint restriction plus confirmed adjusted-mean exclusion',joint['n_records'],joint['n_trials'],f"{joint['beta_per_decade']:.3f}",f"{joint['ci95_cluster'][0]:.3f} to {joint['ci95_cluster'][1]:.3f}",'Trial-cluster normal interval'])
txt+=table('Table 3. Subjective pharmacological contrast slopes per decade',['Analysis','Records','Trials','Hedges g slope','95% interval','Interval method'],rows,'Positive Hedges g favors the experimental group. Negative slopes indicate a decrease in the signed drug–placebo contrast in later-starting trials, not necessarily a decrease in its absolute magnitude. One experimental group requires one usable summary and only one identified EXPERIMENTAL group. Strict fields require one class, an outcome-level participant denominator, and no flagged model-based description. These post hoc comparisons are not confirmatory tests. The last row additionally excludes only NCT01145755/MADRS, identified as an adjusted mean in its final report. Other unresolved sampled records are retained; Supplementary Table S7 provides the corresponding bootstrap interval and other source-informed analyses.')
(BASE/'tables_placebo.txt').write_text(txt,encoding='utf8')
sup=table('Supplementary Table S1. Joint exclusion categories',['Reasons','Records'],[[k or 'Retained',v] for k,v in r['raw_exclusions'].items()],'Counts sum to 3,022. A record may meet several exclusion conditions; each combination is counted once.')
sup+=table('Supplementary Table S2. Reported P-value eligibility audit',['Status','Primary-outcome records'],list(r['pvalue_audit_counts'].items()),'This is an extraction-eligibility audit, not a determination of the validity or success of each original trial.')
sup+=table('Supplementary Table S3. Leave-one-family-out subjective contrasts',['Omitted family','Records','Trials','Slope','95% interval','Interval method'],[trendrow(k,v) for k,v in e['leave_one_family_out'].items()],'Same pharmacological set as Table 3 all; families are omitted individually. The joint strict analysis remains the decisive sensitivity when several uncertainties are addressed together.')
sup+=table('Supplementary Table S4. Same-cohort standardized arm changes',['Arm','Records','Trials','Slope','95% interval','Interval method'],[trendrow(k,v) for k,v in e['same_paired_cohort_arm_changes'].items()],'Each arm uses its own SD and random-effects weights. These slopes do not algebraically decompose the Hedges g contrast slope and do not identify a causal mediation effect.')
sup+=table('Supplementary Table S5. Usable-arm sensitivity definitions',['Analysis','Records','Trials','Slope','95% interval','Interval method'],[trendrow(labels[k],e['subjective_sensitivities'][k]) for k in duplicate_keys],'These three definitions count only usable experimental summaries. The additional group-count condition in main Table 3 confirms that no other identified experimental group was omitted. Both definitions yield identical records and estimates in this dataset.')
status_labels={
 'NUMERICAL_SUMMARY_UNRESOLVED':'Final numerical summary unresolved',
 'FINAL_NUMERICAL_SOURCE_UNAVAILABLE':'Final numerical source unavailable',
 'NUMERIC_FIELDS_CORROBORATED_PARTIAL':'Partial numerical corroboration',
 'NUMERICAL_POPULATION_DISCORDANCE_UNRESOLVED':'Numerical/population discrepancy unresolved',
 'CONFIRMED_ADJUSTED_MEAN_INELIGIBLE':'Confirmed adjusted mean, incompatible with arithmetic standardization'}
# Concise display summaries of the frozen, fully documented source_evidence field.
source_summaries={
 'NCT00711880':'Linked 2007 trial abstract reports adjusted means; raw values not recovered. A retrieved 2009 poster was a different trial.',
 'NCT00962104':'PMDA PDF p25 matches means and SDs; endpoint-specific denominators not separately recovered.',
 'NCT02913521':'Synopsis links DICL001/P130021 and describes week-8 change; final values and sign convention not established.',
 'NCT00485472':'SP905 final synopsis reports LS means different from registry values; no raw mean/SD table.',
 'NCT00401830':'No independent final numerical table recovered.',
 'NCT00109772':'Linked lenalidomide report identified, but selected numerical table not recovered.',
 'NCT02384460':'Article Table 2 matches means and n; rounded SEMs compatible with registry SDs.',
 'NCT04469270':'SAP specifies ITT analysis; selected final population and values not corroborated.',
 'NCT05624528':'Protocol specifies week-8 MADRS and a mixed model; no final numerical table.',
 'NCT00597896':'Article and registry differ in mean change and population information.',
 'NCT01145755':'Final CSR identifies LS means matching rounded registry values.',
 'NCT00154063':'No independent final numerical table recovered.'}
assert set(source_summaries)=={d['nct'] for d in adjudication}
source_rows=[[d['nct']+' / '+d['family'],source_summaries[d['nct']],status_labels[d['source_comparison_status']],
              'Exclude this outcome record only' if d['confirmed_ineligible_for_arithmetic_sensitivity'].lower()=='true' else 'Retain'] for d in adjudication]
sup+=table('Supplementary Table S6. Source comparison in all 12 sampled trials',['Trial / family','Source finding','Status for the selected record','Sensitivity action'],source_rows,'The denominator remains 12 sampled trials. Corroboration is field-specific and does not certify complete record validity. Unavailable or unresolved cases were not replaced, recoded or treated as correct. Complete source URLs, exact selected registry values, document locators and remaining uncertainty are in source_adjudication_12.csv. The two earlier |g|>3 cases were selected separately and are not included in this denominator.')
source_labels={
 'subjective_pharmacological_all':'Subjective pharmacological contrasts, broader eligible set',
 'strict_single_experimental_group':'Strict fields and one experimental group',
 'strict_single_experimental_group_abs_g_at_most_3':'Strict fields, one experimental group and absolute g ≤3',
 'subjective_placebo_strict_fields':'Subjective within-placebo change, strict fields'}
source_trends=[]
for key,label in source_labels.items():
    d=s['results'][key]['source_informed']
    source_trends.append([label,d['n_records'],d['n_trials'],f"{d['beta_per_decade']:.3f}",f"{d['ci95_cluster'][0]:.3f} to {d['ci95_cluster'][1]:.3f}",f"{d['ci95_bootstrap'][0]:.3f} to {d['ci95_bootstrap'][1]:.3f}"])
sup+=table('Supplementary Table S7. Sensitivity after excluding one confirmed adjusted-mean record',['Analysis','Records','Trials','Slope per decade','95% trial-cluster interval','95% trial-bootstrap interval'],source_trends,'All four rows exclude NCT01145755/MADRS/outcome index 0 only; no numeric summary was replaced. The first three slopes are signed Hedges g per decade, with positive g favoring the experimental group; the fourth is standardized within-placebo change per decade. Models otherwise match the earlier analyses. Trial bootstrap intervals use 1,000 replicates. Unresolved source comparisons remain included. These are post hoc sensitivity estimates, with pointwise intervals and no corpus-wide error-rate inference.')
(BASE/'supplementary_tables_placebo.txt').write_text(sup,encoding='utf8')
with (BASE/'results/table3_subjective_sensitivities.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.writer(f);w.writerow(['analysis','records','trials','slope_per_decade','ci95','ci_method']);w.writerows(rows)
print('Three main tables, seven supplementary tables and Table 3 CSV generated from reference, extended and source-informed results.')
