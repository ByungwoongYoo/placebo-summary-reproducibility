# Summary-statistic meaning in placebo-controlled trial analyses

This companion contains code and numerical derivatives for the study of how
summary type, denominator, time point and arm identity affect temporal analyses
of ClinicalTrials.gov results. The analysis is post hoc. The strict joint
sensitivity retains 73 records from 64 trials and its interval includes zero.

## Reproduce the reported calculations

Use Python 3.12 and the packages in `requirements.txt`. From this directory run:

```sh
python reproduce_from_derived.py
python source_adjudicated_sensitivity.py
python build_tables.py
```

The first command recomputes the reference and extended models from
`inputs/record_numeric.csv` and `inputs/contrast_numeric.csv`. It uses the same
model, restrictions, trial resampling and seeds as the full extraction pipeline.
The next commands recompute the source-informed sensitivity and all three main
and seven supplementary tables. The `expected` directory contains the frozen
outputs for comparison. Source-input hashes differ from a full raw extraction
because narrative fields have been removed; numeric model outputs are unchanged.

`reanalyse_placebo.py` and `families.py` also document the source extraction and
matching operations. The compact numeric route starts after extraction; it does
not independently validate the registry's labels or selected summaries. Trial
IDs, outcome indices and group identifiers permit inspection of the public
registry records. The source snapshot and record-table hashes are in
`source_identity.json`.

## Source comparison

The validation files retain all twelve sampled trials, their selected numeric
summaries, public source URLs, document locators and unresolved findings. Only
NCT01145755/MADRS/outcome 0 is excluded in the source-informed sensitivity.
The release omits trial-description paragraphs, local evidence filenames and
reviewer placeholders. It preserves the sampling rows, numeric values and
scientific classifications. The manifest records both original and reduced-file
hashes. No source documents or publisher PDFs are redistributed.

## Data source and reuse

Source: [ClinicalTrials.gov](https://clinicaltrials.gov/), frozen processed snapshot
dated 23 October 2026. This is a historical research dataset, not a current mirror
of the database. Current study pages remain available through the recorded NCT
identifiers. The derivatives select outcome-family records, retain structured
numeric summaries, orient changes, standardize effects and add eligibility and
source-comparison fields. Narrative descriptions and contact information are
excluded from the numerical inputs.

The source's [Terms and Conditions](https://clinicaltrials.gov/about-site/terms-conditions)
govern ClinicalTrials.gov material, including source attribution, processing
dates, modification disclosure and applicable third-party rights. No proprietary
claim is made over the database or its source material. This companion does not
grant a new licence to third-party trial reports. It contains no individual
participant data, clinical-study reports or copyrighted article PDFs.

`CITATION.cff` identifies the author and the version of this companion.
The analysis code is available under the MIT licence; `REUSE.md` describes the
separate conditions for source-derived data and author-generated outputs.
