# Melanoma spatial ICI analysis report

## Design

This analysis asks whether a dedifferentiated melanoma program occurs in recurrent
macrophage-rich/T-cell-poor spatial neighborhoods and whether its compartment-level
score is associated with ICI nonresponse. Spatial discovery is exploratory because
GSE300445 contains only four sections and no GEO response labels. Outcome testing is
reserved for the independent GSE233305 ITX2 cohort and collapses segments to patients.

## Data integrity and single-cell reference

- All 12 declared GEO files passed strict inventory checks.
- GSE115978: 7,186 cells, 23,686 genes, 32 patients.
- 7,072 cells passed baseline QC.
- The atlas uses 32 patient identifiers and patient-equal marker profiles.

## Spatial discovery

- 11,558 tissue spots across 4 sections.
- Five recurrent neighborhoods were summarized with sensitivity analyses at k=4, 6, and 10.
- Two mapping views are retained: nonnegative projected reference profiles and patient-aware marker scores.

## Independent validation

Primary hypothesis: Higher melanoma dedifferentiation score associates with nonresponse in S100B segments.

- ITX2 S100B patients analyzed: 35.
- Mann–Whitney p-value: 0.7784.
- AUC for nonresponse direction: 0.470.
- Result: **not supported** at an unadjusted two-sided 0.05 threshold.

This result is an external retrospective association, not a clinically validated predictor.
The cohort is not used to revise the signature after seeing outcomes.

## Reproduction

Run `bash scripts/run_all.sh`. Raw and generated files remain excluded from Git; code,
manifests, tests, and this report are versioned.
