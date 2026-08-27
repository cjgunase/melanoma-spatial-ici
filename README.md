# Melanoma Spatial ICI

## Single-cell-to-spatial analysis of anti-PD-1 resistance

This project tests whether spatially organized melanoma, T-cell, and myeloid
states provide a more informative view of immune-checkpoint-inhibitor response
than cell abundance or bulk expression alone.

The analysis integrates three public melanoma cohorts:

- **GSE115978:** scRNA-seq reference atlas of melanoma ecosystems.
- **GSE300445:** Visium spatial transcriptomics from four pretreatment primary
  melanomas.
- **GSE233305:** spatial whole-transcriptome profiles from tumor, leukocyte,
  and macrophage compartments in immunotherapy-treated melanoma.

### Planned analytical aims

1. Build a patient-aware single-cell reference of malignant, immune, and
   stromal states.
2. Map those states into Visium tissue sections using complementary reference
   mapping approaches.
3. Identify recurrent tumor-immune spatial neighborhoods and candidate
   cell-cell interactions.
4. Test a concise spatial resistance hypothesis in an independent clinical
   cohort.

### Reproducibility principles

- Raw public data are downloaded from GEO and excluded from Git.
- Every input is recorded with its URL, expected size, and SHA-256 checksum.
- Statistical inference uses the patient, not the individual cell, as the
  experimental unit whenever patient-level replication is available.
- Exploratory associations are not presented as clinically validated
  predictors.
- Generated figures and reports must be reproducible from versioned code.

### Quick start

```bash
uv sync --python 3.11 --extra dev
uv run python scripts/download_data.py
uv run python scripts/inventory_data.py --strict
uv run python scripts/01_ingest_scrna.py
uv run python scripts/02_build_scrna_atlas.py
uv run python scripts/03_map_visium.py
uv run python scripts/04_spatial_neighborhoods.py
uv run python scripts/05_validate_geomx.py
uv run python scripts/06_build_report.py
uv run pytest
```

The complete workflow is also available as `bash scripts/run_all.sh`.
The latest reproducible result is summarized in
[`docs/analysis_report.md`](docs/analysis_report.md).

The computational environment and first analysis entry points will be added as
the corresponding milestones are verified.

### Data policy

This repository contains code and small provenance manifests only. GEO source
files and derived matrices are deliberately excluded from version control.

### Status

**Milestones 1–7 implemented:** the workflow validates all 12 inputs; builds a
patient-aware GSE115978 reference; maps that reference into four GSE300445
Visium sections with complementary approaches; identifies recurrent spatial
neighborhoods with sensitivity checks; and tests one predeclared
dedifferentiated-melanoma resistance score in patient-collapsed GSE233305 ITX2
data. GSE115978 ingestion recovered 7,186 cells, 23,686 genes, and 32 patients;
7,072 cells (98.4%) pass baseline QC. Mitochondrial percentage is explicitly
unavailable because the source matrix contains no mitochondrial genes.

Analysis outputs remain excluded from Git and are regenerated from the commands
above. GSE300445 has only four spatial samples and no response labels in the GEO
record, so its associations remain exploratory. The independent retrospective
validation is not presented as a clinically validated predictor.

### Future directions: Melanoma Spatial Resistance Atlas

The current workflow is a foundation, not a finished clinical predictor. Its
primary six-gene tumor-compartment hypothesis was not supported in independent
ITX2 validation (35 patients; Mann–Whitney p = 0.778; nonresponse AUC = 0.470),
and the four Visium sections showed heterogeneous tumor–immune relationships.
Future work must preserve this null result and must not tune a replacement
signature on the ITX2 primary-test labels.

The next objective is a patient-level **spatial ecosystem fingerprint** that
combines tumor state, immune state, tissue geometry, and mapping uncertainty.
The following work packages are ordered by dependency.

#### F1. Reproduce the published GSE233305 compartment signatures

**Question:** Can this repository reproduce the published CD45, CD68, S100B,
and pseudo-bulk response models before proposing a new model?

**Implementation:**

1. Encode the published training/validation split, response definition, gene
   lists, coefficient signs, preprocessing, and decision thresholds in a
   versioned configuration file. The published S100B model contains positive
   response coefficients for `PSMB8`, `TAX1BP3`, `NOTCH3`, `LCP2`, and `NQO1`,
   and resistance coefficients for `KMT2C`, `OVCA2`, and `MGRN1`.
2. Align GeoMx segments to patients and compartments with explicit assertions;
   collapse repeated regions to patients before inference.
3. Reproduce discovery and validation ROC curves, confusion matrices,
   calibration plots, and bootstrap 95% confidence intervals.
4. Compare reported and reproduced sample counts, exclusions, AUCs, and
   coefficients in `results/07_published_signature_reproduction/`.

**Completion gate:** no new signature development until cohort membership and
published performance are reproduced or every discrepancy is documented. AUC
agreement alone is insufficient if patients or labels differ.

#### F2. Benchmark single-cell-to-spatial mapping

**Question:** Are inferred cell states stable across established mapping or
deconvolution methods?

**Methods:** retain the current projected-reference and marker-score baselines,
then add at least two established approaches:

- [cell2location](https://doi.org/10.1038/s41587-021-01139-4) for Bayesian
  abundance estimation with platform and sensitivity effects;
- [RCTD](https://doi.org/10.1038/s41587-021-00830-w) for robust mixture
  decomposition and cross-platform correction;
- optionally [Tangram](https://doi.org/10.1038/s41592-021-01264-7) for
  probabilistic single-cell-to-space alignment.

**Benchmark design:**

1. Freeze one harmonized gene universe and one patient-balanced GSE115978
   reference before running any method.
2. Hold out 10%–20% of shared genes from mapping and measure their spatial
   reconstruction by per-gene Pearson/Spearman correlation and RMSE.
3. Create pseudo-spots with known cell-type mixtures from held-out patients;
   report abundance correlation, mean absolute error, calibration, and rare
   cell-state recall.
4. Repeat reference construction across patient bootstraps and downsampling
   depths to quantify stability.
5. Report per-spot method agreement, entropy/uncertainty, runtime, peak memory,
   and failure rate. Never use consensus as proof of biological correctness.
6. Overlay inferred S100B/SOX10 melanoma, CD3 T-cell, and CD68 macrophage
   abundance on tissue images for blinded pathology review when an expert is
   available.

**Completion gate:** a cell state may enter the ecosystem model only if it is
recoverable in pseudo-spots, stable to patient/reference resampling, and
supported by at least two mapping views or by orthogonal pathology evidence.

#### F3. Resolve biologically meaningful cell states

**Question:** Which tumor and immune programs should replace broad cell classes?

**Candidate states and programs:**

- melanoma differentiation/plasticity: melanocytic, transitory,
  neural-crest-like, and undifferentiated states;
- macrophages: inflammatory/C1QC-like, SPP1-like, interferon-responsive, and
  hypoxia-associated programs;
- T cells: cytotoxic, memory, exhausted/dysfunctional, regulatory, and cycling
  states;
- microenvironment: IFN-γ response, antigen presentation, hypoxia, ECM/CAF,
  angiogenesis, neutrophil activation, and tertiary-lymphoid/B-cell programs.

Derive scores from published gene sets before outcome analysis. Estimate state
profiles from patient-level pseudobulk or patient-aware mixed models; do not
treat individual cells as independent replicates. Record gene-set provenance,
version, direction, and minimum gene coverage in `config/programs.yaml`.

**Completion gate:** each retained program must have biological provenance,
adequate gene coverage in both reference and spatial assays, patient-level
variation, and robustness to leave-one-patient-out reference construction.

#### F4. Construct interpretable spatial features

**Question:** Which spatial relationships add information beyond abundance?

For each section, derive features at several prespecified spatial scales:

- cell-state abundance and compositional ratios;
- tumor-to-T-cell and tumor-to-macrophage nearest-neighbor distances;
- contact or adjacency enrichment relative to label-preserving spatial
  permutations;
- tumor–immune boundary density and infiltration depth;
- local mixing/segregation, neighborhood entropy, and spatial autocorrelation;
- recurrent neighborhood proportions and transition frequencies;
- ligand–receptor hypotheses restricted to genes expressed in the mapped
  sender and receiver states.

Run analyses at spot-neighborhood sizes such as k = 4, 6, and 10 and at
distance thresholds tied to Visium spot geometry. Retain a feature only when
its direction is stable across scales and it is not explained by section area,
library size, or broad cell abundance alone. Use within-section coordinate
permutations to preserve tissue shape and control spatial autocorrelation.

**Completion gate:** publish a data dictionary containing the formula, units,
scale, missingness rule, biological interpretation, and null model for every
feature.

#### F5. Define and freeze the ecosystem fingerprint

**Question:** Can a small, interpretable feature set summarize resistance
architecture without overfitting?

Use discovery data only to define a concise model. Candidate components may
include melanoma differentiation, macrophage state, T-cell state, boundary
exclusion, and IFN-γ/hypoxia context. Compare each spatial model with two
prespecified baselines: clinical covariates alone and non-spatial cell-state
abundance alone.

Because patient counts are small, prefer penalized logistic regression or a
similarly constrained model over high-capacity classifiers. Use nested,
patient-grouped cross-validation; keep all regions from one patient in the same
fold. Report optimism-corrected AUC, precision-recall AUC, Brier score,
calibration slope/intercept, confidence intervals, and decision thresholds.
Do not select features using ITX2 primary-test outcomes.

Before external validation, freeze:

- eligible patients and exclusion rules;
- response endpoint and time horizon;
- preprocessing and missing-data rules;
- exact feature formulas and coefficients;
- primary metric, statistical test, direction, and success threshold;
- multiplicity plan and sensitivity analyses.

Store this specification in a dated, immutable analysis-plan document and tag
the corresponding commit.

#### F6. Obtain response-annotated spatial validation data

GSE300445 exposes four Visium sections, while its GEO description reports a
broader study with 57 pretreatment bulk samples and outcome-associated spatial
ecosystems. Inspect the associated preprint supplements for an unambiguous
sample-to-response table. If it is not public, request only the minimum
de-identified mapping needed from the study authors; document the request and
do not infer outcomes from figures or sample ordering.

Additional validation may use public response-annotated melanoma spatial
transcriptomic, GeoMx, or spatial-proteomic cohorts, but assay-specific models
must not be treated as directly interchangeable. Define a cross-assay feature
mapping before examining outcomes and require a minimum number of patients in
both outcome groups. Small cohorts are suitable for effect-size estimation and
uncertainty, not definitive classifier claims.

**Completion gate:** run the frozen model once on the untouched validation
cohort. A failed primary test remains a reported result; it does not trigger
outcome-guided model revision on the same cohort.

#### F7. Add pathology review and spatial interpretability

Generate one standardized panel per tissue section containing:

1. the original histology image;
2. QC metrics and tissue mask;
3. mapped cell-state abundance with uncertainty;
4. tumor–immune boundaries and neighborhood assignments;
5. method-disagreement regions;
6. features contributing to the patient-level fingerprint.

Provide the reviewer with randomized section identifiers and a structured
scoring form for melanoma-rich, lymphocyte-rich, macrophage-rich, necrotic,
stromal, and artifact regions. Measure inter-rater agreement if multiple
reviewers are available. Pathology review validates spatial plausibility; it
must not be used to silently relabel samples after outcome analysis.

#### F8. Statistical and reproducibility safeguards

- Use patients as independent units; spots, cells, and regions are nested
  observations.
- Report effect sizes and confidence intervals alongside p-values.
- Control false discovery rate within each declared exploratory family.
- Use spatially constrained permutations for neighborhood tests.
- Run leave-one-patient-out, leave-one-section-out, mapping-method, QC-threshold,
  and spatial-scale sensitivity analyses.
- Record random seeds, package locks, hardware, runtime, peak memory, input
  checksums, and output provenance.
- Add synthetic fixtures for unit tests and one checksum-pinned miniature
  integration dataset for CI. Full GEO execution remains a separate workflow.
- Mark outputs as exploratory, confirmatory, supported, not supported, or not
  testable; never convert missing metadata into inferred labels.

#### F9. Proposed deliverables

**Near-term benchmark/reproducibility paper:**

- exact reproduction of published GSE233305 signatures;
- multi-method spatial-mapping benchmark;
- uncertainty and pathology-alignment analysis;
- reusable, tested workflow and documented null results.

**Later biological/translational paper, conditional on validation:**

- frozen multicomponent ecosystem fingerprint;
- response-annotated discovery and untouched external validation cohorts;
- evidence that spatial features improve on clinical and abundance-only
  baselines;
- interpretable tissue maps and pathology concordance;
- calibrated effect estimates without claims of clinical readiness unless
  prospective validation is performed.

#### F10. Recommended execution order

1. Reproduce GSE233305 published models (F1).
2. Benchmark mapping and quantify uncertainty (F2).
3. Resolve candidate biological states (F3).
4. Implement spatial features and permutation nulls (F4).
5. Write and freeze the model specification (F5).
6. Acquire or identify response-annotated validation data (F6).
7. Run pathology review without outcome labels (F7).
8. Execute the frozen external validation once (F5–F8).
9. Prepare the appropriate deliverable based on the evidence (F9).

### Primary sources

- [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
- [GSE300445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300445)
- [GSE233305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233305)
- [Published GSE233305 compartment signatures](https://doi.org/10.1158/1078-0432.CCR-23-3932)
- [GSE300445-associated preprint](https://doi.org/10.1101/2025.07.07.661465)
