# Analysis plan

## Primary question

Do spatially organized melanoma, T-cell, and myeloid states provide information
about anti-PD-1 resistance beyond cell-state abundance alone?

## Design commitments

- Single-cell discovery will retain patient identifiers throughout processing.
- Differential tests will use patient-level pseudobulk or patient-aware models.
- Spatial analyses will distinguish measured expression, inferred cell-state
  abundance, and hypothesized interactions.
- GSE233305 will be treated as an independent validation cohort, not as an
  additional feature-selection cohort.
- The four-sample Visium cohort is not large enough to train or claim a
  clinically validated response classifier.

## Milestones

1. Data integrity and metadata audit.
2. Reproducible scRNA-seq ingestion and QC report.
3. Patient-aware cell-state atlas.
4. Visium ingestion, QC, and reference mapping.
5. Spatial-neighborhood analysis and sensitivity checks.
6. Independent compartment-level validation.
7. Reproducible report and interview-ready graphical summary.

## Known source-data constraints

- GSE115978 is a processed SMART-seq2 count matrix with published cell labels.
  It contains no `MT-` genes, so mitochondrial-read percentage cannot be used as
  a QC criterion. Filtering is therefore based on detected genes and library
  size, with published read and gene counts retained for concordance checks.
