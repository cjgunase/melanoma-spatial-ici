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
3. identify recurrent tumor-immune spatial neighborhoods and candidate
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
uv sync --python 3.12
uv run python scripts/download_data.py
uv run python scripts/inventory_data.py --strict
uv run python scripts/01_ingest_scrna.py
uv run pytest
```

The computational environment and first analysis entry points will be added as
the corresponding milestones are verified.

### Data policy

This repository contains code and small provenance manifests only. GEO source
files and derived matrices are deliberately excluded from version control.

### Status

**Milestone 1 complete:** all 12 declared source files were downloaded and
validated. GSE115978 ingestion recovered 7,186 cells, 23,686 genes, and 32
patients; 7,072 cells (98.4%) pass the documented baseline detected-gene and
library-size criteria. The source matrix does not contain mitochondrial genes,
so mitochondrial percentage is explicitly unavailable rather than imputed.

Analysis outputs remain excluded from Git and are regenerated from the commands
above. The next milestone is patient-aware normalization, dimensionality
reduction, and cell-state atlas construction.

### Primary sources

- [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
- [GSE300445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300445)
- [GSE233305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233305)
