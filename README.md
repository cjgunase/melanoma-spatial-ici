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

### Primary sources

- [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
- [GSE300445](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300445)
- [GSE233305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233305)
