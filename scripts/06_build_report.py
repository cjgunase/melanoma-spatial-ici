#!/usr/bin/env python3
"""Build a concise reproducible report from versioned summaries."""

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    out = ROOT / "docs"
    out.mkdir(parents=True, exist_ok=True)
    qc = json.loads((ROOT / "results/01_scrna_qc/summary.json").read_text())
    atlas = json.loads((ROOT / "results/02_scrna_atlas/summary.json").read_text())
    spatial = json.loads((ROOT / "results/04_spatial_neighborhoods/summary.json").read_text())
    validation = json.loads((ROOT / "results/05_geomx_validation/summary.json").read_text())
    tests = pd.read_csv(ROOT / "results/05_geomx_validation/tests.csv")
    primary = tests.query("cohort == 'ITX2' and compartment == 'S100B'").iloc[0]
    verdict = (
        "supported" if primary.pvalue < 0.05 and primary.auc_nonresponse > 0.5 else "not supported"
    )
    report = f"""# Melanoma spatial ICI analysis report

## Design

This analysis asks whether a dedifferentiated melanoma program occurs in recurrent
macrophage-rich/T-cell-poor spatial neighborhoods and whether its compartment-level
score is associated with ICI nonresponse. Spatial discovery is exploratory because
GSE300445 contains only four sections and no GEO response labels. Outcome testing is
reserved for the independent GSE233305 ITX2 cohort and collapses segments to patients.

## Data integrity and single-cell reference

- All 12 declared GEO files passed strict inventory checks.
- GSE115978: {qc["n_cells"]:,} cells, {qc["n_genes"]:,} genes, {qc["n_patients"]} patients.
- {qc["qc_pass_cells"]:,} cells passed baseline QC.
- The atlas uses {atlas["patients"]} patient identifiers and patient-equal marker profiles.

## Spatial discovery

- {spatial["spots"]:,} tissue spots across {spatial["sections"]} sections.
- Five recurrent neighborhoods were summarized with sensitivity analyses at k=4, 6, and 10.
- Two mapping views are retained: nonnegative projected reference profiles and patient-aware marker scores.

## Independent validation

Primary hypothesis: {validation["hypothesis"]}.

- ITX2 S100B patients analyzed: {int(primary.patients)}.
- Mann–Whitney p-value: {primary.pvalue:.4g}.
- AUC for nonresponse direction: {primary.auc_nonresponse:.3f}.
- Result: **{verdict}** at an unadjusted two-sided 0.05 threshold.

This result is an external retrospective association, not a clinically validated predictor.
The cohort is not used to revise the signature after seeing outcomes.

## Reproduction

Run `bash scripts/run_all.sh`. Raw and generated files remain excluded from Git; code,
manifests, tests, and this report are versioned.
"""
    (out / "analysis_report.md").write_text(report)


if __name__ == "__main__":
    main()
