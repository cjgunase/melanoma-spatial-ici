#!/usr/bin/env python3
"""Create a validated GSE115978 AnnData object and baseline QC report."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from melanoma_spatial_ici.io import read_gse115978
from melanoma_spatial_ici.qc import calculate_scrna_qc

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "GSE115978_scRNA"
PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results" / "01_scrna_qc"


def main() -> None:
    PROCESSED.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    adata = read_gse115978(
        RAW / "GSE115978_counts.csv.gz",
        RAW / "GSE115978_cell.annotations.csv.gz",
    )
    calculate_scrna_qc(adata)

    summary = {
        "accession": "GSE115978",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_patients": int(adata.obs["samples"].nunique()),
        "cell_types": adata.obs["cell.types"].value_counts().to_dict(),
        "treatment_groups": adata.obs["treatment.group"].value_counts().to_dict(),
        "qc_pass_cells": int(adata.obs["pass_qc"].sum()),
        "qc_pass_fraction": float(adata.obs["pass_qc"].mean()),
        "qc_thresholds": adata.uns["qc_thresholds"],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    columns = [
        "samples",
        "cell.types",
        "treatment.group",
        "Cohort",
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "pass_qc",
    ]
    adata.obs[columns].to_csv(RESULTS / "cell_qc.csv.gz")

    plot_data = adata.obs.reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    sns.histplot(plot_data, x="n_genes_by_counts", hue="pass_qc", bins=60, ax=axes[0])
    sns.histplot(plot_data, x="total_counts", hue="pass_qc", bins=60, ax=axes[1])
    axes[1].set_xscale("log")
    sns.scatterplot(
        plot_data,
        x="no.of.reads",
        y="total_counts",
        hue="pass_qc",
        alpha=0.35,
        s=12,
        linewidth=0,
        ax=axes[2],
    )
    axes[2].set_xscale("log")
    axes[2].set_yscale("log")
    axes[2].set_xlabel("Published read count")
    axes[2].set_ylabel("Count-matrix total")
    for axis in axes:
        axis.grid(False)
    fig.suptitle("GSE115978 baseline single-cell quality metrics", fontweight="bold")
    fig.tight_layout()
    fig.savefig(RESULTS / "qc_distributions.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    adata.write_h5ad(PROCESSED / "GSE115978_ingested_qc.h5ad", compression="gzip")
    print(pd.Series(summary).to_string())


if __name__ == "__main__":
    main()
