#!/usr/bin/env python3
"""Build a patient-aware cell-state reference atlas."""

from pathlib import Path

import matplotlib.pyplot as plt
import scanpy as sc

from melanoma_spatial_ici.analysis import (
    patient_aware_markers,
    patient_celltype_profiles,
    write_json,
)

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    source = ROOT / "data/processed/GSE115978_ingested_qc.h5ad"
    output = ROOT / "data/processed/GSE115978_reference_atlas.h5ad"
    results = ROOT / "results/02_scrna_atlas"
    results.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(source)
    adata = adata[adata.obs["pass_qc"].to_numpy()].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    profiles = patient_celltype_profiles(adata)
    markers = patient_aware_markers(profiles)
    profiles.to_csv(results / "patient_celltype_profiles.csv.gz")
    write_json(results / "patient_aware_markers.json", markers)
    sc.pp.highly_variable_genes(adata, n_top_genes=2500, flavor="cell_ranger", subset=False)
    sc.pp.pca(adata, n_comps=40, mask_var="highly_variable", random_state=17)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=17)
    sc.tl.umap(adata, random_state=17)
    adata.uns["marker_selection_unit"] = "equal-weight patient/cell-type profile"
    adata.write_h5ad(output, compression="gzip")
    sc.pl.umap(adata, color=["cell.types", "samples", "treatment.group"], show=False)
    write_json(
        results / "summary.json",
        {
            "cells": adata.n_obs,
            "genes": adata.n_vars,
            "patients": adata.obs["samples"].nunique(),
            "cell_types": adata.obs["cell.types"].value_counts().to_dict(),
            "inference_unit": "patient",
        },
    )
    plt.savefig(results / "reference_atlas.png", dpi=180, bbox_inches="tight")
    plt.close("all")


if __name__ == "__main__":
    main()
