#!/usr/bin/env python3
"""Ingest Visium sections and map reference cell states with two methods."""

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from melanoma_spatial_ici.analysis import (
    RESISTANCE_GENES,
    extract_visium_archive,
    gene_set_score,
    read_visium_sample,
    write_json,
)

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    raw = ROOT / "data/raw/GSE300445_Visium"
    extracted = raw / "extracted"
    extract_visium_archive(raw / "GSE300445_RAW.tar", extracted)
    markers = json.loads((ROOT / "results/02_scrna_atlas/patient_aware_markers.json").read_text())
    profiles = (
        pd.read_csv(
            ROOT / "results/02_scrna_atlas/patient_celltype_profiles.csv.gz", index_col=[0, 1]
        )
        .groupby(level=1)
        .mean()
    )
    prefixes = sorted(
        {p.name.split("_processed_", 1)[0] for p in extracted.glob("*_matrix.mtx.gz")}
    )
    outputs, summary = [], []
    for prefix in prefixes:
        sample = read_visium_sample(extracted, prefix)
        sample.obs_names = pd.Index([f"{prefix}:{name}" for name in sample.obs_names])
        sc.pp.calculate_qc_metrics(sample, percent_top=None, log1p=False, inplace=True)
        sample = sample[sample.obs["n_genes_by_counts"] >= 100].copy()
        sample.layers["counts"] = sample.X.copy()
        sc.pp.normalize_total(sample, target_sum=1e4)
        sc.pp.log1p(sample)
        dense = pd.DataFrame(sample.X.toarray(), index=sample.obs_names, columns=sample.var_names)
        common = profiles.columns.intersection(sample.var_names)
        variance = profiles.loc[:, common].var(axis=0)
        selected = variance.nlargest(min(500, len(variance))).index
        reference = profiles.loc[:, selected].T.to_numpy()
        coefficients = np.linalg.lstsq(reference, dense.loc[:, selected].T.to_numpy(), rcond=None)[
            0
        ].T
        coefficients = np.clip(coefficients, 0, None)
        coefficients /= np.maximum(coefficients.sum(axis=1, keepdims=True), 1e-12)
        for i, cell_type in enumerate(profiles.index):
            sample.obs[f"map_lsq_{cell_type}"] = coefficients[:, i]
            genes = pd.Index(markers[cell_type]).intersection(sample.var_names)
            sample.obs[f"map_score_{cell_type}"] = gene_set_score(dense, list(genes)).to_numpy()
        sample.obs["resistance_score"] = gene_set_score(dense, RESISTANCE_GENES).to_numpy()
        sample.obsm["spatial"] = sample.obs[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
        outputs.append(sample)
        summary.append(
            {
                "sample": prefix,
                "spots": sample.n_obs,
                "genes": sample.n_vars,
                "median_counts": float(sample.obs["total_counts"].median()),
                "mapping_methods": ["projected_reference_least_squares", "patient_marker_score"],
            }
        )
        sample.write_h5ad(ROOT / f"data/processed/{prefix}_mapped.h5ad", compression="gzip")
    ad.concat(outputs, join="outer", label="section", keys=prefixes, index_unique=None).write_h5ad(
        ROOT / "data/processed/GSE300445_mapped.h5ad", compression="gzip"
    )
    write_json(ROOT / "results/03_visium_mapping/summary.json", summary)


if __name__ == "__main__":
    main()
