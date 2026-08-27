#!/usr/bin/env python3
"""Identify recurrent tumor-immune neighborhoods and sensitivity to k."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans

from melanoma_spatial_ici.analysis import write_json

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    result_dir = ROOT / "results/04_spatial_neighborhoods"
    result_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(ROOT / "data/processed/GSE300445_mapped.h5ad")
    abundance = [c for c in adata.obs if c.startswith("map_lsq_") and not c.endswith("_?")]
    features = adata.obs[abundance].fillna(0)
    model = KMeans(n_clusters=5, random_state=17, n_init=20).fit(features)
    adata.obs["neighborhood"] = pd.Categorical(model.labels_.astype(str))
    centroids = pd.DataFrame(model.cluster_centers_, columns=abundance)
    centroids.index.name = "neighborhood"
    centroids.to_csv(result_dir / "neighborhood_centroids.csv")
    sensitivity = []
    for k in (4, 6, 10):
        for section in adata.obs["section"].unique():
            mask = adata.obs["section"].eq(section).to_numpy()
            coords = adata.obsm["spatial"][mask]
            _, indices = cKDTree(coords).query(coords, k=min(k + 1, len(coords)))
            local = features.loc[mask].to_numpy()[indices[:, 1:]].mean(axis=1)
            frame = pd.DataFrame(local, columns=abundance)
            frame["resistance_score"] = adata.obs.loc[mask, "resistance_score"].to_numpy()
            sensitivity.append(
                {
                    "section": section,
                    "k": k,
                    "resistance_macrophage_r": frame["resistance_score"].corr(
                        frame.get("map_lsq_Macrophage")
                    ),
                    "resistance_tcell_r": frame["resistance_score"].corr(
                        frame.filter(regex=r"map_lsq_T\.").sum(axis=1)
                    ),
                }
            )
    sensitivity_frame = pd.DataFrame(sensitivity)
    sensitivity_frame.to_csv(result_dir / "sensitivity.csv", index=False)
    composition = pd.crosstab(adata.obs["section"], adata.obs["neighborhood"], normalize="index")
    composition.to_csv(result_dir / "section_neighborhood_composition.csv")
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    sns.heatmap(centroids, cmap="mako", ax=axes[0])
    sns.heatmap(composition, cmap="viridis", annot=True, fmt=".2f", ax=axes[1])
    axes[0].set_title("Reference-mapped neighborhood centroids")
    axes[1].set_title("Neighborhood recurrence by section")
    fig.tight_layout()
    fig.savefig(result_dir / "neighborhood_summary.png", dpi=180)
    write_json(
        result_dir / "summary.json",
        {
            "sections": int(adata.obs["section"].nunique()),
            "spots": adata.n_obs,
            "neighborhoods": 5,
            "sensitivity_k": [4, 6, 10],
            "claim_scope": "exploratory spatial recurrence; no response labels available",
        },
    )
    adata.write_h5ad(ROOT / "data/processed/GSE300445_neighborhoods.h5ad", compression="gzip")


if __name__ == "__main__":
    main()
