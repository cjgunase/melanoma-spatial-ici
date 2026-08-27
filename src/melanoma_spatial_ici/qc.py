"""Quality-control calculations for the single-cell reference."""

from __future__ import annotations

import anndata as ad
import scanpy as sc


def calculate_scrna_qc(
    adata: ad.AnnData,
    min_genes: int = 500,
    max_genes: int = 9000,
    min_counts: int = 1000,
    max_pct_mt: float = 20.0,
) -> ad.AnnData:
    """Calculate metrics and flag cells; do not silently discard observations."""
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    base_pass = (
        adata.obs["n_genes_by_counts"].between(min_genes, max_genes)
        & (adata.obs["total_counts"] >= min_counts)
    )
    mitochondrial_genes_present = int(adata.var["mt"].sum())
    adata.obs["pass_qc"] = (
        base_pass & (adata.obs["pct_counts_mt"] <= max_pct_mt)
        if mitochondrial_genes_present
        else base_pass
    )
    adata.uns["qc_thresholds"] = {
        "min_genes": min_genes,
        "max_genes": max_genes,
        "min_counts": min_counts,
        "max_pct_mt": max_pct_mt,
        "mitochondrial_genes_present": mitochondrial_genes_present,
        "mitochondrial_filter_applied": bool(mitochondrial_genes_present),
    }
    return adata
