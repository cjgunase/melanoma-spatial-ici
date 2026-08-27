"""Input readers with explicit orientation and metadata checks."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def read_gse115978(counts_path: Path, annotations_path: Path) -> ad.AnnData:
    """Read the published genes-by-cells count table as cells-by-genes AnnData."""
    annotations = pd.read_csv(annotations_path).set_index("cells", drop=True)
    counts = pd.read_csv(counts_path, index_col=0)
    if not counts.columns.is_unique or not counts.index.is_unique:
        raise ValueError("GSE115978 count identifiers must be unique")
    missing_annotations = counts.columns.difference(annotations.index)
    missing_counts = annotations.index.difference(counts.columns)
    if len(missing_annotations) or len(missing_counts):
        raise ValueError(
            f"count/annotation mismatch: {len(missing_annotations)} cells lack annotations; "
            f"{len(missing_counts)} annotations lack counts"
        )
    annotations = annotations.loc[counts.columns].copy()
    matrix = sparse.csr_matrix(counts.to_numpy(dtype=np.int32, copy=False).T)
    genes = pd.DataFrame(index=counts.index.astype(str))
    genes["mt"] = genes.index.str.upper().str.startswith("MT-")
    adata = ad.AnnData(X=matrix, obs=annotations, var=genes)
    adata.obs_names.name = "cell_id"
    adata.var_names.name = "gene_symbol"
    adata.uns["source_accession"] = "GSE115978"
    adata.uns["matrix_semantics"] = "published integer count matrix"
    return adata
