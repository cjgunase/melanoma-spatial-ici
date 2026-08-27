"""Patient-aware reference, spatial mapping, and validation utilities."""

from __future__ import annotations

import json
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmread

RESISTANCE_GENES = ["AXL", "NGFR", "WNT5A", "EGFR", "FN1", "VIM"]
T_CELL_GENES = ["CD3D", "CD3E", "CD8A", "TRAC", "CCL5"]
MYELOID_GENES = ["CD68", "CD163", "C1QA", "C1QB", "LST1"]


def patient_celltype_profiles(adata: ad.AnnData) -> pd.DataFrame:
    """Return equal-weight patient/cell-type log-normalized mean profiles."""
    groups = adata.obs[["samples", "cell.types"]].astype(str).agg("|".join, axis=1)
    rows = []
    for group in sorted(groups.unique()):
        mask = groups.to_numpy() == group
        mean = np.asarray(adata.X[mask].mean(axis=0)).ravel()
        patient, cell_type = group.split("|", 1)
        rows.append(pd.Series(mean, index=adata.var_names, name=(patient, cell_type)))
    result = pd.DataFrame(rows)
    result.index = pd.MultiIndex.from_tuples(result.index, names=["patient", "cell_type"])
    return result


def patient_aware_markers(profiles: pd.DataFrame, n: int = 30) -> dict[str, list[str]]:
    """Select markers from patient-level profiles, avoiding cell-count weighting."""
    markers: dict[str, list[str]] = {}
    cell_types = profiles.index.get_level_values("cell_type")
    for cell_type in sorted(cell_types.unique()):
        target = profiles.loc[cell_types == cell_type].mean(axis=0)
        other = profiles.loc[cell_types != cell_type].mean(axis=0)
        expressed = target[target > 0.05].index
        markers[cell_type] = (
            (target.loc[expressed] - other.loc[expressed]).nlargest(n).index.tolist()
        )
    return markers


def extract_visium_archive(archive: Path, destination: Path) -> None:
    """Safely extract the declared GEO archive without path traversal."""
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as handle:
        root = destination.resolve()
        for member in handle.getmembers():
            if not (root / member.name).resolve().is_relative_to(root):
                raise ValueError(f"unsafe archive member: {member.name}")
        handle.extractall(destination, filter="data")


def read_visium_sample(directory: Path, prefix: str) -> ad.AnnData:
    """Read one GEO-exported Visium sample and retain in-tissue spots."""
    matrix = mmread(directory / f"{prefix}_processed_matrix.mtx.gz").tocsr().T
    features = pd.read_csv(directory / f"{prefix}_processed_features.tsv.gz", sep="\t", header=None)
    barcodes = pd.read_csv(
        directory / f"{prefix}_processed_barcodes.tsv.gz", sep="\t", header=None
    )[0].astype(str)
    positions = pd.read_csv(directory / f"{prefix}_processed_tissue_positions.csv.gz")
    if "barcode" not in positions.columns:
        positions.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]
    positions = positions.set_index("barcode").reindex(barcodes)
    genes = features.iloc[:, 1].astype(str)
    var = pd.DataFrame(index=pd.Index(genes, name="gene_symbol"))
    var["gene_id"] = features.iloc[:, 0].astype(str).to_numpy()
    var.index = pd.Index(ad.utils.make_index_unique(var.index.astype(str)))
    obs = positions.copy()
    obs.index = pd.Index(barcodes, name="barcode")
    sample = prefix.split("_", 1)[1]
    obs["sample"] = sample
    result = ad.AnnData(X=matrix, obs=obs, var=var)
    return result[result.obs["in_tissue"].fillna(0).astype(int).eq(1)].copy()


def gene_set_score(frame: pd.DataFrame, genes: list[str]) -> pd.Series:
    """Mean within-sample standardized expression for available genes."""
    available = frame.columns.intersection(genes)
    if available.empty:
        return pd.Series(np.nan, index=frame.index)
    values = frame.loc[:, available]
    std = values.std(axis=0).replace(0, np.nan)
    return ((values - values.mean(axis=0)) / std).mean(axis=1)


def read_geomx_workbook(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read normalized expression (segments x genes) and aligned clinical data."""
    book = pd.ExcelFile(path)
    metadata_sheet = "ClinInfo" if "ClinInfo" in book.sheet_names else "SegmentProperties"
    metadata = pd.read_excel(path, sheet_name=metadata_sheet)
    expression = pd.read_excel(path, sheet_name="TargetCountMatrix").set_index("TargetName").T
    expression.index = expression.index.astype(str)
    key = "SegmentDisplayName"
    metadata[key] = metadata[key].astype(str)
    metadata = metadata.drop_duplicates(key).set_index(key).reindex(expression.index)
    if metadata.index.isna().any():
        raise ValueError(f"metadata alignment failed for {path.name}")
    return expression.apply(pd.to_numeric), metadata


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, default=str) + "\n")
