from pathlib import Path

import pandas as pd
import pytest

from melanoma_spatial_ici.io import read_gse115978


def test_reader_preserves_cell_and_gene_orientation(tmp_path: Path):
    counts = pd.DataFrame({"cell_a": [1, 0], "cell_b": [2, 3]}, index=["GENE1", "MT-X"])
    annotations = pd.DataFrame(
        {
            "cells": ["cell_b", "cell_a"],
            "samples": ["p2", "p1"],
            "cell.types": ["T", "Mal"],
        }
    )
    counts_path = tmp_path / "counts.csv.gz"
    annotations_path = tmp_path / "annotations.csv.gz"
    counts.to_csv(counts_path)
    annotations.to_csv(annotations_path, index=False)
    adata = read_gse115978(counts_path, annotations_path)
    assert adata.shape == (2, 2)
    assert adata.obs_names.tolist() == ["cell_a", "cell_b"]
    assert adata.obs["samples"].tolist() == ["p1", "p2"]
    assert adata.var["mt"].tolist() == [False, True]


def test_reader_rejects_metadata_mismatch(tmp_path: Path):
    counts = pd.DataFrame({"cell_a": [1]}, index=["GENE1"])
    annotations = pd.DataFrame({"cells": ["different_cell"]})
    counts_path = tmp_path / "counts.csv.gz"
    annotations_path = tmp_path / "annotations.csv.gz"
    counts.to_csv(counts_path)
    annotations.to_csv(annotations_path, index=False)
    with pytest.raises(ValueError, match="mismatch"):
        read_gse115978(counts_path, annotations_path)
