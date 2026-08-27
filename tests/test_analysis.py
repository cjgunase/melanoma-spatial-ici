import numpy as np
import pandas as pd

from melanoma_spatial_ici.analysis import gene_set_score, patient_aware_markers


def test_gene_set_score_ignores_unavailable_genes():
    frame = pd.DataFrame({"A": [1.0, 2.0, 3.0], "B": [2.0, 2.0, 2.0]})
    score = gene_set_score(frame, ["A", "missing"])
    assert score.notna().all()
    assert np.isclose(score.mean(), 0)


def test_markers_use_patient_profiles():
    index = pd.MultiIndex.from_tuples(
        [("p1", "T"), ("p2", "T"), ("p1", "M"), ("p2", "M")],
        names=["patient", "cell_type"],
    )
    profiles = pd.DataFrame({"TGENE": [4, 5, 0, 0], "MGENE": [0, 0, 3, 4]}, index=index)
    markers = patient_aware_markers(profiles, n=1)
    assert markers == {"M": ["MGENE"], "T": ["TGENE"]}
