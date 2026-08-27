#!/usr/bin/env python3
"""Test the predeclared resistance score in patient-collapsed GeoMx cohorts."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_auc_score

from melanoma_spatial_ici.analysis import (
    RESISTANCE_GENES,
    gene_set_score,
    read_geomx_workbook,
    write_json,
)

ROOT = Path(__file__).resolve().parents[1]


def boolean_outcome(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip().str.lower()
    mapped = text.map(
        {
            "1": True,
            "1.0": True,
            "yes": True,
            "true": True,
            "0": False,
            "0.0": False,
            "no": False,
            "false": False,
        }
    )
    return mapped


def main() -> None:
    raw = ROOT / "data/raw/GSE233305_validation"
    result_dir = ROOT / "results/05_geomx_validation"
    result_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for path in sorted(raw.glob("*.xlsx")):
        expression, metadata = read_geomx_workbook(path)
        cohort = "ITX2" if "ITX2" in path.name else "ITX1"
        compartment = next(x for x in ("CD45", "CD68", "S100B") if x in path.name)
        scores = gene_set_score(expression, RESISTANCE_GENES)
        patient_column = next(
            (c for c in ("CPID", "SPID", "SURGICAL_NUMBER") if c in metadata), None
        )
        if patient_column is None:
            raise ValueError(f"patient identifier unavailable in {path.name}")
        response_column = next(
            (c for c in ("DCB6", "disease_control_", "response") if c in metadata), None
        )
        for segment, score in scores.items():
            records.append(
                {
                    "cohort": cohort,
                    "compartment": compartment,
                    "segment": segment,
                    "patient": metadata.at[segment, patient_column],
                    "resistance_score": score,
                    "response": metadata.at[segment, response_column]
                    if response_column
                    else np.nan,
                }
            )
    segments = pd.DataFrame(records).dropna(subset=["patient", "resistance_score"])
    patients = segments.groupby(["cohort", "compartment", "patient"], as_index=False).agg(
        resistance_score=("resistance_score", "mean"),
        response=("response", "first"),
        n_segments=("segment", "nunique"),
    )
    patients["response_binary"] = boolean_outcome(patients["response"])
    patients.to_csv(result_dir / "patient_scores.csv", index=False)
    tests = []
    for (cohort, compartment), frame in patients.groupby(["cohort", "compartment"]):
        valid = frame.dropna(subset=["response_binary"])
        response_mask = valid["response_binary"].astype(bool)
        responders = valid.loc[response_mask, "resistance_score"]
        nonresponders = valid.loc[~response_mask, "resistance_score"]
        if len(responders) >= 2 and len(nonresponders) >= 2:
            statistic, pvalue = mannwhitneyu(nonresponders, responders, alternative="two-sided")
            auc = roc_auc_score(~response_mask, valid["resistance_score"])
        else:
            statistic = pvalue = auc = np.nan
        tests.append(
            {
                "cohort": cohort,
                "compartment": compartment,
                "patients": len(valid),
                "responders": len(responders),
                "nonresponders": len(nonresponders),
                "mannwhitney_u": statistic,
                "pvalue": pvalue,
                "auc_nonresponse": auc,
                "confirmatory": cohort == "ITX2" and compartment == "S100B",
            }
        )
    tests_frame = pd.DataFrame(tests)
    tests_frame.to_csv(result_dir / "tests.csv", index=False)
    plot = patients[patients["cohort"].eq("ITX2")].dropna(subset=["response_binary"])
    fig, ax = plt.subplots(figsize=(7, 4))
    sns.boxplot(plot, x="compartment", y="resistance_score", hue="response_binary", ax=ax)
    sns.stripplot(
        plot,
        x="compartment",
        y="resistance_score",
        hue="response_binary",
        dodge=True,
        color="black",
        alpha=0.55,
        legend=False,
        ax=ax,
    )
    ax.set_title("Independent ITX2 patient-level resistance score")
    fig.tight_layout()
    fig.savefig(result_dir / "validation_summary.png", dpi=180)
    write_json(
        result_dir / "summary.json",
        {
            "hypothesis": "Higher melanoma dedifferentiation score associates with nonresponse in S100B segments",
            "genes": RESISTANCE_GENES,
            "experimental_unit": "patient",
            "discovery_cohort": "ITX1",
            "independent_validation_cohort": "ITX2",
            "multiplicity": "one predeclared primary ITX2/S100B test; other rows are descriptive",
        },
    )


if __name__ == "__main__":
    main()
