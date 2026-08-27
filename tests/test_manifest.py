import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_dataset_manifest_has_unique_destinations():
    manifest = json.loads((ROOT / "config" / "datasets.json").read_text())
    destinations = []
    for dataset in manifest["datasets"].values():
        for item in dataset["files"]:
            destinations.append((item["directory"], item["name"]))
            assert item["url"].startswith("https://ftp.ncbi.nlm.nih.gov/geo/")
    assert len(destinations) == len(set(destinations))
