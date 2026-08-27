#!/usr/bin/env python3
"""Download versioned public inputs declared in config/datasets.json."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "datasets.json"


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "melanoma-spatial-ici/0.1"})
    with urllib.request.urlopen(request) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle, length=1024 * 1024)
    partial.replace(destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="replace existing files")
    args = parser.parse_args()
    spec = json.loads(MANIFEST.read_text())
    for accession, dataset in spec["datasets"].items():
        for item in dataset["files"]:
            destination = ROOT / item["directory"] / item["name"]
            if destination.exists() and not args.force:
                print(f"SKIP {accession}: {destination.relative_to(ROOT)}")
                continue
            print(f"GET  {accession}: {item['url']}")
            download(item["url"], destination)
    return 0


if __name__ == "__main__":
    sys.exit(main())
