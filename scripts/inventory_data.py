#!/usr/bin/env python3
"""Validate source inputs and write a deterministic checksum inventory."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "config" / "datasets.json"
OUTPUT = ROOT / "data" / "source_inventory.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    spec = json.loads(MANIFEST.read_text())
    records, errors = [], []
    for accession, dataset in spec["datasets"].items():
        for item in dataset["files"]:
            path = ROOT / item["directory"] / item["name"]
            record = {
                "accession": accession,
                "modality": dataset["modality"],
                "role": dataset["role"],
                "path": str(path.relative_to(ROOT)),
                "url": item["url"],
                "exists": path.exists(),
            }
            if path.exists():
                record["bytes"] = path.stat().st_size
                record["sha256"] = sha256(path)
                expected = item.get("expected_bytes")
                record["size_matches"] = expected is None or expected == record["bytes"]
                if not record["size_matches"]:
                    errors.append(f"size mismatch: {record['path']}")
            else:
                errors.append(f"missing: {record['path']}")
            records.append(record)
    OUTPUT.write_text(json.dumps({"files": records}, indent=2) + "\n")
    for record in records:
        status = "OK" if record.get("exists") and record.get("size_matches", True) else "FAIL"
        print(f"{status:4} {record['path']}")
    if errors:
        print("\n" + "\n".join(errors), file=sys.stderr)
    return 1 if args.strict and errors else 0


if __name__ == "__main__":
    sys.exit(main())
