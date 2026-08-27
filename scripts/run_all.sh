#!/usr/bin/env bash
set -euo pipefail
uv run python scripts/download_data.py
uv run python scripts/inventory_data.py --strict
uv run python scripts/01_ingest_scrna.py
uv run python scripts/02_build_scrna_atlas.py
uv run python scripts/03_map_visium.py
uv run python scripts/04_spatial_neighborhoods.py
uv run python scripts/05_validate_geomx.py
uv run python scripts/06_build_report.py
uv run pytest
uv run ruff check .
