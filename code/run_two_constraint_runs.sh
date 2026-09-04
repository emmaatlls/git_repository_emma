#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="/home/emma/Dokumente/thesis/Model_generation_curation"

echo "==> Run 1: exclude_IPPTonly.csv -> Draft_models/excludeIPPTonly"
python "$ROOT_DIR/code/01_carve_draft_models.py" \
  --output-dir "$ROOT_DIR/Draft_models/excludeIPPTonly" \
  --soft-exclude "$ROOT_DIR/exclude_IPPTonly.csv"

python "$ROOT_DIR/code/02_create_MEMOTE_reports.py" \
  --model-dir "$ROOT_DIR/Draft_models/excludeIPPTonly"

echo "==> Run 2: exclude_IPPTandmodel.csv -> Draft_models/excludeModelandIPPT"
python "$ROOT_DIR/code/01_carve_draft_models.py" \
  --output-dir "$ROOT_DIR/Draft_models/excludeModelandIPPT" \
  --soft-exclude "$ROOT_DIR/exclude_IPPTandmodel.csv"

python "$ROOT_DIR/code/02_create_MEMOTE_reports.py" \
  --model-dir "$ROOT_DIR/Draft_models/excludeModelandIPPT"
