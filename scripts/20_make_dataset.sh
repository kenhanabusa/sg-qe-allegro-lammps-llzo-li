#!/usr/bin/env bash
set -euo pipefail
source scripts/_run_id.sh
echo "[20] RUN_ID=$RUN_ID"

mkdir -p dataset
# 優先順位：実データ（QE由来） > 合成データ
if ls dataset/raw/qe_*.extxyz >/dev/null 2>&1; then
  src="$(ls -1 dataset/raw/qe_*.extxyz | head -n 1)"
  cp -av "$src" dataset/train.extxyz
  echo "[20] Using QE-derived dataset: $src"
elif [ -f dataset/raw/synth_llzo_li.extxyz ]; then
  cp -av dataset/raw/synth_llzo_li.extxyz dataset/train.extxyz
  echo "[20] Using synthetic dataset (plumbing test): dataset/raw/synth_llzo_li.extxyz"
else
  echo "[20] ERROR: no dataset found under dataset/raw/"
  exit 1
fi

python3 - <<'PY'
from pathlib import Path
p=Path("dataset/train.extxyz")
assert p.exists() and p.stat().st_size>0
print("[20] dataset/train.extxyz OK")
PY
