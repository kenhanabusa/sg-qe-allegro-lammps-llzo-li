#!/usr/bin/env bash
set -euo pipefail
source scripts/_run_id.sh
echo "[10] RUN_ID=$RUN_ID"

# ここは最終的にQE(AIMD)に置き換える。
# まずは配管テストとして合成extxyzを生成する（後で dataset/raw/qe_llzo_li.extxyz に差し替え）
python3 dataset/scripts/make_synth_llzo_li_extxyz.py

echo "[10] OK (synthetic dataset generated)."
