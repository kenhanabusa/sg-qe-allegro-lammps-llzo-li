#!/usr/bin/env bash
set -euo pipefail
source scripts/_run_id.sh
echo "[30] RUN_ID=$RUN_ID"

# 必須コマンド確認
command -v nequip-train >/dev/null 2>&1 || { echo "[30] ERROR: nequip-train not found (install nequip-allegro/nequip in your env)"; exit 1; }
python3 -c "import yaml" >/dev/null 2>&1 || { echo "[30] ERROR: PyYAML missing"; exit 1; }

# upstreamのtutorial.yaml取得
if [ ! -f tools/allegro_tutorial.yaml ]; then
  rm -rf tools/allegro_upstream
  git clone --depth=1 https://github.com/mir-group/allegro.git tools/allegro_upstream
  cp -av tools/allegro_upstream/configs/tutorial.yaml tools/allegro_tutorial.yaml
fi

# config生成（データパス・type_namesを上書き）
python3 tools/prepare_allegro_config.py tools/allegro_tutorial.yaml allegro/configs/train.yaml dataset/train.extxyz "Li,La,Zr,O"

# 実行
mkdir -p allegro/runs
RUN_DIR="allegro/runs/${RUN_ID}"
mkdir -p "$RUN_DIR"
echo "[30] train dir: $RUN_DIR"

# まずは最小実行（詳細はtrain.yamlに依存）
nequip-train allegro/configs/train.yaml --run-dir "$RUN_DIR"

echo "[30] OK"
