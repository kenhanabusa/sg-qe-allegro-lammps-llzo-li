#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

RUN_ID="$(bash scripts/_run_id.sh)"
echo "[40] RUN_ID=$RUN_ID"

TRAIN_DIR="allegro/runs/${RUN_ID}"

# RUN_ID のディレクトリが無い/空のときは最新を拾う（ズレ耐性）
if [ ! -d "$TRAIN_DIR" ]; then
  if ls -1d allegro/runs/* >/dev/null 2>&1; then
    RUN_ID="$(ls -1dt allegro/runs/* | head -n1 | xargs -n1 basename)"
    TRAIN_DIR="allegro/runs/${RUN_ID}"
    echo "[40] WARN: $ROOT/allegro/runs/${RUN_ID} not found. Use latest RUN_ID=$RUN_ID"
  fi
fi

# checkpoint を選ぶ（last.ckpt 優先）
CKPT=""
if [ -f "$TRAIN_DIR/last.ckpt" ]; then
  CKPT="$TRAIN_DIR/last.ckpt"
else
  CKPT="$(ls -1t "$TRAIN_DIR"/*.ckpt 2>/dev/null | head -n1 || true)"
fi

if [ -z "${CKPT:-}" ]; then
  echo "[40] ERROR: no .ckpt found under: $TRAIN_DIR" >&2
  ls -la "$TRAIN_DIR" || true
  exit 1
fi
echo "[40] CKPT=$CKPT"

command -v nequip-compile >/dev/null 2>&1 || {
  echo "[40] ERROR: nequip-compile not found. Activate env (conda activate nequip2) first." >&2
  exit 1
}

mkdir -p lammps/models

# LAMMPS(pair_style allegro) は .nequip.pth (TorchScript) を受け付ける  [oai_citation:1‡NequIP](https://nequip.readthedocs.io/projects/allegro/en/latest/guide/lammps.html)
OUT="lammps/models/${RUN_ID}.nequip.pth"
echo "[40] compile -> $OUT"
nequip-compile "$CKPT" "$OUT" --device cuda --mode torchscript

# 実行側は固定ファイル名で参照できるようにしておく
ln -sf "$(basename "$OUT")" lammps/models/model.nequip.pth

echo "[40] OK: built"
ls -la lammps/models | sed -n '1,20p'
