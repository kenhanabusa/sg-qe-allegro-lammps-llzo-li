#!/usr/bin/env bash
set -euo pipefail
source scripts/_run_id.sh
echo "[40] RUN_ID=$RUN_ID"

command -v nequip-compile >/dev/null 2>&1 || { echo "[40] ERROR: nequip-compile not found"; exit 1; }

RUN_DIR="allegro/runs/${RUN_ID}"

# ckpt探し（nequipの出力は環境で違うので最大公約数）
ckpt="$(find "$RUN_DIR" -maxdepth 3 -type f \( -name "*.pth" -o -name "*.pt" \) | head -n 1 || true)"
[ -n "$ckpt" ] || { echo "[40] ERROR: no checkpoint found under $RUN_DIR"; exit 1; }

mkdir -p lammps/model
OUT="lammps/model/llzo_li_${RUN_ID}.nequip.pth"

# Allegro docs: TorchScript compile for LAMMPS is nequip-compile --mode torchscript  [oai_citation:1‡NequIP](https://nequip.readthedocs.io/projects/allegro/en/latest/guide/lammps.html)
nequip-compile "$ckpt" "$OUT" --device cuda --mode torchscript

echo "[40] compiled model: $OUT"
