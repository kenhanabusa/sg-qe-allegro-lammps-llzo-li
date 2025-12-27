#!/usr/bin/env bash
set -euo pipefail
source scripts/_run_id.sh
echo "[50] RUN_ID=$RUN_ID"

# lmpバイナリ（後でビルドする）
LMP="${LMP:-$HOME/tools/lammps/bin/lmp}"
if [ ! -x "$LMP" ]; then
  echo "[50] ERROR: LAMMPS not built yet. Build it with:"
  echo "  bash tools/build_lammps_pair_allegro.sh"
  exit 1
fi

# compiled model を探す
model="$(ls -1 lammps/model/llzo_li_${RUN_ID}.nequip.pth 2>/dev/null || true)"
[ -n "$model" ] || { echo "[50] ERROR: compiled model not found. Run: make deploy"; exit 1; }

# dataファイル生成（dataset/train.extxyzの先頭フレームから）
python3 lammps/scripts/extxyz_to_data.py dataset/train.extxyz lammps/in/llzo_li.data

# input中のMODEL差し替え（run dirへコピーしてから実行）
RUN_DIR="lammps/runs/${RUN_ID}"
mkdir -p "$RUN_DIR"/{in,out,model}
cp -av "$model" "$RUN_DIR/model/MODEL.nequip.pth"
cp -av lammps/in/in.allegro "$RUN_DIR/in/in.allegro"
cp -av lammps/in/llzo_li.data "$RUN_DIR/in/llzo_li.data"

# MODELパス差し替え
sed -i 's|lammps/model/MODEL.nequip.pth|model/MODEL.nequip.pth|g' "$RUN_DIR/in/in.allegro"
# dataパス差し替え
sed -i 's|lammps/in/llzo_li.data|in/llzo_li.data|g' "$RUN_DIR/in/in.allegro"
# dump先
mkdir -p "$RUN_DIR/out"
sed -i 's|lammps/out/|out/|g' "$RUN_DIR/in/in.allegro"

echo "[50] running LAMMPS..."
cd "$RUN_DIR"

# pair_nequip_allegro: AllegroはKokkos+MPIで動かせる（例あり） [oai_citation:4‡GitHub](https://github.com/mir-group/pair_nequip_allegro)
"$LMP" -sf kk -k on g 1 -pk kokkos newton on neigh half -in in/in.allegro | tee out/log.lammps

echo "[50] OK: $RUN_DIR"
