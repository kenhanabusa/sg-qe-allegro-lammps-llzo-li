#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_ID="$(bash "$ROOT/scripts/_run_id.sh")"
RUN_DIR="$ROOT/lammps/runs/$RUN_ID"
MODEL="$ROOT/lammps/models/model.nequip.pth"

# LMP は wrapper を使うのが前提（export LMP=~/tools/lammps/bin/lmp）
LMP="${LMP:-lmp}"

echo "[50] RUN_ID=$RUN_ID"
echo "[50] RUN_DIR=$RUN_DIR"
echo "[50] MODEL=$MODEL"
echo "[50] LMP=$LMP"

mkdir -p "$RUN_DIR"

# dataset/train.extxyz が無ければエラーにする（quickstart では先に作られる想定）
test -f "$ROOT/dataset/train.extxyz" || { echo "[50] ERROR: dataset/train.extxyz missing"; exit 1; }

# ASE が使えるかチェック（extxyz -> lammps-data 変換で必要）
python3 - <<'PY' >/dev/null
import ase, ase.io
print("ASE OK", ase.__version__)
PY

# 1) extxyz (first frame) -> LAMMPS data
python3 - <<PY
from ase.io import read, write
atoms = read("$ROOT/dataset/train.extxyz", index=0)
write("$RUN_DIR/data.system", atoms, format="lammps-data", specorder=["Li","La","Zr","O"])
print("[50] Wrote data.system natoms=", len(atoms))
PY

# 2) input (run 0)
cat > "$RUN_DIR/in.allegro" <<'IN'
units metal
atom_style atomic
boundary p p p

read_data data.system

# 必須：typeごとの質量（amu）
# specorder=["Li","La","Zr","O"] → type 1..4
mass 1 6.94
mass 2 138.90547
mass 3 91.224
mass 4 15.999

pair_style allegro
pair_coeff * * model.nequip.pth Li La Zr O

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

thermo 1
thermo_style custom step pe ke etotal temp press

run 0
IN

# モデルは run dir から相対で読む（in.allegro 内は model.nequip.pth を参照）
ln -sf "$MODEL" "$RUN_DIR/model.nequip.pth"

# 実行（ログ保存）
(
  cd "$RUN_DIR"
  "$LMP" -in in.allegro -log log.lammps -echo both > screen.out 2>&1
)

# 成功チェック
grep -n "Loading model" "$RUN_DIR/screen.out" | head -n 5 || true
grep -n "ERROR:" "$RUN_DIR/screen.out" | head -n 5 || true

echo "[50] OK: $RUN_DIR/screen.out"
