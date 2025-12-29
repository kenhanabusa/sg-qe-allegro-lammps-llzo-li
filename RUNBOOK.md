RUNBOOK: QE → Allegro/NequIP → LAMMPS (Li/LLZO small demo)

このRunbookは「初見の人が、できるだけコピペで再現できる」ことを目的にしています。
デモのゴールは「QEで得た energy/forces から学習したモデルを LAMMPS がロードし、run0〜短いMDが動く」ことです。

0) 前提（最初に確認）

0-1) conda 環境（例：nequip2）

conda activate nequip2
command -v nequip-train
command -v nequip-compile
python3 -c "import ase, yaml"

0-2) LAMMPS（pair_style allegro）

export PATH="$HOME/tools/lammps/bin:$PATH"
export LMP="$HOME/tools/lammps/bin/lmp"
lmp -h | head -n 5

0-3) QE（pw.x）と擬ポテンシャル（UPF）

（例：あなたの環境で動いたパス）
PW=/home/dl/bin/pw.x
PSEUDO_DIR=/home/dl/qe_bench/llzo_demo_qe/pseudos

0-4) まず小さい構造を使う（重要）

492 atoms など大きい構造は QE GPU OOM しやすいので、まず 96 atoms の LLZO から始めます。

STRUCT_CIF=/home/dl/llzo_li_interface/llzo_relaxed.cif
（確認）
python3 - <<'PY'
from ase.io import read
a = read("/home/dl/llzo_li_interface/llzo_relaxed.cif")
print("natoms=", len(a))
print("elements=", sorted(set(a.get_chemical_symbols())))
PY

---

1) Quickstart A（QEなし：まず配線チェック）

目的：dataset→train→deploy→lammps が "通る" ことだけ確認

conda activate nequip2
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export LAMMPS_NP=1
export PATH="$HOME/tools/lammps/bin:$PATH"
export LMP="$HOME/tools/lammps/bin/lmp"

rm -f .run_id
make dataset
make train
make deploy
make lammps

RUN_ID="$(bash scripts/_run_id.sh)"
grep -n "NequIP/Allegro: Loading model" "lammps/runs/$RUN_ID/screen.out" | head
grep -n "ERROR:" "lammps/runs/$RUN_ID/screen.out" | head || echo "OK: no ERROR"

---

2) Quickstart B（QEあり：LLZO small → extxyz → train → deploy → lammps）

2-1) 環境変数

conda activate nequip2
export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export NP=1

export PW=/home/dl/bin/pw.x
export PSEUDO_DIR=/home/dl/qe_bench/llzo_demo_qe/pseudos
export STRUCT_CIF=/home/dl/llzo_li_interface/llzo_relaxed.cif

（まず軽め。必要なら上げる）
export ECUTWFC=40
export ECUTRHO=320
export DEGAUSS=0.02
export CONV_THR=1e-6
export MIXING_BETA=0.3

2-2) QE 入力生成（20枚、微小変位）

rm -f qe/inputs/snap_*.in qe/outputs/snap_*.out
rm -f qe/snapshots/snap_*.cif
rm -rf qe/out/*

python3 qe/scripts/make_qe_scf_inputs.py "$STRUCT_CIF" 20 0.01

2-3) QE 実行（SCF）

bash qe/scripts/run_pw_scf.sh

成功チェック（最低これだけ）
grep -n "JOB DONE" qe/outputs/snap_000.out | tail -n 3
grep -n "! *total energy" qe/outputs/snap_000.out | tail -n 3
grep -n "Forces acting on atoms" qe/outputs/snap_000.out | head -n 2

2-4) QE出力 → extxyz 変換

rm -f dataset/raw/qe_llzo_li.extxyz
python3 qe/scripts/qe_out_to_extxyz.py

（frames確認）
python3 - <<'PY'
from ase.io import read
frames = read("dataset/raw/qe_llzo_li.extxyz", index=":")
print("frames =", len(frames), "natoms =", len(frames[0]))
PY

2-5) 学習用データとして採用

cp -av dataset/raw/qe_llzo_li.extxyz dataset/train.extxyz

2-6) train → deploy → lammps

rm -f .run_id
export RUN_ID="$(date +%Y%m%d-%H%M%S)"
echo "RUN_ID=$RUN_ID"

bash scripts/30_train_allegro.sh
make deploy

export PATH="$HOME/tools/lammps/bin:$PATH"
export LMP="$HOME/tools/lammps/bin/lmp"
export LAMMPS_NP=1
make lammps

（成功判定）
grep -n "NequIP/Allegro: Loading model" "lammps/runs/$RUN_ID/screen.out" | head
grep -n "ERROR:" "lammps/runs/$RUN_ID/screen.out" | head || echo "OK: no ERROR"

---

3) 追加：短いMD（200 step）を回す（任意）

RUN_DIR="lammps/runs/$RUN_ID"

cat > "$RUN_DIR/in.md_short" <<'IN'
units metal
atom_style atomic
boundary p p p

read_data data.system

mass 1 6.94
mass 2 138.90547
mass 3 91.224
mass 4 15.999

pair_style allegro
pair_coeff * * model.nequip.pth Li La Zr O

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

timestep 0.0001
velocity all create 10.0 12345 mom yes rot yes dist gaussian
fix int all nvt temp 10.0 10.0 0.1

thermo 10
dump d all custom 10 traj.lammpstrj id type x y z fx fy fz
run 200

write_data data.after
IN

( cd "$RUN_DIR" && "$LMP" -in in.md_short > screen.md_short.out 2>&1 )
tail -n 40 "$RUN_DIR/screen.md_short.out"
ls -la "$RUN_DIR/traj.lammpstrj" "$RUN_DIR/data.after"

---

4) Troubleshooting（よくある）

QE が GPU OOM する
	•	構造が大きい（492 atomsなど）→ まず 96 atoms を使う
	•	ECUTWFC/ECUTRHO を下げる
	•	NP=1 でまず通す
	•	最終手段：CPUで回す（遅いが通る）
export CUDA_VISIBLE_DEVICES=""
export NP=1

mpirun が HCOLL / IB で落ちる
	•	NP=1 でまず通す
	•	qe/scripts/run_pw_scf.sh は pw.x の libmpi.so から prefix を推定して mpirun に反映する想定（環境差を吸収）

extxyz 変換で total energy が見つからない
	•	qe/outputs/snap_000.out に JOB DONE と ! total energy があるか確認
	•	out が途中で落ちてないか確認（OOMやMPIエラー）

trainで RUN_ID がずれる
	•	rm -f .run_id で新規 run にする
	•	export RUN_ID=… で固定してから train→deploy→lammps を通す

