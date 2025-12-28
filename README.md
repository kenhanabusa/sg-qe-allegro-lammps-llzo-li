# QE → Allegro → LAMMPS (Li/LLZO small demo)

このリポジトリは、**Quantum ESPRESSO (QE)** で得た energy/forces を元に **Allegro/NequIP** を学習し、**LAMMPS（pair_style allegro）** でモデルをロードして計算できることを確認するための「最小〜小さめデモ」です。

> ⚠️ 注意: デモ用途の小規模データで学習します。科学的に精度保証するものではありません。  
> 本格用途では、データ数・設定・検証（テストセット/外挿検知/MD安定性評価など）が必要です。

---

## 全体像（パイプライン）

1. **QE (pw.x)**: 構造に微小変位を入れたスナップショットを SCF 計算 → energy/forces を得る  
2. **extxyz**: QE 出力を `dataset/raw/*.extxyz` にまとめる  
3. **Allegro/NequIP**: `dataset/train.extxyz` を学習 → `allegro/runs/<RUN_ID>/last.ckpt`  
4. **deploy**: ckpt を TorchScript にコンパイル → `lammps/models/*.nequip.pth`  
5. **LAMMPS**: pair_style allegro でモデルロード（run 0 / short MD）

---

## RUN_ID の考え方（重要）

このリポジトリは、実行物を `RUN_ID` でまとめます（dataset/train 〜 deploy 〜 lammps が同じIDで揃うようにするため）。

- 既定では `scripts/_run_id.sh` が RUN_ID を決めます
- RUN_ID は `.run_id` に保存され、同じ run を連続実行すると同じIDを使います（再現性が上がる）

新しい run を始めたい場合：

```bash
rm -f .run_id
```

または、自分で固定したい場合：

```bash
export RUN_ID=20251228-xxxxxx
```

---

## Quickstart A（最短：QEなし / synthetic dataset）

まずは「配管（dataset → train → deploy → LAMMPS load）」が通ることを確認します。

### 1) 環境（例）

```bash
cd ~/work/sg-qe-allegro-lammps-llzo-li || exit 1
source ~/miniforge3/etc/profile.d/conda.sh
conda activate nequip2

export CUDA_VISIBLE_DEVICES=0
export OMP_NUM_THREADS=1
export LAMMPS_NP=1

# LAMMPS wrapper
export PATH="$HOME/tools/lammps/bin:$PATH"
export LMP="$HOME/tools/lammps/bin/lmp"
```

### 2) 一気通貫

```bash
rm -f .run_id
make dataset
make train
make deploy
make lammps
```

### 3) 成功確認（例）

```bash
RUN_ID="$(bash scripts/_run_id.sh)"
grep -n "NequIP/Allegro: Loading model" "lammps/runs/$RUN_ID/screen.out" | head
grep -n "ERROR:" "lammps/runs/$RUN_ID/screen.out" | head || echo "OK: no ERROR"
```

---

## Quickstart B（QEあり：LLZO 96 atoms → extxyz → train）

### 0) 前提（例：この環境で動いた設定）

```bash
# QE
export PW=/home/dl/bin/pw.x
export PSEUDO_DIR=/home/dl/qe_bench/llzo_demo_qe/pseudos

# まずは 96 atoms の LLZO（492 atoms はOOMしやすい）
export STRUCT_CIF=/home/dl/llzo_li_interface/llzo_relaxed.cif

# まずは 1MPI で確実に
export NP=1
export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=0
```

### 1) QE入力生成 → SCF → extxyz 化

```bash
cd ~/work/sg-qe-allegro-lammps-llzo-li || exit 1

# 掃除
rm -f qe/inputs/snap_*.in qe/outputs/snap_*.out
rm -f qe/snapshots/snap_*.cif
rm -rf qe/out/*

# 例：20枚（disp=0.01 Å）
python3 qe/scripts/make_qe_scf_inputs.py "$STRUCT_CIF" 20 0.01
bash qe/scripts/run_pw_scf.sh

rm -f dataset/raw/qe_llzo_li.extxyz
python3 qe/scripts/qe_out_to_extxyz.py
cp -av dataset/raw/qe_llzo_li.extxyz dataset/train.extxyz

python3 - <<'PY'
from ase.io import read
frames = read("dataset/train.extxyz", index=":")
print("frames =", len(frames), "natoms =", len(frames[0]))
print("info keys example =", list(frames[0].info.keys()))
PY
```

### 2) 学習 → deploy → LAMMPS load

```bash
# 新しいrunにしたい場合
rm -f .run_id

make train
make deploy
make lammps

RUN_ID="$(bash scripts/_run_id.sh)"
grep -n "NequIP/Allegro: Loading model" "lammps/runs/$RUN_ID/screen.out" | head
grep -n "ERROR:" "lammps/runs/$RUN_ID/screen.out" | head || echo "OK: no ERROR"
```

---

## （任意）短いMDを回す例

`make lammps` は「モデルがロードできるか」をまず確認することを優先しています。  
MDを回す場合は run dir で入力を書いて実行します。

```bash
RUN_ID="$(bash scripts/_run_id.sh)"
RUN_DIR="lammps/runs/$RUN_ID"

cat > "$RUN_DIR/in.md_short" <<'IN'
units metal
atom_style atomic
boundary p p p

read_data data.system

# type 1..4 = Li La Zr O
mass 1 6.94
mass 2 138.90547
mass 3 91.224
mass 4 15.999

pair_style allegro
pair_coeff * * model.nequip.pth Li La Zr O

neighbor 2.0 bin
neigh_modify every 1 delay 0 check yes

timestep 0.0001   # 0.1 fs（安全側）

velocity all create 10.0 12345 mom yes rot yes dist gaussian
fix int all nvt temp 10.0 10.0 0.1

thermo 10
dump d1 all custom 20 traj.lammpstrj id type x y z fx fy fz
run 200
undump d1
unfix int
write_data data.after
IN

( cd "$RUN_DIR" && "$LMP" -in in.md_short > screen.md_short.out 2>&1 )

grep -n "NequIP/Allegro: Loading model" "$RUN_DIR/screen.md_short.out" | head
grep -n "ERROR:" "$RUN_DIR/screen.md_short.out" | head || echo "OK: no ERROR"
ls -la "$RUN_DIR/traj.lammpstrj" "$RUN_DIR/data.after"
```

---

## トラブルシューティング（よくあるやつ）

### QE が GPU OOM する
- 構造が大きすぎる（492 atoms など） → **96 atoms を使う**  
- `ECUTWFC/ECUTRHO` を下げる、`NP=1` でまず通す  
- どうしてもなら CPU で回す：
  ```bash
  export CUDA_VISIBLE_DEVICES=""
  export NP=1
  ```

### mpirun が HCOLL / IB で落ちる
- `NP=1` で通す（まずこれ）
- `qe/scripts/run_pw_scf.sh` 側で prefix 検出＆回避を入れてある想定（PR で差分管理）

### extxyz 変換が「total energy が無い」などで落ちる
- QE の out に `JOB DONE` と `! total energy` があるか確認
- `qe/outputs/snap_XXX.out` が生成されているか確認

### LAMMPS が「Masses not set」や「pair_coeff が…」で落ちる
- `mass 1..4` と `pair_coeff * * model.nequip.pth Li La Zr O` を入れる
- `model.nequip.pth` が run dir から見えるか（symlink でもOK）

---

## Runbook / Notes
公開Runbook（Quickstart / 再現手順）は Notes 側にまとめる想定です。

- 例：`https://notes.server-gear.com/solutions/founding5/`（該当セクションからリンク）

