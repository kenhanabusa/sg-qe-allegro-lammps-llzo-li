import numpy as np
from pathlib import Path

OUT = Path("dataset/raw/synth_llzo_li.extxyz")
OUT.parent.mkdir(parents=True, exist_ok=True)

# Li/La/Zr/O の「見た目だけ」小さめ界面っぽいボックス（合成データ＝配管テスト用）
type_names = ["Li", "La", "Zr", "O"]
counts = {"Li": 40, "La": 8, "Zr": 8, "O": 80}  # 合計136 atoms（小さめ）
N = sum(counts.values())
cell = np.diag([16.0, 16.0, 18.0])  # Å

# 多体っぽく見えるように、ゆるい擬似ポテンシャル（単なるデモ用）
rng = np.random.default_rng(123)
pos0 = rng.random((N, 3)) @ cell

# 簡単な「斥力＋弱い引力」風の擬似エネルギー/力（物理意味は無し：配管確認用）
def energy_forces(pos):
    # speciesごとの係数（適当）
    eps = np.array([0.02, 0.05, 0.05, 0.03])  # Li,La,Zr,O
    sig = np.array([2.0, 2.8, 2.6, 2.2])
    # type index
    types = []
    for t in type_names:
        types += [t]*counts[t]
    ti = np.array([type_names.index(t) for t in types], dtype=int)

    E = 0.0
    F = np.zeros_like(pos)
    # O(N^2)でも小さいのでOK
    for i in range(N):
        ri = pos[i]
        for j in range(i+1, N):
            rj = pos[j]
            dr = ri - rj
            r = np.linalg.norm(dr) + 1e-9
            # mixing
            eij = np.sqrt(eps[ti[i]] * eps[ti[j]])
            sij = 0.5*(sig[ti[i]] + sig[ti[j]])
            x = sij / r
            # LJ-ish
            e = 4*eij*(x**12 - x**6)
            E += e
            # force
            dEdr = 4*eij*(-12*x**12/r + 6*x**6/r)
            fij = -dEdr * (dr / r)
            F[i] += fij
            F[j] -= fij
    return float(E), F

# 30フレーム生成（小さく）
n_frames = 30
with OUT.open("w") as f:
    types = []
    for t in type_names:
        types += [t]*counts[t]
    for k in range(n_frames):
        pos = pos0 + 0.15*rng.standard_normal(pos0.shape)  # 少し揺らす
        E, F = energy_forces(pos)
        f.write(f"{N}\n")
        # extxyzヘッダ
        f.write(f"Lattice=\"{cell[0,0]:.6f} 0 0  0 {cell[1,1]:.6f} 0  0 0 {cell[2,2]:.6f}\" ")
        f.write("Properties=species:S:1:pos:R:3:force:R:3 ")
        f.write(f"energy={E:.12f} pbc=\"T T T\"\n")
        for i in range(N):
            x,y,z = pos[i]
            fx,fy,fz = F[i]
            f.write(f"{types[i]} {x:.8f} {y:.8f} {z:.8f} {fx:.8f} {fy:.8f} {fz:.8f}\n")

print(f"Wrote {OUT} (synthetic; replace with QE-derived data later)")
