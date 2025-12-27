#!/usr/bin/env python3
import sys
import numpy as np
from ase import Atoms
from ase.io import write
from ase.calculators.singlepoint import SinglePointCalculator

def main(out_path: str):
    rng = np.random.default_rng(0)

    # Li/LLZOっぽい最小構成（配管テスト用）
    comp = [("Li", 16), ("La", 8), ("Zr", 8), ("O", 32)]
    symbols = [el for el, n in comp for _ in range(n)]
    n = len(symbols)

    cell = np.diag([12.0, 12.0, 12.0])  # Å
    diag = np.diag(cell)

    frames = []
    nframes = 20
    for i in range(nframes):
        pos = rng.random((n, 3)) * diag
        a = Atoms(symbols=symbols, positions=pos, cell=cell, pbc=True)

        # ★ここが重要：calc 経由で energy/forces を持たせる（extxyzに確実に出る）
        F = rng.normal(scale=0.05, size=(n, 3))     # eV/Å っぽい雰囲気の値
        E = float(rng.normal(loc=-1.0, scale=0.2))  # eV っぽい雰囲気の値
        a.calc = SinglePointCalculator(a, energy=E, forces=F)

        # 念のため total_energy も info に残す（あなたのpatch側が total_energy を参照しててもOK）
        a.info["total_energy"] = E

        frames.append(a)

    write(out_path, frames, format="extxyz", write_results=True)
    print(f"Wrote {out_path} (frames={nframes}, atoms={n})")

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "dataset/raw/synth_llzo_li.extxyz"
    main(out)
