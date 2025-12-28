#!/usr/bin/env python3
import re, sys
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator

RY_TO_EV = 13.605693122994
BOHR_TO_ANG = 0.529177210903
RY_BOHR_TO_EV_ANG = RY_TO_EV / BOHR_TO_ANG

def parse_qe_energy_forces(out_text: str, natoms: int):
    m = re.findall(r'!\s+total energy\s+=\s+([-\d\.Ee+]+)\s+Ry', out_text)
    if not m:
        raise ValueError("No '! total energy = ... Ry' found")
    e_ry = float(m[-1])

    idx = out_text.rfind("Forces acting on atoms")
    if idx < 0:
        raise ValueError("No 'Forces acting on atoms' block found")
    tail = out_text[idx:].splitlines()
    forces = []
    for line in tail:
        mm = re.search(r'force\s*=\s*([-\d\.Ee+]+)\s+([-\d\.Ee+]+)\s+([-\d\.Ee+]+)', line)
        if mm:
            forces.append([float(mm.group(1)), float(mm.group(2)), float(mm.group(3))])
            if len(forces) == natoms:
                break
    if len(forces) != natoms:
        raise ValueError(f"Forces parsed {len(forces)} != natoms {natoms}")

    forces = np.array(forces) * RY_BOHR_TO_EV_ANG
    e_ev = e_ry * RY_TO_EV
    return e_ev, forces

def main():
    snaps = Path("qe/snapshots")
    outs = Path("qe/outputs")
    out_extxyz = Path("dataset/raw/qe_llzo_li.extxyz")

    files = sorted(snaps.glob("snap_*.cif"))
    if not files:
        raise SystemExit("No snapshots found under qe/snapshots/")
    frames = []
    for snap in files:
        base = snap.stem
        out = outs / f"{base}.out"
        if not out.exists():
            raise SystemExit(f"Missing QE output: {out}")
        a = read(snap)
        txt = out.read_text(errors="ignore")
        e, f = parse_qe_energy_forces(txt, len(a))
        a.calc = SinglePointCalculator(a, energy=e, forces=f)
        a.info["energy"] = e
        a.info["total_energy"] = e
        frames.append(a)

    out_extxyz.parent.mkdir(parents=True, exist_ok=True)
    # avoid ASE extxyz KeyError: keys from calculator may already exist in atoms.info/arrays
    for a in frames:
        # these are commonly written from calculator results
        for k in ('energy','free_energy','stress','virial'):
            a.info.pop(k, None)
        # if any were stored as arrays too, remove to avoid collisions
        for k in ('forces','stress','virial'):
            if k in a.arrays:
                del a.arrays[k]

    write(out_extxyz, frames, format="extxyz", write_results=True)
    print(f"OK: wrote {out_extxyz} (frames={len(frames)})")

if __name__ == "__main__":
    main()
