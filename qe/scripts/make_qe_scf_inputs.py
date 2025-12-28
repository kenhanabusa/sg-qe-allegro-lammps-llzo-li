#!/usr/bin/env python3
import os, sys
import numpy as np
from pathlib import Path
from ase.io import read, write
from ase.io.espresso import write_espresso_in

def pick_pp(pseudo_dir: Path, element: str) -> str:
    pats = [f"{element}*.UPF", f"{element}*.upf"]
    for pat in pats:
        hits = sorted(pseudo_dir.glob(pat))
        if hits:
            return hits[0].name
    raise SystemExit(f"Cannot find pseudopotential for {element} in {pseudo_dir}")

def main():
    cif = Path(sys.argv[1])
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    disp = float(sys.argv[3]) if len(sys.argv) > 3 else 0.02

    out_snap = Path("qe/snapshots"); out_in = Path("qe/inputs")
    out_snap.mkdir(parents=True, exist_ok=True); out_in.mkdir(parents=True, exist_ok=True)

    pseudo_dir = Path(os.environ["PSEUDO_DIR"]).expanduser()
    elems = ["Li","La","Zr","O"]
    pseudos = {e: pick_pp(pseudo_dir, e) for e in elems}

    ecutwfc = float(os.environ.get("ECUTWFC", "60"))
    ecutrho = float(os.environ.get("ECUTRHO", "480"))
    degauss = float(os.environ.get("DEGAUSS", "0.01"))
    conv_thr = float(os.environ.get("CONV_THR", "1e-6"))
    mixing_beta = float(os.environ.get("MIXING_BETA", "0.3"))

    base = read(cif)
    kpts = (1,1,1)

    rng = np.random.default_rng(0)
    for i in range(n):
        a = base.copy()
        a.positions = a.positions + disp * rng.standard_normal(a.positions.shape)

        snap = out_snap / f"snap_{i:03d}.cif"
        write(snap, a)

        prefix = f"snap_{i:03d}"
        input_data = dict(
            control=dict(
                calculation="scf",
                prefix=prefix,
                pseudo_dir=str(pseudo_dir),
                outdir="./qe/out",
                tprnfor=True,
                tstress=False,
                verbosity="high",
            ),
            system=dict(
                ibrav=0,
                ecutwfc=ecutwfc,
                ecutrho=ecutrho,
                occupations="smearing",
                smearing="mv",
                degauss=degauss,
                nosym=True,
                noinv=True,
            ),
            electrons=dict(
                conv_thr=conv_thr,
                mixing_beta=mixing_beta,
            ),
        )

        inp = out_in / f"{prefix}.in"
        with inp.open("w") as f:
            write_espresso_in(
                f, a,
                input_data=input_data,
                pseudopotentials=pseudos,
                kpts=kpts,
            )

    print(f"OK: wrote {n} inputs into qe/inputs (disp={disp} Å)")
    print("pseudos:", pseudos)

if __name__ == "__main__":
    main()
