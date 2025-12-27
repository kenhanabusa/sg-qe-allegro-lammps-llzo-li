from pathlib import Path
import sys
import numpy as np

extxyz = Path(sys.argv[1])
out = Path(sys.argv[2])

# 超簡易extxyz parser（最初のフレームだけ）
lines = extxyz.read_text().splitlines()
n = int(lines[0].strip())
hdr = lines[1]
atoms = lines[2:2+n]

# Lattice抽出（超簡易）
import re
m = re.search(r'Lattice="([^"]+)"', hdr)
if not m:
  raise SystemExit("No Lattice in extxyz header")
lat = list(map(float, m.group(1).split()))
# 3x3 row-major
ax = lat[0:3]; bx = lat[3:6]; cx = lat[6:9]
lx = ax[0]; ly = bx[1]; lz = cx[2]

species = []
pos = []
for a in atoms:
  sp, x,y,z, *_ = a.split()
  species.append(sp)
  pos.append([float(x),float(y),float(z)])
pos = np.array(pos)

types = {"Li":1,"La":2,"Zr":3,"O":4}
t = [types[s] for s in species]

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w") as f:
  f.write("LAMMPS data (from extxyz first frame)\n\n")
  f.write(f"{n} atoms\n")
  f.write(f"{len(types)} atom types\n\n")
  f.write(f"0.0 {lx:.6f} xlo xhi\n")
  f.write(f"0.0 {ly:.6f} ylo yhi\n")
  f.write(f"0.0 {lz:.6f} zlo zhi\n\n")
  f.write("Masses\n\n")
  f.write("1 6.94\n2 138.905\n3 91.224\n4 15.999\n\n")
  f.write("Atoms # atomic\n\n")
  for i,(ti,pi) in enumerate(zip(t,pos), start=1):
    f.write(f"{i} {ti} {pi[0]:.8f} {pi[1]:.8f} {pi[2]:.8f}\n")
print(f"Wrote {out}")
