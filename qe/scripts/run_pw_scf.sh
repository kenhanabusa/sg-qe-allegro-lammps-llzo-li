#!/usr/bin/env bash
set -euo pipefail

PW="${PW:-pw.x}"
NP="${NP:-1}"

OUTDIR="qe/outputs"
mkdir -p "$OUTDIR" qe/out

# pw.x がリンクしている OpenMPI prefix を自動検出（LAMMPS wrapper と同じ作戦）
LIBMPI="$(ldd "$PW" | awk '/libmpi\.so/ {print $3; exit}')"
if [ -n "${LIBMPI:-}" ]; then
  OMPI_BASE="$(dirname "$(dirname "$LIBMPI")")"
  export PATH="$OMPI_BASE/bin:$PATH"
  export LD_LIBRARY_PATH="$OMPI_BASE/lib:$OMPI_BASE/lib64:${LD_LIBRARY_PATH:-}"
  export OPAL_PREFIX="$OMPI_BASE"
  MPIRUN=( mpirun --prefix "$OMPI_BASE" -np "$NP" )
else
  MPIRUN=( mpirun -np "$NP" )
fi

# NP>1 で HCOLL/ML がIBを要求して落ちる環境があるので、その時だけ除外（NP=1なら不要）
if [ "$NP" -gt 1 ]; then
  MPIRUN+=( --mca coll ^hcoll,ml )
fi

for inpf in qe/inputs/snap_*.in; do
  base="$(basename "$inpf" .in)"
  outf="$OUTDIR/${base}.out"
  echo "[run] $inpf -> $outf"
  "${MPIRUN[@]}" "$PW" -in "$inpf" > "$outf" 2>&1
done

echo "OK: SCF done"
