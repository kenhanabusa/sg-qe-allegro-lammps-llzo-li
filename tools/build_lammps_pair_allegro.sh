#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

# LAMMPS + pair_nequip_allegro をビルドして ~/tools/lammps に入れる
PREFIX="$HOME/tools/lammps"
SRC_LAMMPS="$HOME/tools/src_lammps"
SRC_PAIR="$HOME/tools/src_pair_nequip_allegro"
BUILD="$HOME/tools/build_lammps_pair_allegro"

mkdir -p "$HOME/tools"

# pair_nequip_allegro は 2025-09-10 以降のLAMMPSが必要、と明記  [oai_citation:5‡GitHub](https://github.com/mir-group/pair_nequip_allegro)
rm -rf "$SRC_LAMMPS" "$SRC_PAIR" "$BUILD"
git clone --depth=1 https://github.com/lammps/lammps "$SRC_LAMMPS"
git clone --depth=1 https://github.com/mir-group/pair_nequip_allegro "$SRC_PAIR"

# patch
cd "$SRC_PAIR"
./patch_lammps.sh "$SRC_LAMMPS"

# torchのcmake prefix（README推奨） [oai_citation:6‡GitHub](https://github.com/mir-group/pair_nequip_allegro)
PY=python3
CMAKE_PREFIX_PATH="$($PY -c 'import torch;print(torch.utils.cmake_prefix_path)' 2>/dev/null || true)"
if [ -z "$CMAKE_PREFIX_PATH" ]; then
  echo "[build] ERROR: torch not available in current python env."
  echo "[build] Activate the env that has PyTorch + nequip-allegro, then rerun:"
  echo "  bash tools/build_lammps_pair_allegro.sh"
  exit 1
fi

mkdir -p "$BUILD"
cd "$BUILD"

cmake "$SRC_LAMMPS/cmake" \
  -D CMAKE_CUDA_COMPILER:FILEPATH="${CUDA_NVCC_EXECUTABLE}" \
  -D CUDAToolkit_ROOT:PATH="${CUDA_HOME}" \
  -D Torch_DIR:PATH="${TORCH_DIR}" \
  -D MKL_INCLUDE_DIR:PATH="${CONDA_PREFIX}/include" \

  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_INSTALL_PREFIX="$PREFIX" \
  -D BUILD_MPI=ON \
  -D PKG_KOKKOS=ON \
  -D Kokkos_ENABLE_CUDA=ON \
  -D Kokkos_ENABLE_OPENMP=ON \
  -D Kokkos_ARCH_AMPERE80=ON \
  -D CMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"

cmake --build . -j
cmake --install .

echo "[build] OK: $PREFIX/bin/lmp"
