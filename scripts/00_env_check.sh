#!/usr/bin/env bash
set -euo pipefail
echo "[env] hostname=$(hostname)"
echo "[env] pwd=$(pwd)"
command -v nvidia-smi >/dev/null && nvidia-smi -L | head || true
command -v apptainer >/dev/null && apptainer --version || true
