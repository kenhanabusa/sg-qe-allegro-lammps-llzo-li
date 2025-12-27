#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RID_FILE="$ROOT/.run_id"
if [ -n "${RUN_ID:-}" ]; then
  echo "$RUN_ID" > "$RID_FILE"
fi
if [ -f "$RID_FILE" ]; then
  RUN_ID="$(cat "$RID_FILE")"
else
  RUN_ID="$(date +%Y%m%d-%H%M%S)"
  echo "$RUN_ID" > "$RID_FILE"
fi
export RUN_ID
