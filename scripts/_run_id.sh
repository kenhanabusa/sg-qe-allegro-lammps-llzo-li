#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 1) RUN_ID 環境変数があればそれを返す
if [ -n "${RUN_ID:-}" ]; then
  echo "$RUN_ID"
  exit 0
fi

# 2) 既に学習結果があるなら最新の run ディレクトリ名を返す
if ls -1d "$ROOT/allegro/runs/"* >/dev/null 2>&1; then
  ls -1dt "$ROOT/allegro/runs/"* | head -n1 | xargs -n1 basename
  exit 0
fi

# 3) 何も無ければ新規作成
date +%Y%m%d-%H%M%S
