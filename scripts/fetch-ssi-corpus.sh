#!/usr/bin/env bash
# Fetch jyviko's exact-ground-truth surface/surface intersection corpus
# (https://github.com/jyviko/ssi-corpus, MIT OR Apache-2.0) at a pinned commit into
# corpus/ssi/ (gitignored; never commit case data). Idempotent: re-running leaves an
# existing checkout at the pinned commit.
#
# Usage: scripts/fetch-ssi-corpus.sh [dest]      (default <repo>/corpus/ssi)
# Used by: cargo run --profile difftest -p difftest --features cpp -- --ssi-corpus
set -euo pipefail

URL="https://github.com/jyviko/ssi-corpus.git"
COMMIT="b04c8a6"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
DEST="${1:-$ROOT/corpus/ssi}"

if [ -d "$DEST/.git" ]; then
  have="$(git -C "$DEST" rev-parse HEAD)"
  case "$have" in
    "$COMMIT"*) echo "ssi-corpus already at $COMMIT in $DEST"; exit 0 ;;
  esac
  git -C "$DEST" fetch --quiet origin
else
  mkdir -p "$(dirname "$DEST")"
  git clone --quiet "$URL" "$DEST"
fi
git -C "$DEST" -c advice.detachedHead=false checkout --quiet "$COMMIT"
echo "ssi-corpus at $(git -C "$DEST" rev-parse --short=7 HEAD) in $DEST"
