#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

python3 "$REPO_DIR/scripts/10_prepare_trna_domain_bed.py" \
  --header \
  --out "$REPO_DIR/hg38-tRNA-domains.bed"
