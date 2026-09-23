#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 ALIGNMENT.phy OUTPUT_PREFIX" >&2
  exit 1
fi

alignment="$(realpath "$1")"
output_prefix="$2"
threads="${THREADS:-10}"
seed="${IQTREE_SEED:-20260923}"
mkdir -p "$(dirname "$output_prefix")"

if command -v iqtree3 >/dev/null 2>&1; then
  iqtree_exe="iqtree3"
elif command -v iqtree2 >/dev/null 2>&1; then
  iqtree_exe="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  iqtree_exe="iqtree"
else
  echo "IQ-TREE executable not found (tried iqtree3, iqtree2, and iqtree)" >&2
  exit 1
fi

"$iqtree_exe" \
  -s "$alignment" \
  -st DNA \
  -m GTR+ASC \
  -B 1000 \
  --alrt 1000 \
  -T "$threads" \
  -seed "$seed" \
  --prefix "$output_prefix"

echo "IQ-TREE result: ${output_prefix}.treefile"
