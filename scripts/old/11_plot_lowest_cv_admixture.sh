#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 4 || $# -gt 6 ]]; then
  echo "Usage: $0 TREEFILE ADMIXTURE_DIR METADATA.csv OUTPUT.pdf [K_MIN] [K_MAX]" >&2
  exit 1
fi

tree_file="$1"
admixture_dir="$2"
metadata_file="$3"
output_file="$4"
k_min="${5:-16}"
k_max="${6:-22}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PLOT_MODE=lowest_cv \
Rscript "${script_dir}/10_plot_tree_admixture.R" \
  "$tree_file" \
  "$admixture_dir" \
  "$metadata_file" \
  "$output_file" \
  "$k_min" "$k_max"

