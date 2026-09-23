#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 CLEANED_380.vcf.gz OUTPUT_DIR" >&2
  exit 1
fi

input_vcf="$(realpath "$1")"
output_dir="$(mkdir -p "$2" && realpath "$2")"
threads="${THREADS:-10}"
k_min="${K_MIN:-2}"
k_max="${K_MAX:-19}"
replicates="${REPLICATES:-2}"
base_seed="${BASE_SEED:-20260923}"
work_dir="${output_dir}/prepared"
prefix="${work_dir}/dryas_cleaned"
mkdir -p "$work_dir"

# Convert ordinary gzip to BGZF for indexed access.
gzip -cd "$input_vcf" | bgzip -c > "${prefix}.vcf.gz"
tabix -f -p vcf "${prefix}.vcf.gz"
bcftools query -l "${prefix}.vcf.gz" > "${output_dir}/sample_order.txt"

plink --allow-extra-chr --double-id \
  --vcf "${prefix}.vcf.gz" \
  --make-bed \
  --out "$prefix"

# Physical 50-kb windows, step size 5 variants, r^2 threshold 0.2.
plink --allow-extra-chr \
  --bfile "$prefix" \
  --indep-pairwise 50kb 5 0.2 \
  --out "${prefix}_ld"

plink --allow-extra-chr \
  --bfile "$prefix" \
  --extract "${prefix}_ld.prune.in" \
  --make-bed \
  --out "${prefix}_ld"

# ADMIXTURE treats chromosome code 0 as unlinked/nonstandard sequence data.
awk '{$1="0"; print}' "${prefix}_ld.bim" > "${prefix}_ld.bim.tmp"
mv "${prefix}_ld.bim.tmp" "${prefix}_ld.bim"

for rep in $(seq 1 "$replicates"); do
  rep_dir="${output_dir}/rep${rep}"
  mkdir -p "$rep_dir"
  ln -sf "${prefix}_ld.bed" "${rep_dir}/dryas_ld.bed"
  ln -sf "${prefix}_ld.bim" "${rep_dir}/dryas_ld.bim"
  ln -sf "${prefix}_ld.fam" "${rep_dir}/dryas_ld.fam"
  for k in $(seq "$k_min" "$k_max"); do
    seed=$((base_seed + rep * 1000 + k))
    (
      cd "$rep_dir"
      admixture --cv=10 -s "$seed" -j"$threads" dryas_ld.bed "$k" \
        > "K${k}.rep${rep}.log" 2>&1
      mv "dryas_ld.${k}.Q" "K${k}.rep${rep}.Q"
      mv "dryas_ld.${k}.P" "K${k}.rep${rep}.P"
    )
  done
done

{
  echo -e "K\treplicate\tCV_error"
  for log in "${output_dir}"/rep*/K*.log; do
    k=$(basename "$log" | sed -E 's/K([0-9]+)\.rep[0-9]+\.log/\1/')
    rep=$(basename "$log" | sed -E 's/K[0-9]+\.rep([0-9]+)\.log/\1/')
    cv=$(grep -Eo 'CV error \(K=[0-9]+\): [0-9.eE+-]+' "$log" | awk '{print $NF}' | tail -1)
    echo -e "${k}\t${rep}\t${cv}"
  done
} > "${output_dir}/cross_validation.tsv"

echo "ADMIXTURE results: ${output_dir}"

