#!/usr/bin/env python3
"""Apply initial SNP/DP/MAC filters and cache ingroup genotypes."""

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-dp", type=int, default=15)
    parser.add_argument("--max-dp", type=int, default=150)
    parser.add_argument("--min-mac", type=int, default=4)
    parser.add_argument("--initial-capacity", type=int, default=200000)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(args.metadata, encoding="utf-8-sig")

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                break
        else:
            raise RuntimeError("VCF header line not found")

    order = pd.DataFrame({"sampleID": samples, "vcf_index": range(len(samples))})
    md = order.merge(metadata, on="sampleID", how="left", validate="one_to_one")
    required = {"s384", "s380", "Plate", "geoID", "spID"}
    missing = required.difference(md.columns)
    if missing:
        raise RuntimeError(f"Missing or empty metadata fields: {sorted(missing)}")
    if md["s384"].isna().any():
        unmatched = md.loc[md["s384"].isna(), "sampleID"].tolist()
        raise RuntimeError(f"VCF samples absent from metadata: {unmatched}")
    md.to_csv(outdir / "s384_vcf_metadata.csv", index=False)

    ingroup_idx = np.flatnonzero(md["s380"].notna().to_numpy())
    if len(ingroup_idx) != 380:
        raise RuntimeError(f"Expected 380 ingroup samples; found {len(ingroup_idx)}")
    n_ing = len(ingroup_idx)
    min_called10 = int(np.ceil(0.10 * n_ing))
    min_called80 = int(np.ceil(0.80 * n_ing))

    capacity = args.initial_capacity
    genotypes = np.lib.format.open_memmap(
        outdir / "genotypes_call10.int8.npy",
        mode="w+",
        dtype=np.int8,
        shape=(capacity, n_ing),
    )
    genotypes[:] = -1

    called10 = np.zeros(n_ing, dtype=np.int64)
    het10 = np.zeros(n_ing, dtype=np.int64)
    depth10 = np.zeros(n_ing, dtype=np.int64)
    called80 = np.zeros(n_ing, dtype=np.int64)
    het80 = np.zeros(n_ing, dtype=np.int64)
    depth80 = np.zeros(n_ing, dtype=np.int64)
    sites = []
    counts = {
        "total_records": 0,
        "biallelic_snp_records": 0,
        "pass_call10_mac4_dp15_150": 0,
        "pass_call80_mac4_dp15_150": 0,
    }

    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            counts["total_records"] += 1
            fields = line.rstrip().split("\t")
            chrom, pos, variant_id, ref, alt = fields[:5]
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            counts["biallelic_snp_records"] += 1
            fmt = fields[8].split(":")
            gt_i, dp_i = fmt.index("GT"), fmt.index("DP")

            genotype = np.full(n_ing, -1, dtype=np.int8)
            depths = np.zeros(n_ing, dtype=np.int16)
            for out_i, vcf_i in enumerate(ingroup_idx):
                values = fields[9 + vcf_i].split(":")
                try:
                    gt, dp = values[gt_i], int(values[dp_i])
                except (ValueError, IndexError):
                    continue
                if not args.min_dp <= dp <= args.max_dp or gt in {"./.", ".|.", "."}:
                    continue
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) != 2 or any(x not in {"0", "1"} for x in alleles):
                    continue
                genotype[out_i] = int(alleles[0]) + int(alleles[1])
                depths[out_i] = dp

            called = genotype >= 0
            n_called = int(called.sum())
            if n_called < min_called10:
                continue
            ac = int(genotype[called].sum())
            mac = min(ac, 2 * n_called - ac)
            if mac < args.min_mac:
                continue

            row = counts["pass_call10_mac4_dp15_150"]
            if row >= capacity:
                raise RuntimeError("Increase --initial-capacity")
            genotypes[row] = genotype
            pass80 = n_called >= min_called80
            sites.append(
                (chrom, int(pos), variant_id, ref, alt, n_called,
                 n_called / n_ing, ac, mac, int(pass80))
            )
            counts["pass_call10_mac4_dp15_150"] += 1
            called10 += called
            het10 += genotype == 1
            depth10 += depths
            if pass80:
                counts["pass_call80_mac4_dp15_150"] += 1
                called80 += called
                het80 += genotype == 1
                depth80 += depths

    genotypes.flush()
    sites_df = pd.DataFrame(
        sites,
        columns=["CHROM", "POS", "ID", "REF", "ALT", "N_CALLED",
                 "CALL_RATE", "AC", "MAC", "PASS_CALL80"],
    )
    sites_df.to_csv(outdir / "sites_call10.tsv.gz", sep="\t", index=False)

    sample_qc = md.iloc[ingroup_idx].copy().reset_index(drop=True)
    sample_qc["called_call10"] = called10
    sample_qc["missing_rate_call10"] = 1 - called10 / len(sites_df)
    sample_qc["heterozygosity_call10"] = het10 / np.maximum(called10, 1)
    sample_qc["mean_depth_call10"] = depth10 / np.maximum(called10, 1)
    n80 = counts["pass_call80_mac4_dp15_150"]
    sample_qc["called_call80"] = called80
    sample_qc["missing_rate_call80"] = 1 - called80 / n80
    sample_qc["heterozygosity_call80"] = het80 / np.maximum(called80, 1)
    sample_qc["mean_depth_call80"] = depth80 / np.maximum(called80, 1)
    sample_qc.to_csv(outdir / "sample_qc.csv", index=False)

    counts.update(
        n_vcf_samples=len(samples),
        n_ingroup_samples=n_ing,
        min_called_call10=min_called10,
        min_called_call80=min_called80,
        genotype_matrix_capacity=capacity,
        genotype_matrix_used_rows=len(sites_df),
    )
    with open(outdir / "pass1_summary.json", "w") as handle:
        json.dump(counts, handle, indent=2)
    print(json.dumps(counts, indent=2))


if __name__ == "__main__":
    main()
