#!/usr/bin/env python3
"""Test Plate 1 heterozygosity excess within matched taxon-locality strata."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2


def bh_fdr(p_values):
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1].clip(max=1)
    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output


def add_cmh(a, b, c, d, numerator, variance):
    n = a + b + c + d
    valid = n > 1
    expected = np.zeros_like(n, dtype=float)
    var = np.zeros_like(n, dtype=float)
    expected[valid] = (a[valid] + b[valid]) * (a[valid] + c[valid]) / n[valid]
    var[valid] = (
        (a[valid] + b[valid]) * (c[valid] + d[valid])
        * (a[valid] + c[valid]) * (b[valid] + d[valid])
        / (n[valid] ** 2 * (n[valid] - 1))
    )
    numerator += np.where(var > 0, a - expected, 0)
    variance += var


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-dir", required=True)
    parser.add_argument("--chunk-size", type=int, default=5000)
    args = parser.parse_args()
    qc_dir = Path(args.qc_dir)

    with open(qc_dir / "pass1_summary.json") as handle:
        summary = json.load(handle)
    n_sites = summary["genotype_matrix_used_rows"]
    genotypes = np.load(qc_dir / "genotypes_call10.int8.npy", mmap_mode="r")
    sites = pd.read_csv(qc_dir / "sites_call10.tsv.gz", sep="\t")
    samples = pd.read_csv(qc_dir / "sample_qc.csv")

    samples["plate1"] = samples["Plate"].astype(str).eq("1")
    samples["stratum"] = (
        samples["geoID"].fillna("NA") + "|" + samples["spID"].fillna("NA")
    )
    eligible = samples["Plate"].astype(str).isin(["1", "2", "3", "4"])
    strata = []
    for name, group in samples[eligible].groupby("stratum"):
        p1 = group.index[group["plate1"]].to_numpy()
        other = group.index[~group["plate1"]].to_numpy()
        if len(p1) and len(other):
            strata.append((name, p1, other))
    if not strata:
        raise RuntimeError("No matched Plate 1 versus Plate 2-4 strata found")

    with open(qc_dir / "matched_strata.json", "w") as handle:
        json.dump(
            [{"stratum": name, "plate1_indices": p1.tolist(),
              "other_indices": other.tolist()} for name, p1, other in strata],
            handle,
            indent=2,
        )

    numerator = np.zeros(n_sites)
    variance = np.zeros(n_sites)
    diff_sum = np.zeros(n_sites)
    weight_sum = np.zeros(n_sites)
    informative = np.zeros(n_sites, dtype=np.int16)
    positive = np.zeros(n_sites, dtype=np.int16)
    negative = np.zeros(n_sites, dtype=np.int16)

    for start in range(0, n_sites, args.chunk_size):
        stop = min(start + args.chunk_size, n_sites)
        block = np.asarray(genotypes[start:stop], dtype=np.int16)
        size = stop - start
        num = np.zeros(size)
        var = np.zeros(size)
        diffs = np.zeros(size)
        weights = np.zeros(size)
        info = np.zeros(size, dtype=np.int16)
        pos = np.zeros(size, dtype=np.int16)
        neg = np.zeros(size, dtype=np.int16)

        for _, p1_idx, other_idx in strata:
            p1 = block[:, p1_idx]
            other = block[:, other_idx]
            p1_called = (p1 >= 0).sum(axis=1)
            other_called = (other >= 0).sum(axis=1)
            p1_het = (p1 == 1).sum(axis=1)
            other_het = (other == 1).sum(axis=1)
            add_cmh(
                p1_het, p1_called - p1_het,
                other_het, other_called - other_het,
                num, var,
            )
            valid = (p1_called > 0) & (other_called > 0)
            difference = np.zeros(size)
            difference[valid] = (
                p1_het[valid] / p1_called[valid]
                - other_het[valid] / other_called[valid]
            )
            weight = np.zeros(size)
            weight[valid] = (
                p1_called[valid] * other_called[valid]
                / (p1_called[valid] + other_called[valid])
            )
            diffs += weight * difference
            weights += weight
            info += valid
            pos += valid & (difference > 0)
            neg += valid & (difference < 0)

        section = slice(start, stop)
        numerator[section] = num
        variance[section] = var
        diff_sum[section] = diffs
        weight_sum[section] = weights
        informative[section] = info
        positive[section] = pos
        negative[section] = neg
        print(f"Processed {stop:,}/{n_sites:,} SNPs", flush=True)

    z2 = np.divide(
        numerator ** 2, variance,
        out=np.zeros_like(numerator), where=variance > 0,
    )
    result = sites.copy()
    result["het_p"] = chi2.sf(z2, 1)
    result["het_q"] = bh_fdr(result["het_p"])
    result["het_rate_diff_P1_minus_other"] = np.divide(
        diff_sum, weight_sum, out=np.zeros(n_sites), where=weight_sum > 0
    )
    result["het_info"] = informative
    result["het_pos"] = positive
    result["het_neg"] = negative
    result.to_csv(qc_dir / "locus_batch_metrics.tsv.gz", sep="\t", index=False)
    print(f"Matched strata: {len(strata)}")


if __name__ == "__main__":
    main()
