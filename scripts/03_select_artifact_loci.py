#!/usr/bin/env python3
"""Select batch-associated SNPs and expand them to complete ipyrad loci."""

import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--qc-dir", required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-het-excess", type=float, default=0.10)
    parser.add_argument("--min-strata", type=int, default=5)
    parser.add_argument("--min-positive-fraction", type=float, default=0.70)
    args = parser.parse_args()
    qc_dir = Path(args.qc_dir)

    metrics = pd.read_csv(qc_dir / "locus_batch_metrics.tsv.gz", sep="\t")
    metrics["LOCUS"] = metrics["ID"].str.extract(r"^(loc\d+)")
    non_tied = (metrics["het_pos"] + metrics["het_neg"]).clip(lower=1)
    keep = (
        (metrics["het_q"] < args.fdr)
        & (metrics["het_rate_diff_P1_minus_other"] > args.min_het_excess)
        & (metrics["het_info"] >= args.min_strata)
        & (metrics["het_pos"] / non_tied >= args.min_positive_fraction)
        & metrics["LOCUS"].notna()
    )
    flagged = metrics.loc[keep].sort_values("het_q")
    loci = sorted(flagged["LOCUS"].unique())
    pd.Series(loci).to_csv(
        qc_dir / "exclude_loci_strict_FDR05.txt", index=False, header=False
    )
    flagged.to_csv(qc_dir / "flagged_artifact_sites.tsv.gz", sep="\t", index=False)

    conservative = keep & (metrics["het_q"] < 0.01)
    loci01 = sorted(metrics.loc[conservative, "LOCUS"].unique())
    pd.Series(loci01).to_csv(
        qc_dir / "exclude_loci_conservative_FDR01.txt", index=False, header=False
    )
    summary = {
        "flagged_snps_fdr05": int(keep.sum()),
        "excluded_loci_fdr05": len(loci),
        "excluded_loci_fdr01": len(loci01),
    }
    with open(qc_dir / "artifact_locus_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
