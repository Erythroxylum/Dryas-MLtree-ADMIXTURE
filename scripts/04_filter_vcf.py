#!/usr/bin/env python3
"""Write the cleaned 384-sample tree and 380-sample ADMIXTURE VCFs."""

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


BASE_INDEX = {"C": 0, "A": 1, "T": 2, "G": 3}  # ipyrad CATG order
HETEROZYGOTES = {"0/1", "1/0", "0|1", "1|0"}


def mask_genotype(value, fmt, ref, alt, min_dp, max_dp, min_balance):
    fields = value.split(":")
    gt_i, dp_i = fmt.index("GT"), fmt.index("DP")
    gt = fields[gt_i]
    try:
        depth = int(fields[dp_i])
    except (ValueError, IndexError):
        depth = 0
    valid = gt not in {"./.", ".|.", "."} and min_dp <= depth <= max_dp

    if valid and gt in HETEROZYGOTES:
        try:
            counts = [int(x) for x in fields[fmt.index("CATG")].split(",")]
            ref_count = counts[BASE_INDEX[ref]]
            alt_count = counts[BASE_INDEX[alt]]
            balance = alt_count / (ref_count + alt_count)
            valid = min_balance <= balance <= 1 - min_balance
        except (ValueError, IndexError, ZeroDivisionError, KeyError):
            valid = False

    alleles = gt.replace("|", "/").split("/")
    if not valid or len(alleles) != 2 or any(x not in {"0", "1"} for x in alleles):
        fields[gt_i] = "./."
        return ":".join(fields), None, None
    return ":".join(fields), int(alleles[0]) + int(alleles[1]), depth


def update_info(info, n_called, depth):
    values, seen = [], set()
    for item in info.split(";"):
        key = item.split("=", 1)[0]
        if key == "NS":
            values.append(f"NS={n_called}")
            seen.add("NS")
        elif key == "DP":
            values.append(f"DP={depth}")
            seen.add("DP")
        else:
            values.append(item)
    if "NS" not in seen:
        values.append(f"NS={n_called}")
    if "DP" not in seen:
        values.append(f"DP={depth}")
    return ";".join(values)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--qc-dir", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-dp", type=int, default=15)
    parser.add_argument("--max-dp", type=int, default=150)
    parser.add_argument("--min-balance", type=float, default=0.25)
    parser.add_argument("--min-mac", type=int, default=4)
    args = parser.parse_args()

    qc_dir = Path(args.qc_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    metadata = pd.read_csv(qc_dir / "s384_vcf_metadata.csv")
    ingroup = [i for i, value in enumerate(metadata["s380"].notna()) if value]
    excluded = set(
        pd.read_csv(
            qc_dir / "exclude_loci_strict_FDR05.txt", header=None
        )[0].astype(str)
    )

    configurations = [
        {
            "name": "s384_batchclean_ab25_locusFDR05_mac4_dp15-150_miss10.vcf.gz",
            "subset": False,
            "min_called": 38,
        },
        {
            "name": "s380_batchclean_ab25_locusFDR05_mac4_dp15-150_miss80.vcf.gz",
            "subset": True,
            "min_called": 304,
        },
    ]
    outputs = []
    counts = {}
    for config in configurations:
        handle = gzip.open(outdir / config["name"], "wt", compresslevel=6)
        outputs.append((handle, config))
        counts[config["name"]] = 0

    with gzip.open(args.vcf, "rt") as source:
        for line in source:
            if line.startswith("##"):
                for handle, _ in outputs:
                    handle.write(line)
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                for handle, config in outputs:
                    call_rate = 0.8 if config["subset"] else 0.1
                    handle.write(
                        "##batch_cleaning=<"
                        f"excluded_ipyrad_loci={len(excluded)},"
                        f"heterozygote_alt_balance={args.min_balance}-"
                        f"{1 - args.min_balance},DP={args.min_dp}-{args.max_dp},"
                        f"MAC={args.min_mac},ingroup_call_rate={call_rate}>\n"
                    )
                    samples = [header[9 + i] for i in ingroup] if config["subset"] else header[9:]
                    handle.write("\t".join(header[:9] + samples) + "\n")
                continue

            fields = line.rstrip().split("\t")
            ref, alt = fields[3], fields[4]
            if (
                "," in alt or len(ref) != 1 or len(alt) != 1
                or ref not in BASE_INDEX or alt not in BASE_INDEX
            ):
                continue
            match = re.match(r"^(loc\d+)", fields[2])
            locus = match.group(1) if match else None
            if locus in excluded:
                continue

            fmt = fields[8].split(":")
            masked, genotype, depth = [], [], []
            for value in fields[9:]:
                new_value, dosage, dp = mask_genotype(
                    value, fmt, ref, alt,
                    args.min_dp, args.max_dp, args.min_balance,
                )
                masked.append(new_value)
                genotype.append(dosage)
                depth.append(dp)

            called_ingroup = [i for i in ingroup if genotype[i] is not None]
            n_called = len(called_ingroup)
            ac = sum(genotype[i] for i in called_ingroup)
            mac = min(ac, 2 * n_called - ac) if n_called else 0
            if mac < args.min_mac:
                continue

            for handle, config in outputs:
                if n_called < config["min_called"]:
                    continue
                chosen = ingroup if config["subset"] else range(len(masked))
                chosen_called = [i for i in chosen if genotype[i] is not None]
                out_fields = fields[:9]
                out_fields[7] = update_info(
                    out_fields[7],
                    len(chosen_called),
                    sum(depth[i] for i in chosen_called),
                )
                handle.write("\t".join(out_fields + [masked[i] for i in chosen]) + "\n")
                counts[config["name"]] += 1

    for handle, _ in outputs:
        handle.close()
    pd.DataFrame(
        [{"file": name, "n_sites": count} for name, count in counts.items()]
    ).to_csv(outdir / "cleaned_vcf_counts.csv", index=False)
    for name, count in counts.items():
        print(f"{name}: {count:,} SNPs")


if __name__ == "__main__":
    main()
