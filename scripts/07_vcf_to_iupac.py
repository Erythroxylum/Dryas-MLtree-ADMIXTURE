#!/usr/bin/env python3
"""Convert a biallelic SNP VCF into a relaxed PHYLIP IUPAC alignment."""

import argparse
import gzip


IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    samples = None
    sequences = None
    n_sites = 0
    with open_text(args.vcf) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if any(any(char.isspace() for char in name) for name in samples):
                    raise RuntimeError("PHYLIP sample names cannot contain whitespace")
                sequences = [bytearray() for _ in samples]
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            ref, alt = fields[3].upper(), fields[4].upper()
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            fmt = fields[8].split(":")
            gt_i = fmt.index("GT")
            heterozygote = IUPAC.get(frozenset((ref, alt)), "N")
            for sequence, value in zip(sequences, fields[9:]):
                try:
                    gt = value.split(":")[gt_i].replace("|", "/")
                except IndexError:
                    gt = "./."
                if gt == "0/0":
                    base = ref
                elif gt == "1/1":
                    base = alt
                elif gt in {"0/1", "1/0"}:
                    base = heterozygote
                else:
                    base = "N"
                sequence.extend(base.encode("ascii"))
            n_sites += 1

    if not samples or n_sites == 0:
        raise RuntimeError("No samples or biallelic SNPs found")
    with open(args.output, "w") as output:
        output.write(f"{len(samples)} {n_sites}\n")
        for sample, sequence in zip(samples, sequences):
            output.write(f"{sample} {sequence.decode('ascii')}\n")
    print(f"Wrote {len(samples)} samples and {n_sites:,} SNPs to {args.output}")


if __name__ == "__main__":
    main()

