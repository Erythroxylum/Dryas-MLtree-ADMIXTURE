#!/usr/bin/env python3
"""Convert a biallelic SNP VCF into a relaxed PHYLIP IUPAC alignment."""

import argparse
import gzip

import numpy as np


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
    parser.add_argument(
        "--keep-constant-compatible",
        action="store_true",
        help=(
            "Retain columns with only one unambiguous nucleotide after IUPAC "
            "encoding. By default these are removed for compatibility with +ASC."
        ),
    )
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
    retained_sites = n_sites
    if not args.keep_constant_compatible:
        alignment = np.frombuffer(b"".join(sequences), dtype="S1").reshape(
            len(samples), n_sites
        )
        observed_states = np.zeros(n_sites, dtype=np.uint8)
        for bit, base in enumerate((b"A", b"C", b"G", b"T")):
            observed_states |= (alignment == base).any(axis=0).astype(np.uint8) << bit
        state_counts = np.fromiter(
            (int(value).bit_count() for value in observed_states),
            dtype=np.uint8,
            count=n_sites,
        )
        keep = state_counts >= 2
        retained_sites = int(keep.sum())
        if retained_sites == 0:
            raise RuntimeError("No unambiguously variable sites remain")
        sequences = [
            np.frombuffer(sequence, dtype=np.uint8)[keep].tobytes()
            for sequence in sequences
        ]

    with open(args.output, "w") as output:
        output.write(f"{len(samples)} {retained_sites}\n")
        for sample, sequence in zip(samples, sequences):
            output.write(f"{sample} {bytes(sequence).decode('ascii')}\n")
    print(
        f"Wrote {len(samples)} samples and {retained_sites:,} SNPs to {args.output}"
    )
    if retained_sites != n_sites:
        print(
            f"Removed {n_sites - retained_sites:,} constant-compatible columns "
            "after IUPAC encoding for +ASC"
        )


if __name__ == "__main__":
    main()
