#!/usr/bin/env python3
"""
build_mspi_fragments.py

Reconstructs an in-silico MspI genome digest matching the structure of your
original MspI_fragments_hg38.bed. I don't have your original script — this
was reverse-engineered from the RESULT file you uploaded, by inspecting its
cut pattern against a synthetic test sequence. Please sanity-check the output
against your original file before relying on it (see "Verifying" below).

What MspI does: MspI recognises CCGG and cuts between the two Cs (C^CGG).
For each chromosome, this script finds every CCGG occurrence (both strands
are equivalent since CCGG is palindromic) and cuts the chromosome into
fragments at each site, using the chromosome start/end as the outer
boundaries.


Requirements:
    pip install pyfaidx --break-system-packages

Usage:
    python build_mspi_fragments.py \\
        --fasta /path/to/hg38.fa \\
        --out MspI_fragments_hg38.bed

Runtime: whole-genome digest of hg38 takes a few minutes and ~1-2 GB RAM
with pyfaidx's indexed access (it doesn't load the whole genome into memory
at once).
"""

import argparse
import re

STANDARD_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def digest_chromosome(seq: str):
    """Return a list of (start, end) 0-based half-open fragment coordinates
    for one chromosome sequence, cut at every CCGG site (C^CGG)."""
    seq = str(seq).upper()

    # (?=CCGG) is a zero-width lookahead so overlapping CCGG occurrences
    # (e.g. within CCGGCCGG) are still all detected — re.finditer with a
    # plain "CCGG" pattern would skip overlapping matches.
    cut_positions = {0, len(seq)}
    for m in re.finditer(r"(?=CCGG)", seq):
        cut_positions.add(m.start() + 1)  # cut is between the 1st and 2nd C

    cuts = sorted(cut_positions)
    return list(zip(cuts[:-1], cuts[1:]))


def main():
    parser = argparse.ArgumentParser(description="In-silico MspI genome digest")
    parser.add_argument("--fasta", required=True, help="Path to hg38 reference FASTA (indexed or not; pyfaidx will index it)")
    parser.add_argument("--out", required=True, help="Output BED path")
    parser.add_argument("--chroms", nargs="*", default=STANDARD_CHROMS,
                         help="Chromosomes to include (default: chr1-22, X, Y, M — matches the uploaded result file)")
    args = parser.parse_args()

    try:
        from pyfaidx import Fasta
    except ImportError:
        raise SystemExit("pyfaidx is required: pip install pyfaidx --break-system-packages")

    genome = Fasta(args.fasta)

    with open(args.out, "w") as out:
        out.write("chr\tfrag_start\tfrag_end\tfrag_id\n")
        for chrom in args.chroms:
            if chrom not in genome:
                print(f"WARNING: {chrom} not found in FASTA, skipping")
                continue
            seq = genome[chrom][:].seq
            fragments = digest_chromosome(seq)
            for start, end in fragments:
                out.write(f"{chrom}\t{start}\t{end}\t{chrom}:{start}-{end}\n")
            print(f"{chrom}: {len(fragments)} fragments")

    print(f"\nDone. Wrote fragments to {args.out}")


if __name__ == "__main__":
    main()
