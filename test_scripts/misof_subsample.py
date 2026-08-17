#!/usr/bin/env python3
"""
Build a random partial NEXUS alignment from misofproteinalignment.nex by
picking individual alignment columns (base pairs) uniformly at random,
without replacement, until a target total length is reached, then writing
just those columns (in ascending column order, so the result is still a
valid alignment) for every taxon.

misofproteinalignment.nex's own "matrix" block is a plain NEXUS block: one
line per taxon, "<name><TAB><full ungapped-length sequence>", terminated by
a line containing just ";". Everything after that (the "begin SETS;" block
of CHARSET/TAXSET definitions) is irrelevant here and is never parsed.

Output is NEXUS, not FASTA, and deliberately keeps the same explicit
"format datatype=protein" declaration the input file has. IQ-TREE's NEXUS
reader (Alignment::extractDataBlock, alignment/alignment.cpp) trusts that
declaration directly; a FASTA file has no such field, so IQ-TREE instead
falls back to content-based sniffing (Alignment::detectSequenceType), which
requires >90% of a sequence's letters to be one of the 20 unambiguous amino
acids. misofproteinalignment.nex is ~32% ambiguous 'X' overall (a real
phylogenomic supermatrix, not simulated data), well past that threshold, so
EVERY random column subsample of it -- not just unlucky draws -- reads as
"Unknown sequence type" as FASTA. Keeping the explicit NEXUS declaration
sidesteps the heuristic entirely, exactly as the input file itself does.

Usage:
    python3 misof_subsample.py <input.nex> <n_sites> <output.nex> [--seed N]
"""
import argparse
import random
import sys


def read_matrix(nexus_path):
    """Return (taxa_in_order, {name: sequence}) from the file's matrix block."""
    taxa = []
    sequences = {}
    in_matrix = False
    with open(nexus_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if not in_matrix:
                if stripped.lower() == "matrix":
                    in_matrix = True
                continue
            if stripped == ";":
                break
            if not stripped:
                continue
            name, seq = stripped.split(None, 1)
            taxa.append(name)
            sequences[name] = seq
    return taxa, sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_nexus")
    parser.add_argument("n_sites", type=int)
    parser.add_argument("output_nexus")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    rng = random.Random(args.seed)

    taxa, sequences = read_matrix(args.input_nexus)
    if not taxa:
        sys.exit(f"error: no taxa found in matrix block of {args.input_nexus}")

    lengths = {len(seq) for seq in sequences.values()}
    if len(lengths) != 1:
        sys.exit(f"error: taxa have inconsistent sequence lengths: {sorted(lengths)}")
    nchar = lengths.pop()

    if args.n_sites > nchar:
        sys.exit(f"error: requested {args.n_sites} sites but alignment only has {nchar}")

    columns = sorted(rng.sample(range(nchar), args.n_sites))

    with open(args.output_nexus, "w") as out:
        out.write("#NEXUS\n")
        out.write("begin DATA;\n")
        out.write(f"\tdimensions ntax={len(taxa)} nchar={len(columns)};\n")
        out.write("\tformat datatype=protein missing=? gap=-;\n")
        out.write("\tmatrix\n")
        for name in taxa:
            seq = sequences[name]
            partial = "".join(seq[i] for i in columns)
            out.write(f"\t{name}\t{partial}\n")
        out.write("\t;\n")
        out.write("end;\n")

    print(f"wrote {len(taxa)} taxa x {len(columns)} sites to {args.output_nexus}", file=sys.stderr)


if __name__ == "__main__":
    main()
