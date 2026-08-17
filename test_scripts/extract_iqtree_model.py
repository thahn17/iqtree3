#!/usr/bin/env python3
"""
Turn a real iqtree3 run's fitted model parameters into a literal IQ-TREE
model spec string usable by spr_topology_test's own "model <spec>"
--hillclimb flag (tree/spr_topology_test.cpp), so a search there can start
from ML-estimated parameters an external iqtree3 run already fit -- e.g.
LG+FO's estimated amino acid frequencies -- instead of either spr_topology_
test's own un-fit starting values or having no way to supply them at all.

Reads the "Model of substitution: <name>" and "State frequencies:" /
"pi(X) = <value>" lines from a .iqtree report file (protein: 20 lines, DNA:
4), and prints "<base>+F{f1,f2,...}" to stdout, with the base model name's
own "+F"/"+FO"/"+FC" suffix (if any) stripped and the frequencies in the
exact state order the report printed them in -- alignment/alignment.cpp's
symbols_protein ("ARNDCQEGHILKMFPSTWYV") for protein, A/C/G/T for DNA, both
of which are the same order a "+F{...}" spec expects. This is the same
literal-frequency syntax IQ-TREE itself accepts on its own -m flag; nothing
here is spr_topology_test-specific except which flag the result gets passed
to.

Usage:
    python3 extract_iqtree_model.py <run.iqtree>

Example, feeding the result straight into spr_topology_test:
    MODEL_SPEC="$(python3 test_scripts/extract_iqtree_model.py nni_run.iqtree)"
    build/spr_topology_test --hillclimb real.fa 10 10000 fast quiet notree \\
        starttree nni_run.treefile model "$MODEL_SPEC" record
"""
import argparse
import re
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("iqtree_report")
    args = parser.parse_args()

    with open(args.iqtree_report) as f:
        text = f.read()

    model_match = re.search(r"^Model of substitution:\s*(\S+)", text, re.MULTILINE)
    if not model_match:
        sys.exit(f"error: no 'Model of substitution:' line found in {args.iqtree_report}")
    base_model = re.sub(r"\+F[OC]?$", "", model_match.group(1))

    freqs = re.findall(r"^\s*pi\(([A-Za-z*]+)\)\s*=\s*([0-9.eE+-]+)", text, re.MULTILINE)
    if not freqs:
        sys.exit(f"error: no 'pi(X) = ...' frequency lines found in {args.iqtree_report} "
                "-- was this run's -m string missing '+F'/'+FO'/'+FC'?")

    values = [v for _, v in freqs]
    print(f"{base_model}+F{{{','.join(values)}}}")


if __name__ == "__main__":
    main()
