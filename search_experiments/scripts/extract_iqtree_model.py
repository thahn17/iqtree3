#!/usr/bin/env python3
"""
Turn a real iqtree3 run's fitted model parameters -- WHICHEVER model
ModelFinder (-m MFP/TEST/...) actually selected, not a fixed one chosen up
front -- into a literal IQ-TREE model spec string usable by
spr_topology_test's own "model <spec>" --hillclimb flag
(search_experiments/spr_topology_test.cpp), so a search there can start from the same
ML-estimated parameters the external iqtree3 run already fit, instead of
either spr_topology_test's own un-fit starting values or having no way to
supply them at all.

Reads the "Model of substitution: <name>" line (e.g. "LG+F+I+G4",
"WAG+FQ+R6", "Q.insect+FO+G4", "JTT" with no +F/+G/+R at all -- ANY
combination ModelFinder can return) and pulls in whichever of the
following are actually present, each carried over as a FIXED value via
the same "+X{...}" syntax IQ-TREE's own getNameParams methods use to
serialize a fitted model (model/rate{gamma,invar,free,...}.cpp):

  frequency  +F / +FO / +FC -> "pi(X) = ..." lines (State frequencies:)
             -> +F{f1,...,f20 (or 4 for DNA)}, same state order the
             report prints (alignment/alignment.cpp's symbols_protein,
             "ARNDCQEGHILKMFPSTWYV", for protein; A/C/G/T for DNA)
             +FQ (equal frequencies) -> kept bare, "+FQ": parameter-free,
             nothing to fit
             no +F suffix at all -> nothing added: the base matrix's own
             built-in frequencies, also nothing to fit
  invariable +I -> "Proportion of invariable sites: <p>" -> +I{p}
  sites
  gamma      +G<n> -> "Gamma shape alpha: <a>" -> +G<n>{a} (n taken from
                    the model name's own "+G<n>" digits, default 4 if
                    bare "+G")
  FreeRate   +R<n> -> "Site proportion and rates:  (p1,r1) (p2,r2) ..."
                    -> +R<n>{p1,r1,p2,r2,...} (n likewise from "+R<n>")

+I, +G and +R compose freely (I+G and I+R both occur; G+R does not, in
standard IQ-TREE); whichever combination the model name shows, that
combination -- in the SAME order the name shows it, always I before G/R,
matching RateGammaInvar::getNameParams/RateFreeInvar::getNameParams's own
"RateInvar::getNameParams() + RateGamma-or-Free::getNameParams()" -- is
what gets reproduced. A component present in the model name but whose
matching report line can't be found (shouldn't normally happen, but
malformed/truncated reports do) falls back to that component's bare,
unfit tag (e.g. "+G4" with no {alpha}) rather than failing outright --
ModelFactory will just re-estimate that one piece fresh, which is no
worse than not having this script's fit-carryover at all.

The base substitution matrix name is taken as everything before the
model name's first '+' (PhyloTree::getModelName()'s own
"<subst><freq><rate>" concatenation order -- e.g. "LG+F+I+G4" -- is
exactly what "Model of substitution:" prints, so hunting for a
"+F"-family suffix anchored to the END of the string, like an earlier
version of this script did, breaks the moment any rate suffix follows
it).

Usage:
    python3 extract_iqtree_model.py <run.iqtree>

Example, feeding the result straight into spr_topology_test:
    MODEL_SPEC="$(python3 search_experiments/scripts/extract_iqtree_model.py nni_run.iqtree)"
    build/spr_topology_test --hillclimb real.fa 10 10000 fast quiet notree \\
        starttree nni_run.treefile model "$MODEL_SPEC" gtr record
"""
import argparse
import re
import sys


def extract_freq_spec(full_model, text):
    """Return (freq_spec, base_model) -- freq_spec is "", "+FQ", or a fixed
    "+F{...}" built from this report's own "pi(X) = ..." lines (falling
    back to the bare matched tag, e.g. "+F", if none are found)."""
    freq_match = re.search(r"^([^+]+)(\+F([OCQ]?))?", full_model)
    base_model = freq_match.group(1)
    freq_tag = freq_match.group(2) or ""
    freq_kind = freq_match.group(3) or ""

    if not freq_tag:
        return "", base_model
    if freq_kind == "Q":
        return "+FQ", base_model  # equal frequencies: parameter-free, nothing to fit

    freqs = re.findall(r"^\s*pi\(([A-Za-z*]+)\)\s*=\s*([0-9.eE+-]+)", text, re.MULTILINE)
    if not freqs:
        return freq_tag, base_model  # e.g. bare "+F"/"+FO": couldn't find values to fit
    values = [v for _, v in freqs]
    return f"+F{{{','.join(values)}}}", base_model


def extract_rate_spec(full_model, text):
    """Return the +I/+G<n>/+R<n> suffix (any subset, I always first),
    each fixed to this report's own fitted value where found."""
    has_i = re.search(r"\+I(?:\+|$)", full_model) is not None
    gamma_match = re.search(r"\+G(\d*)(?:\+|$)", full_model)
    free_match = re.search(r"\+R(\d*)(?:\+|$)", full_model)

    spec = ""
    if has_i:
        pinvar_match = re.search(r"^Proportion of invariable sites:\s*([0-9.eE+-]+)", text, re.MULTILINE)
        spec += f"+I{{{pinvar_match.group(1)}}}" if pinvar_match else "+I"

    if gamma_match:
        ncategory = gamma_match.group(1) or "4"
        alpha_match = re.search(r"^Gamma shape alpha:\s*([0-9.eE+-]+)", text, re.MULTILINE)
        spec += f"+G{ncategory}{{{alpha_match.group(1)}}}" if alpha_match else f"+G{ncategory}"
    elif free_match:
        ncategory = free_match.group(1) or "4"
        rates_match = re.search(r"^Site proportion and rates:\s*(.+)$", text, re.MULTILINE)
        pairs = re.findall(r"\(([^,()]+),([^,()]+)\)", rates_match.group(1)) if rates_match else []
        if pairs:
            flat = ",".join(v for pair in pairs for v in pair)
            spec += f"+R{ncategory}{{{flat}}}"
        else:
            spec += f"+R{ncategory}"

    return spec


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
    full_model = model_match.group(1)

    freq_spec, base_model = extract_freq_spec(full_model, text)
    rate_spec = extract_rate_spec(full_model, text)

    print(f"{base_model}{freq_spec}{rate_spec}")


if __name__ == "__main__":
    main()
