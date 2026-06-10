"""Prune a large reference tree (e.g. GTDB bac120) down to a focal genome set.

For the bacterial co-occurrence / niche tests we need a species tree, but computing
one de novo over a bacterial reference set is heavy and unnecessary when the genomes
are GTDB representatives: the defensible, tractable choice is to PRUNE the published
GTDB bac120 reference tree to the focal genome_ids. Tip labels are matched exactly
(GTDB accession == genome_id); a low overlap is warned about loudly so a label-format
mismatch is visible rather than silently producing a tiny tree. Empty / zero-overlap
input yields a trivial '();' tree so the DAG stays valid.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


def load_keep(path):
    import pandas as pd
    if not path or not os.path.exists(path):
        return set()
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
        col = "genome_id" if "genome_id" in df.columns else df.columns[0]
        return {str(g).strip() for g in df[col].dropna() if str(g).strip()}
    except Exception:
        return set()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-tree", required=True)
    ap.add_argument("--keep", required=True, help="TSV with a genome_id column")
    ap.add_argument("--out-tree", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tree)), exist_ok=True)
    keep = load_keep(args.keep)

    if not args.in_tree or not os.path.exists(args.in_tree) or os.path.getsize(args.in_tree) == 0 \
            or not keep:
        with open(args.out_tree, "w") as fh:
            fh.write("();\n")
        utils.log("prune_tree: empty reference tree or empty keep set; wrote trivial tree.")
        return

    try:
        import dendropy
    except Exception:
        with open(args.out_tree, "w") as fh:
            fh.write("();\n")
        utils.log("prune_tree: dendropy unavailable; wrote trivial tree.")
        return

    try:
        tree = dendropy.Tree.get(path=args.in_tree, schema="newick",
                                 preserve_underscores=True)
    except Exception as exc:
        with open(args.out_tree, "w") as fh:
            fh.write("();\n")
        utils.log(f"prune_tree: could not parse reference tree ({exc}); trivial tree.")
        return

    tip_labels = {lf.taxon.label for lf in tree.leaf_node_iter() if lf.taxon}
    # match exactly; also try stripping a GTDB 'GB_'/'RS_' accession prefix both ways
    present = keep & tip_labels
    if len(present) < 0.5 * len(keep):
        # try matching ignoring common GTDB accession prefixes
        norm = lambda s: s.replace("GB_", "").replace("RS_", "")
        tip_norm = {norm(t): t for t in tip_labels}
        present = {tip_norm[norm(g)] for g in keep if norm(g) in tip_norm}
    if not present:
        with open(args.out_tree, "w") as fh:
            fh.write("();\n")
        utils.log(f"prune_tree: NONE of {len(keep)} genome_ids match the reference tree "
                  f"tip labels (format mismatch?); wrote trivial tree.")
        return
    if len(present) < 0.5 * len(keep):
        utils.log(f"prune_tree: WARNING only {len(present)}/{len(keep)} genome_ids matched "
                  f"the reference tree -- check that genome_ids are GTDB accessions.")

    try:
        tree.retain_taxa_with_labels(list(present))
        tree.write(path=args.out_tree, schema="newick", unquoted_underscores=True,
                   suppress_rooting=True)
        n = sum(1 for _ in tree.leaf_node_iter())
        utils.log(f"prune_tree: kept {n} tips of {len(keep)} requested -> {args.out_tree}")
    except Exception as exc:
        with open(args.out_tree, "w") as fh:
            fh.write("();\n")
        utils.log(f"prune_tree: retain failed ({exc}); wrote trivial tree.")


if __name__ == "__main__":
    main()
