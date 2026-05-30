"""Root a phylogeny at a defensible position before any comparative / ancestral use.

Why this exists
---------------
IQ-TREE emits an effectively unrooted tree (arbitrary basal placement). Several
downstream analyses are root-DEPENDENT and give wrong/undefined results on an
arbitrarily-rooted tree: stochastic mapping of gain/loss (phytools::make.simmap →
the independent-acquisition count and the methanogen-enrichment test), phylogenetic
logistic regression, Pagel's correlated-evolution test, the D statistic, and
GeneRax's UndatedDTL (which needs a rooted species tree). This step roots the tree
once, defensibly, so every consumer inherits the same justified root.

Rooting policy (in order of defensibility)
------------------------------------------
1. outgroup  : root so a user-designated outgroup is sister to the rest. The
   outgroup is given either as explicit genome_ids (`--outgroup`) or as a GTDB
   taxon (`--outgroup-taxon`, resolved to tips via the taxonomy table). This is
   the gold standard and is preferred. If the outgroup tips cannot be resolved in
   the tree, we DO NOT silently invent a root — we fall back to midpoint and say so.
2. midpoint  : a heuristic root at the midpoint of the longest tip-to-tip path.
   Reasonable when no outgroup is available, but it assumes roughly clock-like
   deepest divergence, which rate variation (e.g. fast-evolving DPANN/Nanoarchaeota)
   can violate. Use only as a fallback; an outgroup is preferred.

Robustness: trivial/empty trees are passed through unchanged; any failure falls
back to passing the input through (never breaks the DAG). The rooting actually
applied is written to --log and stderr so it is auditable / can go in the methods.
"""
from __future__ import annotations
import argparse
import os
import sys


def _log(msg: str, logpath: str | None = None) -> None:
    print(f"[root_tree] {msg}", file=sys.stderr)
    if logpath:
        try:
            with open(logpath, "a") as fh:
                fh.write(msg + "\n")
        except Exception:
            pass


def _read_text(path: str) -> str:
    try:
        with open(path) as fh:
            return fh.read().strip()
    except Exception:
        return ""


def resolve_outgroup(explicit: str, taxon: str, taxonomy: str, tip_labels: set) -> set:
    """Return the set of tree tips that constitute the outgroup."""
    og = set()
    if explicit:
        for g in explicit.split(","):
            g = g.strip()
            if g:
                og.add(g)
    if taxon and taxonomy and os.path.exists(taxonomy):
        tx = taxon.strip().lower()
        try:
            with open(taxonomy) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                gi = header.index("user_genome") if "user_genome" in header else 0
                ci = header.index("classification") if "classification" in header else 1
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) <= max(gi, ci):
                        continue
                    if tx in f[ci].lower():
                        og.add(f[gi].strip())
        except Exception as exc:
            _log(f"WARNING: could not parse taxonomy for outgroup taxon ({exc})")
    # keep only tips actually in the tree
    return {g for g in og if g in tip_labels}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--in-tree", required=True)
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--method", default="outgroup", help="outgroup | midpoint")
    ap.add_argument("--outgroup", default="", help="comma-separated genome_ids")
    ap.add_argument("--outgroup-taxon", default="", help="GTDB taxon (resolved via --taxonomy)")
    ap.add_argument("--taxonomy", default="", help="gtdbtk summary (user_genome, classification)")
    ap.add_argument("--log", default="", help="path to write the rooting actually applied")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_tree)), exist_ok=True)
    if args.log:
        os.makedirs(os.path.dirname(os.path.abspath(args.log)), exist_ok=True)
        open(args.log, "w").close()

    raw = _read_text(args.in_tree)
    # trivial / empty tree -> pass through unchanged
    if raw in ("", "();", ";") or raw.count(",") < 2:
        with open(args.out_tree, "w") as fh:
            fh.write((raw or "();") + "\n")
        _log(f"method=passthrough(trivial_tree) tips<3", args.log)
        return

    try:
        import dendropy
    except Exception as exc:
        with open(args.out_tree, "w") as fh:
            fh.write(raw + "\n")
        _log(f"method=passthrough(dendropy_unavailable:{exc})", args.log)
        return

    try:
        tree = dendropy.Tree.get(data=raw, schema="newick", preserve_underscores=True)
        # Mark as rooted so rerooting produces a real bifurcating root rather than
        # collapsing back to the basal trifurcation IQ-TREE emits.
        tree.is_rooted = True
        tip_labels = {lf.taxon.label for lf in tree.leaf_node_iter() if lf.taxon}
        n = len(tip_labels)
        method = (args.method or "outgroup").strip().lower()
        applied = None

        if method == "outgroup":
            og = resolve_outgroup(args.outgroup, args.outgroup_taxon, args.taxonomy, tip_labels)
            if 1 <= len(og) < n:
                labels = sorted(og)
                if len(labels) == 1:
                    node = next(lf for lf in tree.leaf_node_iter()
                                if lf.taxon and lf.taxon.label == labels[0])
                else:
                    node = tree.mrca(taxon_labels=labels)
                edge = getattr(node, "edge", None)
                if node is None or node is tree.seed_node or edge is None or edge.length is None:
                    tree.reroot_at_midpoint(update_bipartitions=True)
                    applied = "midpoint(outgroup_unrootable)"
                else:
                    # Root ON the outgroup edge (split it) so the root is bifurcating:
                    # one child is the outgroup clade, the other is everything else.
                    half = edge.length / 2.0
                    # tree.is_rooted is True (set above), so the new basal bifurcation
                    # is kept rather than collapsed. (kwarg-free for dendropy 4.x/5.x.)
                    tree.reroot_at_edge(edge, length1=half, length2=half,
                                        update_bipartitions=True, suppress_unifurcations=True)
                    under = {lf.taxon.label for lf in node.leaf_iter() if lf.taxon}
                    mono = "monophyletic" if under == set(labels) else "non-monophyletic@MRCA"
                    applied = f"outgroup({len(labels)} tips;{mono})"
            else:
                tree.reroot_at_midpoint(update_bipartitions=True)
                applied = "midpoint(outgroup_unresolved)" if (args.outgroup or args.outgroup_taxon) \
                    else "midpoint(no_outgroup_configured)"
        elif method == "midpoint":
            tree.reroot_at_midpoint(update_bipartitions=True)
            applied = "midpoint"
        else:
            tree.reroot_at_midpoint(update_bipartitions=True)
            applied = f"midpoint(unknown_method:{method})"

        tree.is_rooted = True
        newick = tree.as_string(schema="newick", suppress_rooting=False,
                                unquoted_underscores=True).strip()
        with open(args.out_tree, "w") as fh:
            fh.write(newick + "\n")
        _log(f"method={applied}; tips={n}", args.log)
    except Exception as exc:
        with open(args.out_tree, "w") as fh:
            fh.write(raw + "\n")
        _log(f"method=passthrough(rooting_failed:{exc})", args.log)


if __name__ == "__main__":
    main()
