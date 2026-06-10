"""Cross-domain Prometheus phylogeny: origin and bacteria<->archaea transfer direction.

Reads the combined cross-domain gene tree (archaeal + bacterial validated Prometheus
proteins + the Pseudomonas Prometheus reference + DNA2/NurA-HerA paralog outgroup) and
summarises, from TOPOLOGY ALONE (no DTL model, deliberately -- a single rooted species
tree spanning both domains is contested, so we do not over-model):

  * tips per domain (Archaea / Bacteria, with Prometheus_ref counted as Bacteria);
  * minimum number of cross-domain transitions (Fitch small-parsimony on the
    Bacteria/Archaea character over the Prometheus ingroup, rooted on the paralog
    outgroup);
  * whether the archaeal sequences are monophyletic, and the number of maximal
    archaeal-only clades = the count of independent archaeal acquisitions;
  * a polarity call: if the ingroup root state is Bacterial and the archaeal sequences
    form derived/nested clades, the direction is bacteria->archaea (the expected result
    given Prometheus was characterised in bacteria); otherwise 'ambiguous'.

This answers the question the archaea-only design cannot: where the archaeal copies
came from. It is a robustness-oriented, assumption-light complement to within-domain
reconciliation, NOT a significance test. Empty/trivial input -> header-only + note.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SUMMARY_HEADER = "metric\tvalue\tnote"
TIPS_HEADER = "tip\tdomain"


def load_id_set(fasta_path):
    return set(utils.fasta_ids(fasta_path)) if fasta_path and os.path.exists(fasta_path) else set()


def load_ref_labels(path):
    labels = {}
    if not path or not os.path.exists(path):
        return labels
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            if parts[1].strip() in ("Prometheus_ref", "Paralog_ref"):
                labels[parts[0].strip()] = parts[1].strip()
    return labels


def fitch(tree, state_of):
    """2-state Fitch small-parsimony over a rooted dendropy tree.

    state_of: dict tip_label -> 'A' or 'B' (ingroup tips only; tree already pruned to
    them). Returns (n_changes, root_state_set). Postorder union/intersection.
    """
    sets = {}
    n_changes = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            lbl = node.taxon.label if node.taxon else None
            sets[node] = {state_of.get(lbl)} if lbl in state_of else set()
        else:
            child_sets = [sets[c] for c in node.child_nodes() if sets.get(c)]
            if not child_sets:
                sets[node] = set()
                continue
            inter = set.intersection(*child_sets)
            if inter:
                sets[node] = inter
            else:
                sets[node] = set.union(*child_sets)
                n_changes += 1
    root_set = sets.get(tree.seed_node, set())
    return n_changes, root_set


def count_archaeal_clades(tree, domain_of):
    """Number of MAXIMAL archaeal-only clades (= independent archaeal acquisitions).

    A node is archaeal-only if every descendant tip is Archaea. A clade is maximal when
    the node is archaeal-only but its parent is not. Counting those gives the number of
    distinct archaeal groups on the tree.
    """
    all_arch = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            lbl = node.taxon.label if node.taxon else None
            all_arch[node] = (domain_of.get(lbl) == "Archaea")
        else:
            all_arch[node] = all(all_arch[c] for c in node.child_nodes()) and bool(node.child_nodes())
    n = 0
    for node in tree.preorder_node_iter():
        if not all_arch.get(node):
            continue
        parent = node.parent_node
        if parent is None or not all_arch.get(parent):
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--arch-proteins", required=True, help="archaeal validated proteins FASTA")
    ap.add_argument("--bact-proteins", required=True, help="bacterial validated proteins FASTA")
    ap.add_argument("--ref-labels", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-tips", required=True)
    args = ap.parse_args()

    for p in (args.out, args.out_tips):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    def bail(note):
        utils.write_empty(args.out, SUMMARY_HEADER)
        with open(args.out, "a") as fh:
            fh.write(f"status\tincomplete\t{note}\n")
        utils.write_empty(args.out_tips, TIPS_HEADER)
        utils.log(f"crossdomain_analysis: {note}")

    txt = ""
    if os.path.exists(args.tree) and os.path.getsize(args.tree) > 0:
        with open(args.tree) as fh:
            txt = fh.read().strip()
    if not txt or txt in ("();", ";") or txt.count(",") < 2:
        bail("trivial/empty cross-domain tree")
        return
    try:
        import dendropy
    except Exception:
        bail("dendropy unavailable")
        return

    arch_ids = load_id_set(args.arch_proteins)
    bact_ids = load_id_set(args.bact_proteins)
    ref_labels = load_ref_labels(args.ref_labels)

    try:
        tree = dendropy.Tree.get(data=txt, schema="newick", preserve_underscores=True)
    except Exception as exc:
        bail(f"could not parse tree ({exc})")
        return

    # domain per tip: Prometheus_ref counts as Bacteria (it IS bacterial Prometheus);
    # Paralog_ref is the outgroup (excluded from the ingroup character).
    domain_of = {}
    for lf in tree.leaf_node_iter():
        if not lf.taxon:
            continue
        lbl = lf.taxon.label
        if lbl in ref_labels:
            domain_of[lbl] = "Paralog" if ref_labels[lbl] == "Paralog_ref" else "Bacteria"
        elif lbl in arch_ids:
            domain_of[lbl] = "Archaea"
        elif lbl in bact_ids:
            domain_of[lbl] = "Bacteria"
        else:
            domain_of[lbl] = "Unknown"

    with open(args.out_tips, "w") as fh:
        fh.write(TIPS_HEADER + "\n")
        for lbl, dom in sorted(domain_of.items()):
            fh.write(f"{lbl}\t{dom}\n")

    n_arch = sum(1 for d in domain_of.values() if d == "Archaea")
    n_bact = sum(1 for d in domain_of.values() if d == "Bacteria")
    paralogs = [l for l, d in domain_of.items() if d == "Paralog"]

    if n_arch == 0 or n_bact == 0:
        bail(f"need both domains in the tree (archaea={n_arch}, bacteria={n_bact})")
        return

    # ---- root on the paralog outgroup (else midpoint) -------------------- #
    tree.is_rooted = True
    try:
        if len(paralogs) >= 1:
            if len(paralogs) == 1:
                mrca = next(lf for lf in tree.leaf_node_iter()
                            if lf.taxon and lf.taxon.label == paralogs[0])
            else:
                mrca = tree.mrca(taxon_labels=paralogs)
            edge = getattr(mrca, "edge", None)
            if mrca is not None and mrca is not tree.seed_node and edge is not None and edge.length:
                tree.reroot_at_edge(edge, length1=edge.length / 2.0, length2=edge.length / 2.0,
                                    update_bipartitions=True, suppress_unifurcations=True)
                rooting = f"outgroup({len(paralogs)} paralog tips)"
            else:
                tree.reroot_at_midpoint(update_bipartitions=True); rooting = "midpoint(outgroup_unrootable)"
        else:
            tree.reroot_at_midpoint(update_bipartitions=True); rooting = "midpoint(no_paralog)"
    except Exception as exc:
        rooting = f"unrooted({exc})"

    # ---- ingroup = Archaea + Bacteria tips (drop paralog + unknown) ------ #
    ingroup_labels = [l for l, d in domain_of.items() if d in ("Archaea", "Bacteria")]
    ing = tree.extract_tree_with_taxa_labels(ingroup_labels)
    ing.is_rooted = True
    state_of = {l: ("A" if domain_of[l] == "Archaea" else "B") for l in ingroup_labels}

    n_changes, root_set = fitch(ing, state_of)
    n_arch_clades = count_archaeal_clades(ing, domain_of)
    archaeal_mono = (n_arch_clades == 1)

    if root_set == {"B"}:
        direction = "bacteria->archaea"
        dir_note = "ingroup root is bacterial; archaeal sequences are derived/nested"
    elif root_set == {"A"}:
        direction = "archaea->bacteria_or_archaeal_root"
        dir_note = "ingroup root resolves archaeal (unexpected; check outgroup/rooting)"
    else:
        direction = "ambiguous"
        dir_note = f"ingroup root state unresolved ({sorted(root_set)})"

    rows = [
        ("n_archaeal_tips", n_arch, "validated archaeal Prometheus on the cross-domain tree"),
        ("n_bacterial_tips", n_bact, "validated bacterial Prometheus + Pseudomonas reference"),
        ("n_crossdomain_transitions", n_changes,
         f"Fitch minimum Bacteria<->Archaea changes on the ingroup; rooting={rooting}"),
        ("n_archaeal_clades", n_arch_clades,
         "maximal archaeal-only clades = independent archaeal acquisitions"),
        ("archaeal_monophyletic", archaeal_mono,
         "True = a single archaeal clade (one acquisition or vertical within archaea)"),
        ("ingroup_root_state", "/".join(sorted(root_set)) or "NA", "Fitch state at the ingroup root"),
        ("inferred_direction", direction, dir_note),
    ]
    with open(args.out, "w") as fh:
        fh.write(SUMMARY_HEADER + "\n")
        for m, v, note in rows:
            fh.write(f"{m}\t{v}\t{note}\n")

    utils.log(f"crossdomain: archaea={n_arch} bacteria={n_bact} transitions={n_changes} "
              f"archaeal_clades={n_arch_clades} direction={direction}")


if __name__ == "__main__":
    main()
