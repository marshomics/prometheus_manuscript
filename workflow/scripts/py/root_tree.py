"""Root (and finalise) the species tree before any comparative / ancestral use.

Why this exists
---------------
IQ-TREE emits an effectively unrooted tree (a basal trifurcation at an arbitrary
place). Several downstream analyses are root-DEPENDENT and give wrong or undefined
results on an arbitrarily rooted tree: stochastic mapping of gains/losses
(phytools::make.simmap -> the independent-acquisition count and the
methanogen-enrichment test), phylogenetic logistic regression, Pagel's test, the D
statistic, and GeneRax's UndatedDTL (which needs a rooted species tree). This step
roots the tree once, defensibly, so every consumer inherits the same justified root.

What the canonical output guarantees
------------------------------------
1. Tips are analysed genome_ids. IQ-TREE 2 rewrites every character outside
   [A-Za-z0-9_-.|/] in a sequence name to '_' (verified on IQ-TREE 2.3.6). With
   --keep-tips the original genome_id is restored for every tip whose rewritten form
   maps back to exactly one analysed genome. A tip that is neither an analysed genome
   nor a GTDB reference genome (--ref-taxonomy) is an error.
2. The root is BIFURCATING, so ape::is.rooted() is TRUE. The R consumers midpoint-root
   any tree that ape reports as unrooted; a basal trifurcation would silently replace
   the configured root with a midpoint root.
3. Branch-support labels (IQ-TREE "SH-aLRT/UFBoot") stay attached to their
   bipartitions. DendroPy 4.6.1 leaves internal-node labels on the nodes when it
   reroots, which moves support values onto the wrong splits (reproduced on an IQ-TREE
   2.3.6 tree in tests/scaffolds; the general problem: Czech, Huerta-Cepas & Stamatakis
   2017, Mol Biol Evol 34:1535-1542, doi:10.1093/molbev/msx055, verified via CrossRef).
   Every label is therefore keyed by its
   bipartition before rooting and re-assigned afterwards. The split on which the root
   sits is represented by both children of the root; both carry its support value.
4. The Newick is written without a '[&R]' token and without quoting labels such as
   '95.3/100' (ape keeps quote characters as part of a label). No root branch length.

Rooting policy (in order of defensibility)
------------------------------------------
1. outgroup : the outgroup is given as explicit tip ids (--outgroup) and/or as GTDB
   taxa (--outgroup-taxon, comma-separated). A taxon is matched EXACTLY against the
   rank tokens of the GTDB classification ('p__Altiarchaeota' does not match
   'p__Altiarchaeota_A', which GTDB treats as a different lineage); a bare name without
   a rank prefix matches that name at any rank. Taxa are resolved through --taxonomy
   (the analysed genomes) and --ref-taxonomy (GTDB reference genomes, tree.scope=domain).
   The outgroup must be present and MONOPHYLETIC in the unrooted tree, i.e. one edge
   separates exactly the outgroup tips from all other tips; the root is placed at the
   midpoint of that edge. If a configured outgroup cannot be used (absent, equal to all
   tips, or not monophyletic) the script FAILS (--on-outgroup-failure fail, the
   default) and lists the ingroup tips that intrude into the outgroup clade. Only when
   explicitly requested (--on-outgroup-failure midpoint) does it fall back to midpoint
   rooting, and the fallback is recorded.
   If no outgroup is configured at all, midpoint rooting is applied and recorded as
   'midpoint(no_outgroup_configured)'.
2. midpoint : root at the midpoint of the longest tip-to-tip path. A heuristic: it
   assumes a roughly clock-like deepest divergence, which rate variation (e.g.
   fast-evolving DPANN) can violate. Use it only when no outgroup is available.

Pruning and membership
----------------------
--keep-tips lists the analysed genomes. After rooting, tips that are not analysed
genomes (GTDB reference genomes when tree.scope=domain) are pruned, so every consumer
sees only analysed genomes; --out-full-tree keeps the unpruned rooted tree.
--out-membership writes one row per analysed genome (in_msa, in_tree, status) so a
genome missing from the tree is never silent (GTDB-Tk --min_perc_aa filter, bac120
assignment, no markers, identify failure: --gtdbtk-dropped; or an MSA too small for
a tree).

Failure policy
--------------
An unreadable tree, missing dendropy, an unknown method, an unknown tip label or an
unusable configured outgroup are errors (exit 1). Only a genuinely trivial tree
('();' or fewer than 3 tips, because the MSA held fewer than 4 genomes) is passed
through unchanged, and that is recorded in --log.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

IQTREE_UNSAFE = re.compile(r"[^A-Za-z0-9_\-.|/]")
# Characters that force quoting of an unquoted Newick label (Newick/PHYLIP rules).
NEWICK_UNSAFE = re.compile(r"[\s()\[\]':;,]")
RANK_TOKEN = re.compile(r"^[dpcofgs]__", re.IGNORECASE)
VALID_METHODS = ("outgroup", "midpoint")


class RootingError(Exception):
    """A configured rooting request that cannot be honoured."""


def _say(msg: str) -> None:
    print(f"[root_tree] {msg}", file=sys.stderr)


def iqtree_name(name: str) -> str:
    """The sequence name IQ-TREE 2 uses for `name` (unsafe characters -> '_')."""
    return IQTREE_UNSAFE.sub("_", name)


# --------------------------------------------------------------------------- #
# small readers
# --------------------------------------------------------------------------- #
def read_id_list(path: str) -> list:
    """One id per line (first whitespace token); blank / '#' lines and a leading
    'genome_id' header line are skipped. Order preserved, duplicates removed."""
    out, seen = [], set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            tok = line.strip().split("\t")[0].strip()
            if not tok or tok.startswith("#"):
                continue
            if i == 0 and tok == "genome_id":
                continue
            if tok not in seen:
                seen.add(tok)
                out.append(tok)
    return out


def read_fasta_ids(path: str) -> list:
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                tok = line[1:].strip().split()
                if tok:
                    ids.append(tok[0])
    return ids


def read_taxonomy(path: str) -> dict:
    """id -> GTDB classification. Accepts a headed table (user_genome|genome_id +
    classification) or a headerless two-column file (GTDB style)."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        return out
    first = lines[0].split("\t")
    id_col, cls_col, start = 0, 1, 0
    low = [c.strip().lower() for c in first]
    if "classification" in low or "user_genome" in low or "genome_id" in low:
        start = 1
        for name in ("user_genome", "genome_id"):
            if name in low:
                id_col = low.index(name)
                break
        if "classification" in low:
            cls_col = low.index("classification")
        else:
            raise ValueError(f"taxonomy {path}: header has no 'classification' column")
    for ln in lines[start:]:
        f = ln.split("\t")
        if len(f) <= max(id_col, cls_col):
            continue
        gid = f[id_col].strip()
        if gid:
            out[gid] = f[cls_col].strip()
    return out


def read_dropped(path: str) -> dict:
    """gtdbtk_dropped.tsv (genome_id<TAB>reason, no header) -> {id: 'r1; r2'}."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t", 1)
            if not f[0].strip() or f[0].startswith("#") or f[0] == "genome_id":
                continue
            reason = f[1].strip() if len(f) > 1 else "unspecified"
            gid = f[0].strip()
            out[gid] = f"{out[gid]}; {reason}" if gid in out else reason
    return out


def resolve_taxa(spec: str, taxonomies: list) -> tuple:
    """Ids whose classification contains one of the taxa in `spec` (comma list).

    Exact token matching: 'p__X' must equal a rank token of the classification
    (case-insensitive); a bare 'X' must equal the name part of any rank token.
    Returns (ids, {taxon: n_matched}).
    """
    ids, counts = set(), {}
    taxa = [t.strip() for t in (spec or "").split(",") if t.strip()]
    for taxon in taxa:
        want = taxon.casefold()
        ranked = bool(RANK_TOKEN.match(taxon))
        n = 0
        for table in taxonomies:
            for gid, cls in table.items():
                for tok in cls.split(";"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    if ranked:
                        hit = tok.casefold() == want
                    else:
                        hit = (tok.split("__", 1)[-1].casefold() == want)
                    if hit:
                        ids.add(gid)
                        n += 1
                        break
        counts[taxon] = n
    return ids, counts


# --------------------------------------------------------------------------- #
# tree helpers (dendropy for parsing / rerooting / pruning; own split bookkeeping
# and writer)
# --------------------------------------------------------------------------- #
def leaf_masks(tree, index: dict) -> dict:
    """id(node) -> bitmask of the leaves below it (postorder, iterative)."""
    masks = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            masks[id(nd)] = 1 << index[nd.taxon.label]
        else:
            m = 0
            for ch in nd.child_node_iter():
                m |= masks[id(ch)]
            masks[id(nd)] = m
    return masks


def canonical_split(mask: int, full: int) -> int:
    """Orientation-free key of a bipartition: the side WITHOUT leaf 0."""
    return (full ^ mask) if (mask & 1) else mask


def record_support(tree, index: dict, full: int) -> dict:
    """{canonical split: label} for every labelled internal (non-seed) node."""
    masks = leaf_masks(tree, index)
    out = {}
    for nd in tree.preorder_internal_node_iter():
        if nd is tree.seed_node:
            continue
        lab = nd.label
        if lab is None or str(lab).strip() == "":
            continue
        out[canonical_split(masks[id(nd)], full)] = str(lab)
    return out


def reassign_support(tree, index: dict, full: int, support: dict) -> int:
    """Re-attach support labels by bipartition after rooting; return #labelled nodes."""
    masks = leaf_masks(tree, index)
    n = 0
    for nd in tree.preorder_internal_node_iter():
        if nd is tree.seed_node:
            nd.label = None
            continue
        lab = support.get(canonical_split(masks[id(nd)], full))
        nd.label = lab
        n += lab is not None
    return n


def _fmt_len(x) -> str:
    return "%.10f" % float(x)


def _fmt_label(s) -> str:
    s = "" if s is None else str(s)
    if s and NEWICK_UNSAFE.search(s):
        return "'" + s.replace("'", "''") + "'"
    return s


def to_newick(tree) -> str:
    """Iterative Newick writer: no [&R], no root branch length, labels unquoted
    unless Newick requires quoting, branch lengths with 10 decimals (IQ-TREE style)."""
    out = []
    stack = [(tree.seed_node, False)]
    while stack:
        nd, done = stack.pop()
        if nd is None:              # separator between siblings
            out.append(",")
            continue
        kids = nd.child_nodes()
        if kids and not done:
            stack.append((nd, True))
            out.append("(")
            for i, ch in enumerate(reversed(kids)):
                stack.append((ch, False))
                if i < len(kids) - 1:
                    stack.append((None, None))
            continue
        if kids:
            out.append(")")
            if nd is not tree.seed_node and nd.label:
                out.append(_fmt_label(nd.label))
        else:
            out.append(_fmt_label(nd.taxon.label if nd.taxon else ""))
        if nd is not tree.seed_node and nd.edge.length is not None:
            out.append(":" + _fmt_len(nd.edge.length))
    return "".join(out) + ";"


def tip_labels(tree) -> list:
    return [lf.taxon.label for lf in tree.leaf_node_iter() if lf.taxon is not None]


def root_degree(tree) -> int:
    return len(tree.seed_node.child_nodes())


def _neighbours(tree):
    """Undirected adjacency {id(node): [(node, length), ...]} (length None -> error)."""
    adj = {}
    for nd in tree.preorder_node_iter():
        adj.setdefault(id(nd), [])
        if nd.parent_node is not None:
            if nd.edge.length is None:
                raise ValueError("midpoint rooting needs branch lengths on every edge")
            L = float(nd.edge.length)
            adj[id(nd)].append((nd.parent_node, L))
            adj.setdefault(id(nd.parent_node), []).append((nd, L))
    return adj


def _farthest(start, adj):
    """(farthest leaf, its distance, predecessor map) from `start`; ties -> smallest label."""
    dist, pred = {id(start): 0.0}, {id(start): None}
    nodes = {id(start): start}
    stack = [start]
    while stack:
        nd = stack.pop()
        for nb, L in adj[id(nd)]:
            if id(nb) not in dist:
                dist[id(nb)] = dist[id(nd)] + L
                pred[id(nb)] = nd
                nodes[id(nb)] = nb
                stack.append(nb)
    leaves = [n for n in nodes.values() if n.is_leaf() and n is not start]
    dmax = max(dist[id(n)] for n in leaves)
    tied = [n for n in leaves if dist[id(n)] >= dmax - 1e-12]
    best = min(tied, key=lambda n: n.taxon.label)   # deterministic tie-break
    return best, dist[id(best)], pred


def reroot_midpoint(tree) -> float:
    """Root at the midpoint of the longest tip-to-tip path (bifurcating, always).

    DendroPy 4.6.1's reroot_at_midpoint can turn a tip into the root node or leave
    a multifurcating root when the midpoint falls on an existing node, so the
    midpoint is located here and the root is placed with reroot_at_edge. When the
    midpoint coincides with a node, the root sits on the adjacent path edge at
    distance 0 from that node (a zero-length branch keeps the root bifurcating).
    Returns the tree diameter.
    """
    adj = _neighbours(tree)
    first = sorted(tree.leaf_node_iter(), key=lambda n: n.taxon.label)[0]
    a, _, _ = _farthest(first, adj)
    b, diameter, pred = _farthest(a, adj)
    if diameter <= 0:
        raise ValueError("cannot midpoint-root a tree whose longest path has length 0")
    path = [b]
    while pred[id(path[-1])] is not None:
        path.append(pred[id(path[-1])])
    path.reverse()                                  # a ... b
    half, acc = diameter / 2.0, 0.0
    for u, v in zip(path[:-1], path[1:]):
        L = next(l for nb, l in adj[id(u)] if nb is v)
        if acc + L > half or v is b:
            s = min(max(half - acc, 0.0), L)        # distance from u to the root
            if v.parent_node is u:                  # edge belongs to v (head v, tail u)
                edge, tail_len, head_len = v.edge, s, L - s
            else:                                   # edge belongs to u (head u, tail v)
                edge, tail_len, head_len = u.edge, L - s, s
            tree.reroot_at_edge(edge, length1=tail_len, length2=head_len,
                                update_bipartitions=False, suppress_unifurcations=True)
            return diameter
        acc += L
    raise RuntimeError("midpoint not found on the longest path")


def find_split_node(tree, index: dict, full: int, og_mask: int):
    """Non-seed node whose edge separates exactly the outgroup from the rest."""
    masks = leaf_masks(tree, index)
    comp = full ^ og_mask
    for nd in tree.preorder_node_iter():
        if nd is tree.seed_node:
            continue
        m = masks[id(nd)]
        if m == og_mask or m == comp:
            return nd
    return None


def outgroup_intruders(tree, index: dict, outgroup: set) -> list:
    """Ingroup tips inside the smallest clade holding the outgroup when the tree is
    viewed from an ingroup tip (diagnostic for a non-monophyletic outgroup)."""
    import copy
    t = copy.deepcopy(tree)
    ingroup = sorted(set(tip_labels(t)) - outgroup)
    if not ingroup:
        return []
    anchor = next(lf for lf in t.leaf_node_iter() if lf.taxon.label == ingroup[0])
    t.is_rooted = True
    t.reroot_at_edge(anchor.edge, length1=None, length2=None,
                     update_bipartitions=False, suppress_unifurcations=True)
    node = t.mrca(taxon_labels=sorted(outgroup))
    below = {lf.taxon.label for lf in node.leaf_iter()}
    return sorted(below - outgroup)


# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-tree", required=True)
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--method", default="outgroup", help="outgroup | midpoint")
    ap.add_argument("--outgroup", default="", help="comma-separated tip ids (genome_ids)")
    ap.add_argument("--outgroup-taxon", default="",
                    help="GTDB taxon/taxa (comma list), exact rank-token match")
    ap.add_argument("--taxonomy", default="",
                    help="analysed-genome taxonomy (user_genome, classification)")
    ap.add_argument("--ref-taxonomy", default="",
                    help="GTDB reference taxonomy (id<TAB>classification), scope=domain")
    ap.add_argument("--on-outgroup-failure", default="fail", choices=("fail", "midpoint"),
                    help="configured outgroup unusable: fail (default) or midpoint fallback")
    ap.add_argument("--keep-tips", default="",
                    help="analysed genome ids (one per line): restore IQ-TREE-renamed "
                         "labels, validate tips, prune non-analysed tips after rooting")
    ap.add_argument("--msa", default="", help="FASTA the tree was inferred from (membership)")
    ap.add_argument("--gtdbtk-dropped", default="",
                    help="genome_id<TAB>reason for genomes GTDB-Tk left out of the MSA")
    ap.add_argument("--out-full-tree", default="",
                    help="rooted tree BEFORE pruning to --keep-tips")
    ap.add_argument("--out-membership", default="",
                    help="per analysed genome: in_msa, in_tree, status")
    ap.add_argument("--log", default="", help="path to write the rooting actually applied")
    return ap


def write_text(path: str, text: str) -> None:
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def write_membership(path, expected, msa_ids, final_tips, dropped, trivial) -> dict:
    """Write the membership table; return {status: count}."""
    counts = {}
    rows = ["genome_id\tin_msa\tin_tree\tstatus"]
    msa_set, tip_set = set(msa_ids), set(final_tips)
    for g in expected:
        in_msa, in_tree = g in msa_set, g in tip_set
        if in_tree:
            status = "in_tree"
        elif in_msa and trivial:
            status = "msa_too_small_for_tree(<4_genomes)"
        elif in_msa:
            raise RuntimeError(f"genome {g} is in the MSA but not in the inferred tree")
        elif g in dropped:
            status = "dropped_by_gtdbtk: " + dropped[g]
        else:
            status = "not_in_gtdbtk_msa: no reason reported"
        key = status.split(":")[0]
        counts[key] = counts.get(key, 0) + 1
        rows.append(f"{g}\t{in_msa}\t{in_tree}\t{status}")
    write_text(path, "\n".join(rows) + "\n")
    return counts


def main(argv=None) -> int:
    ap = build_parser()
    args = ap.parse_args(argv)
    if args.out_membership and not (args.keep_tips and args.msa):
        ap.error("--out-membership needs --keep-tips and --msa")
    method = (args.method or "").strip().lower()
    if method not in VALID_METHODS:
        _say(f"ERROR: unknown --method {args.method!r} (expected one of {VALID_METHODS})")
        return 1

    for p in (args.out_tree, args.out_full_tree, args.out_membership, args.log):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)) or ".", exist_ok=True)

    with open(args.in_tree) as fh:          # unreadable input -> exception -> exit 1
        raw = fh.read().strip()
    expected = read_id_list(args.keep_tips) if args.keep_tips else []
    msa_ids = read_fasta_ids(args.msa) if args.msa else []
    dropped = read_dropped(args.gtdbtk_dropped) if args.gtdbtk_dropped else {}
    info = {"requested_method": method, "outgroup_ids": args.outgroup.strip(),
            "outgroup_taxon": args.outgroup_taxon.strip(),
            "on_outgroup_failure": args.on_outgroup_failure}

    import dendropy  # a missing package is an error, not a pass-through

    trivial = raw in ("", ";", "();")
    tree = None
    if not trivial:
        tree = dendropy.Tree.get(data=raw, schema="newick", preserve_underscores=True)
        if len(tip_labels(tree)) < 3:
            trivial = True
    if trivial:
        # nonfatal: genuinely trivial input (MSA had < 4 genomes); recorded below
        write_text(args.out_tree, (raw or "();") + "\n")
        write_text(args.out_full_tree, (raw or "();") + "\n")
        counts = {}
        if args.out_membership:
            counts = write_membership(args.out_membership, expected, msa_ids, [], dropped,
                                      trivial=True)
        n = 0 if tree is None else len(tip_labels(tree))
        lines = [f"method=none(trivial_tree;tips={n})", f"n_tips_input={n}",
                 f"n_tips_output={n}"] + [f"{k}={v}" for k, v in info.items()]
        lines += [f"membership_{k}={v}" for k, v in sorted(counts.items())]
        write_text(args.log, "\n".join(lines) + "\n")
        _say(lines[0])
        return 0

    # ---- tip labels: restore IQ-TREE renames, validate ------------------------
    ref_tax = read_taxonomy(args.ref_taxonomy) if args.ref_taxonomy else {}
    user_tax = read_taxonomy(args.taxonomy) if args.taxonomy else {}
    relabelled = 0
    if expected:
        expected_set = set(expected)
        back = {}
        for g in list(expected) + sorted(ref_tax):
            back.setdefault(iqtree_name(g), []).append(g)
        for lf in tree.leaf_node_iter():
            lab = lf.taxon.label
            if lab in expected_set or lab in ref_tax:
                continue
            orig = back.get(lab, [])
            if len(orig) == 1:
                if NEWICK_UNSAFE.search(orig[0]):
                    raise ValueError(f"genome_id {orig[0]!r} cannot be written as an unquoted "
                                     "Newick label; rename it in the input metadata")
                lf.taxon.label = orig[0]
                relabelled += 1
            elif len(orig) > 1:
                raise ValueError(f"tip {lab!r} is the IQ-TREE rename of several ids "
                                 f"{orig}; make the ids distinct")
            else:
                raise ValueError(f"tip {lab!r} is neither an analysed genome (--keep-tips) "
                                 "nor a GTDB reference genome (--ref-taxonomy)")
    tips = tip_labels(tree)
    if len(set(tips)) != len(tips):
        raise ValueError("duplicate tip labels in the input tree")
    n_in = len(tips)

    # ---- bipartition bookkeeping ----------------------------------------------
    if root_degree(tree) == 2:
        # an already-rooted input: unroot so that every edge is one bipartition
        tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    index = {lab: i for i, lab in enumerate(sorted(tips))}
    full = (1 << len(index)) - 1
    support = record_support(tree, index, full)
    tree.is_rooted = True  # keep the new basal bifurcation when rerooting

    # ---- rooting ----------------------------------------------------------------
    applied, fallback = None, ""
    if method == "outgroup":
        configured = bool(args.outgroup.strip() or args.outgroup_taxon.strip())
        if not configured:
            reroot_midpoint(tree)
            applied = "midpoint(no_outgroup_configured)"
            _say("WARNING: rooting_method=outgroup but no tree.outgroup / tree.outgroup_taxon "
                 "configured: midpoint-rooted. Set an outgroup for a defensible root.")
        else:
            requested = {g.strip() for g in args.outgroup.split(",") if g.strip()}
            tax_ids, tax_counts = resolve_taxa(args.outgroup_taxon, [user_tax, ref_tax])
            info["outgroup_taxon_matches"] = ",".join(f"{k}:{v}" for k, v in tax_counts.items())
            wanted = requested | tax_ids
            tipset = set(tips)
            og = wanted & tipset
            missing = sorted(requested - tipset)
            if missing:
                info["outgroup_ids_not_in_tree"] = ",".join(missing[:50])
            try:
                if not og:
                    raise RootingError("outgroup_not_in_tree: no configured outgroup tip "
                                       f"is in the tree (requested ids: {sorted(requested)[:10]}; "
                                       f"taxon matches: {tax_counts})")
                if og == tipset:
                    raise RootingError("outgroup_is_all_tips: the outgroup contains every tip")
                og_mask = 0
                for g in og:
                    og_mask |= 1 << index[g]
                node = find_split_node(tree, index, full, og_mask)
                if node is None:
                    intr = outgroup_intruders(tree, index, og)
                    info["outgroup_intruders"] = ",".join(intr[:50])
                    raise RootingError(
                        f"outgroup_not_monophyletic: {len(og)} outgroup tips do not form a "
                        f"clade; {len(intr)} ingroup tip(s) nest inside the outgroup clade "
                        f"(first: {intr[:10]})")
                e = node.edge
                half = None if e.length is None else e.length / 2.0
                tree.reroot_at_edge(e, length1=half, length2=half,
                                    update_bipartitions=False, suppress_unifurcations=True)
                applied = f"outgroup({len(og)} tips;monophyletic)"
                info["outgroup_tips"] = ",".join(sorted(og)[:200])
                info["outgroup_monophyletic"] = "True"
            except RootingError as exc:
                fallback = str(exc)
                info["outgroup_monophyletic"] = "False"
                if args.on_outgroup_failure != "midpoint":
                    msg = (f"ERROR: configured outgroup cannot root the tree -> {fallback}. "
                           "Fix tree.outgroup / tree.outgroup_taxon, or set "
                           "tree.outgroup_on_failure: midpoint to accept a labelled midpoint root.")
                    lines = ["method=failed", f"failure={fallback}"] + \
                            [f"{k}={v}" for k, v in info.items()]
                    write_text(args.log, "\n".join(lines) + "\n")
                    _say(msg)
                    return 1
                reroot_midpoint(tree)
                applied = "midpoint(fallback:" + fallback.split(":")[0] + ")"
                _say(f"WARNING: {fallback}; midpoint-rooted as requested "
                     "(--on-outgroup-failure midpoint).")
    else:
        reroot_midpoint(tree)
        applied = "midpoint"

    if root_degree(tree) != 2:
        raise RuntimeError(f"rooting produced a root of degree {root_degree(tree)} (expected 2)")
    tree.seed_node.edge.length = None
    n_support = reassign_support(tree, index, full, support)
    full_newick = to_newick(tree)
    write_text(args.out_full_tree, full_newick + "\n")

    # ---- prune to the analysed genomes ------------------------------------------
    n_pruned = 0
    if expected:
        keep = set(expected)
        drop = [t for t in tips if t not in keep]
        n_pruned = len(drop)
        if drop:
            if len(tips) - len(drop) < 2:
                raise ValueError("fewer than 2 analysed genomes remain after pruning references")
            tree.prune_taxa_with_labels(drop, update_bipartitions=False,
                                        suppress_unifurcations=True)
            if root_degree(tree) != 2:
                raise RuntimeError("pruning produced a non-bifurcating root")
            tree.seed_node.edge.length = None
            tree.seed_node.label = None
    final_tips = tip_labels(tree)
    write_text(args.out_tree, to_newick(tree) + "\n")

    counts = {}
    if args.out_membership:
        counts = write_membership(args.out_membership, expected, msa_ids, final_tips,
                                  dropped, trivial=False)
    kids = tree.seed_node.child_nodes()
    sizes = sorted(len(list(k.leaf_iter())) for k in kids)
    lines = [f"method={applied}; tips={len(final_tips)}",
             f"method_applied={applied}", f"fallback_reason={fallback}",
             f"n_tips_input={n_in}", f"n_reference_tips_pruned={n_pruned}",
             f"n_tips_output={len(final_tips)}", f"relabelled_iqtree_tips={relabelled}",
             f"support_labels_reassigned={n_support}",
             f"root_child_clade_sizes={sizes[0]},{sizes[1]}"]
    lines += [f"{k}={v}" for k, v in info.items()]
    lines += [f"membership_{k}={v}" for k, v in sorted(counts.items())]
    write_text(args.log, "\n".join(lines) + "\n")
    _say(lines[0] + (f" (pruned {n_pruned} reference tips)" if n_pruned else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
