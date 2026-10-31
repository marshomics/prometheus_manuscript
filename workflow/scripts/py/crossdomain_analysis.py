"""Cross-domain Prometheus phylogeny: origin and bacteria<->archaea transfer direction.

Reads the combined cross-domain ML gene tree (archaeal + bacterial validated Prometheus
proteins + the Prometheus references + the DNA2/NurA-HerA paralog outgroup; IQ-TREE
node labels 'SH-aLRT/UFBoot') and reports, from TOPOLOGY ALONE (no DTL model -- a
single rooted species tree spanning both domains is contested, so we do not
over-model), a parsimony reading of where the archaeal copies came from. Every number is
topology-based and descriptive; none is a significance test.

Rooting. The tree is rooted on the edge that separates the Paralog_ref tips from all
other tips (the Prometheus ingroup). If no edge does (a paralog nested inside the
ingroup, or candidates nested among the paralogs), the largest paralog-only clade is
used, `outgroup_monophyletic` is False and the direction is not called.

Character. Archaea = 'A', Bacteria = 'B' (validated bacterial loci; Prometheus_ref tips
count as Bacteria by default -- they are the characterised bacterial ProA -- or are
treated as missing with --ref-domain exclude). Paralog tips carry no state.

Metrics (crossdomain_summary.tsv; metric, value, note):
  * n_crossdomain_transitions: Fitch minimum number of A<->B changes on the
    outgroup-rooted ingroup; ..._supported: the same after collapsing branches below the
    support thresholds (soft polytomies; a lower bound robust to unsupported resolution);
  * n_archaeal_clades: maximal archaeal-only clades (= independent archaeal acquisitions
    IF the polarity is bacteria->archaea and there were no archaea->bacteria transfers);
    ..._supported: the minimum number of archaeal groups compatible with the supported
    topology; ..._nested_supported: clades whose attachment inside a bacterial grade is
    itself supported (see crossdomain_archaeal_clades.tsv);
  * ingroup_root_state: Fitch state set at the ingroup root of the ML tree;
  * inferred_direction: bacteria->archaea (or archaea->bacteria) only when ALL hold --
    the outgroup splits off cleanly, the ML root state is B (A), it stays B (A) on the
    support-collapsed tree, in >= --robust-fraction of the UFBoot replicate trees and of
    the balanced subsamples (the larger domain randomly thinned to the size of the
    smaller, floor --subsample-min, --n-subsample replicates, seeded), and -- without
    replicate trees -- the
    outgroup split itself is supported; otherwise 'ambiguous' with the failed criteria in
    the note. The balanced subsampling is the control for the bacteria:archaea sampling
    imbalance; midpoint rooting and excluding the references are reported as sensitivity.
Optional --codon adds a synonymous-saturation diagnostic for the cross-domain dN/dS
(NG86 pairwise dS within Archaea, within Bacteria and across domains).
Empty/trivial input (no tree, or one domain absent) -> header + status rows only.
"""
from __future__ import annotations
import argparse
import math
import os
import random
import re
import statistics
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SUMMARY_HEADER = "metric\tvalue\tnote"
TIPS_HEADER = "tip\tdomain\tsource\tstate\tarchaeal_clade"
CLADES_HEADER = ("clade_id\tn_tips\tstem_support\tparent_support\tsister_n_archaea\t"
                 "sister_n_bacteria\touter_n_archaea\touter_n_bacteria\tnested_in_bacteria\t"
                 "nested_supported\texample_tip")
SATURATION_HEADER = ("pair_class\tn_pairs\tn_saturated\tfrac_saturated\tmedian_dS\t"
                     "median_dN\tmedian_codons\tnote")
METHOD_NOTE = ("topology-based: Fitch parsimony on the outgroup-rooted ML gene tree; "
               "no DTL model, not a significance test")


class AnalysisError(Exception):
    pass


# --------------------------------------------------------------------------- #
# Newick -> undirected tree
# --------------------------------------------------------------------------- #
def parse_newick(txt):
    """Minimal newick parser. Returns (parent, children, label, length) lists; node 0 is
    the root. Handles quoted labels ('' escapes), [comments], internal labels."""
    s = txt.strip()
    parent, children, label, length = [None], [[]], [None], [None]

    def new(par):
        parent.append(par)
        children.append([])
        label.append(None)
        length.append(None)
        children[par].append(len(parent) - 1)
        return len(parent) - 1

    cur, i, n = 0, 0, len(s)
    delim = set("(),:;[")
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
        elif c == "[":
            j = s.find("]", i)
            if j < 0:
                raise AnalysisError("unterminated [comment] in newick")
            i = j + 1
        elif c == "(":
            cur = new(cur)
            i += 1
        elif c == ",":
            if parent[cur] is None:
                raise AnalysisError("malformed newick (',' at top level)")
            cur = new(parent[cur])
            i += 1
        elif c == ")":
            if parent[cur] is None:
                raise AnalysisError("malformed newick (unbalanced ')')")
            cur = parent[cur]
            i += 1
        elif c == ":":
            j = i + 1
            while j < n and s[j] not in delim and not s[j].isspace():
                j += 1
            try:
                length[cur] = float(s[i + 1:j])
            except ValueError:
                raise AnalysisError(f"bad branch length '{s[i + 1:j]}'")
            i = j
        elif c == ";":
            break
        elif c == "'":
            j, buf = i + 1, []
            while j < n:
                if s[j] == "'" and j + 1 < n and s[j + 1] == "'":
                    buf.append("'")
                    j += 2
                elif s[j] == "'":
                    break
                else:
                    buf.append(s[j])
                    j += 1
            label[cur] = "".join(buf)
            i = j + 1
        else:
            j = i
            while j < n and s[j] not in delim and not s[j].isspace():
                j += 1
            label[cur] = s[i:j]
            i = j
    return parent, children, label, length


def parse_support(lbl):
    """IQ-TREE 'alrt/ufboot' (or a single value) -> (alrt, ufboot); None when absent."""
    if lbl is None or str(lbl).strip() == "":
        return None
    parts = [p.split(":")[0] for p in str(lbl).split("/")]
    vals = []
    for p in parts:
        try:
            vals.append(float(p))
        except ValueError:
            return None
    if len(vals) >= 2:
        return (vals[-2], vals[-1])
    return (None, vals[0]) if vals else None


class UTree:
    """Undirected tree: adjacency {node: {nbr: length}}, leaf labels, edge supports."""

    def __init__(self, txt):
        parent, children, label, length = parse_newick(txt)
        self.adj = {i: {} for i in range(len(parent))}
        self.support = {}
        self.label = {}
        for v in range(len(parent)):
            p = parent[v]
            if p is not None:
                self.adj[v][p] = self.adj[p][v] = length[v] or 0.0
                if children[v]:
                    sup = parse_support(label[v])
                    if sup is not None:
                        self.support[frozenset((v, p))] = sup
            if not children[v]:
                if not label[v]:
                    raise AnalysisError("unlabelled leaf in tree")
                self.label[v] = label[v]
        root = 0
        if len(self.adj[root]) == 2 and root not in self.label:   # suppress a rooted bifurcation
            a, b = list(self.adj[root])
            # both root edges are ONE bipartition; keep the more conservative support
            sups = [s for s in (self.support.pop(frozenset((root, a)), None),
                                self.support.pop(frozenset((root, b)), None)) if s is not None]
            sup = min(sups, key=lambda s: (s[1] if s[1] is not None else 1e9,
                                           s[0] if s[0] is not None else 1e9)) if sups else None
            la, lb = self.adj[root][a], self.adj[root][b]
            del self.adj[a][root], self.adj[b][root], self.adj[root]
            self.adj[a][b] = self.adj[b][a] = la + lb
            if sup is not None:
                self.support[frozenset((a, b))] = sup
        self.leaves = sorted(self.label, key=lambda v: self.label[v])
        labels = [self.label[v] for v in self.leaves]
        if len(set(labels)) != len(labels):
            raise AnalysisError("duplicated tip labels in tree")
        self.bit = {v: 1 << k for k, v in enumerate(self.leaves)}
        self.full = (1 << len(self.leaves)) - 1
        self.by_label = {self.label[v]: v for v in self.leaves}

    def orient(self, root, away):
        """Children lists + postorder of the subtree hanging from `root` (parent `away`)."""
        par, kids, order, stack = {root: away}, {}, [], [root]
        while stack:
            x = stack.pop()
            order.append(x)
            kids[x] = [y for y in self.adj[x] if y != par[x]]
            for y in kids[x]:
                par[y] = x
                stack.append(y)
        order.reverse()
        return par, kids, order

    def edge_masks(self):
        """[(child, parent, leafmask_below_child)] for every edge (arbitrary rooting)."""
        start = next(v for v in self.adj if v not in self.label) if len(self.adj) > 1 else 0
        par, kids, order = self.orient(start, None)
        mask, out = {}, []
        for x in order:
            mask[x] = self.bit.get(x, 0)
            for y in kids[x]:
                mask[x] |= mask[y]
            if par[x] is not None:
                out.append((x, par[x], mask[x]))
        return out

    def mask_of(self, labels):
        m = 0
        for lb in labels:
            m |= self.bit[self.by_label[lb]]
        return m

    def sup(self, u, v):
        return self.support.get(frozenset((u, v)))


def outgroup_root(tree, out_labels):
    """(ingroup_root, away, monophyletic, intruders) for rooting on the outgroup."""
    omask = tree.mask_of(out_labels)
    imask = tree.full ^ omask
    best = None
    for c, p, m in tree.edge_masks():
        if m == imask:
            return c, p, True, []
        if m == omask:
            return p, c, True, []
        for side_mask, ing, aw in ((m, p, c), (tree.full ^ m, c, p)):
            if side_mask and not side_mask & imask:
                size = bin(side_mask).count("1")
                if best is None or size > best[0]:
                    best = (size, ing, aw, side_mask)
    if best is None:
        return None, None, False, []
    _, ing, aw, side_mask = best
    intr = [tree.label[v] for v in tree.leaves if tree.bit[v] & omask and not tree.bit[v] & side_mask]
    return ing, aw, False, intr


# --------------------------------------------------------------------------- #
# Parsimony on an oriented (sub)tree
# --------------------------------------------------------------------------- #
def combine(results):
    """Fitch rule over a node's children [(stateset, cost)], soft for polytomies (exact
    for a 2-state character): one change iff both states occur as singleton sets."""
    cost = sum(c for _, c in results)
    sets = [s for s, _ in results if s]
    if not sets:
        return frozenset(), cost
    singles = {next(iter(s)) for s in sets if len(s) == 1}
    if len(singles) == 2:
        return frozenset("AB"), cost + 1
    if len(singles) == 1:
        return frozenset(singles), cost
    return frozenset("AB"), cost


def oriented(tree, root, away, collapse_pred=None):
    """Oriented subtree (kids, postorder) with optional contraction of weak edges."""
    par, kids, order = tree.orient(root, away)
    if collapse_pred is None:
        return par, kids, order
    eff = {}
    for x in order:
        out, stack = [], list(kids[x])
        while stack:
            y = stack.pop()
            if kids[y] and collapse_pred(y, par[y]):
                stack.extend(kids[y])
            else:
                out.append(y)
        eff[x] = out
    # postorder over the effective tree
    order2, stack = [], [root]
    while stack:
        x = stack.pop()
        order2.append(x)
        stack.extend(eff[x])
    order2.reverse()
    return par, eff, order2


def fitch(tree, kids, order, states):
    res = {}
    for x in order:
        if not kids[x]:
            s = states.get(tree.label.get(x))
            res[x] = (frozenset(s) if s else frozenset(), 0)
        else:
            res[x] = combine([res[y] for y in kids[x]])
    return res


def counts(tree, kids, order, states):
    nA, nB = {}, {}
    for x in order:
        if not kids[x]:
            s = states.get(tree.label.get(x))
            nA[x], nB[x] = int(s == "A"), int(s == "B")
        else:
            nA[x] = sum(nA[y] for y in kids[x])
            nB[x] = sum(nB[y] for y in kids[x])
    return nA, nB


def state_str(s):
    return "/".join(sorted(s)) if s else "NA"


# --------------------------------------------------------------------------- #
def midpoint_state(tree, states):
    """Fitch root state with the tree midpoint-rooted on the character tips."""
    tips = [v for v in tree.leaves if states.get(tree.label[v])]
    if len(tips) < 2:
        return frozenset()

    def far(src):
        dist, prev, stack = {src: 0.0}, {src: None}, [src]
        while stack:
            x = stack.pop()
            for y, ln in tree.adj[x].items():
                if y not in dist:
                    dist[y], prev[y] = dist[x] + ln, x
                    stack.append(y)
        best = max(tips, key=lambda v: (dist[v], -v))
        return best, dist, prev

    a, _, _ = far(tips[0])
    b, dist, prev = far(a)
    half = dist[b] / 2.0
    x = b
    while prev[x] is not None and dist[prev[x]] > half:
        x = prev[x]
    y = prev[x]
    if y is None:
        return frozenset()
    sides = []
    for r, aw in ((x, y), (y, x)):
        _, kids, order = tree.orient(r, aw)
        sides.append(fitch(tree, kids, order, states)[r])
    return combine(sides)[0]


def analyse_tree(tree, states, out_labels):
    """Outgroup-rooted ML Fitch: returns dict or None when the tree has no outgroup."""
    if not out_labels:
        return None
    ing, aw, mono, intr = outgroup_root(tree, out_labels)
    if ing is None:
        return None
    _, kids, order = tree.orient(ing, aw)
    res = fitch(tree, kids, order, states)
    return {"root": ing, "away": aw, "mono": mono, "intruders": intr,
            "state": res[ing][0], "cost": res[ing][1], "kids": kids, "order": order}


# --------------------------------------------------------------------------- #
def saturation(codon_path, domain_of, n_pairs, seed, ds_max):
    """NG86 pairwise dN/dS by pair class; saturated = dS undefined (p_S >= 0.75) or > ds_max."""
    rows = []
    recs = [(n, s.upper()) for n, s in utils.read_fasta(codon_path)] if codon_path else []
    recs = [(n, s) for n, s in recs if domain_of.get(n) in ("Archaea", "Bacteria")]
    if len(recs) < 2:
        return [("all", 0, 0, "", "", "", "", "codon alignment empty or < 2 Prometheus "
                 "sequences; not assessable")]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds
    stops = {"TAA", "TAG", "TGA"}
    acgt = re.compile(r"^[ACGT]{3}$")
    rng = random.Random(seed)
    groups = {"Archaea": [r for r in recs if domain_of[r[0]] == "Archaea"],
              "Bacteria": [r for r in recs if domain_of[r[0]] == "Bacteria"]}
    classes = (("within_archaea", groups["Archaea"], groups["Archaea"]),
               ("within_bacteria", groups["Bacteria"], groups["Bacteria"]),
               ("cross_domain", groups["Archaea"], groups["Bacteria"]))
    for cname, g1, g2 in classes:
        if cname.startswith("within"):
            pairs = [(i, j) for i in range(len(g1)) for j in range(i + 1, len(g1))]
        else:
            pairs = [(i, j) for i in range(len(g1)) for j in range(len(g2))]
        if not pairs:
            rows.append((cname, 0, 0, "", "", "", "", "no pairs"))
            continue
        if len(pairs) > n_pairs:
            pairs = rng.sample(pairs, n_pairs)
        ds, dn, ncod, nsat, nused = [], [], [], 0, 0
        for i, j in pairs:
            s1, s2 = g1[i][1], g2[j][1]
            L = min(len(s1), len(s2)) // 3 * 3
            c1, c2 = [], []
            for k in range(0, L, 3):
                a, b = s1[k:k + 3], s2[k:k + 3]
                if acgt.match(a) and acgt.match(b) and a not in stops and b not in stops:
                    c1.append(a)
                    c2.append(b)
            if len(c1) < 50:
                continue
            nused += 1
            ncod.append(len(c1))
            try:
                d_n, d_s = cal_dn_ds(CodonSeq("".join(c1)), CodonSeq("".join(c2)), method="NG86")
            except (ValueError, ZeroDivisionError, KeyError):
                d_n, d_s = -1, -1
            if d_s is None or d_s < 0 or d_s > ds_max or not math.isfinite(d_s):
                nsat += 1
            else:
                ds.append(d_s)
                if d_n is not None and d_n >= 0 and math.isfinite(d_n):
                    dn.append(d_n)
        if nused == 0:
            rows.append((cname, 0, 0, "", "", "", "", "no pair with >= 50 shared codons"))
            continue
        rows.append((cname, nused, nsat, round(nsat / nused, 4),
                     round(statistics.median(ds), 4) if ds else "",
                     round(statistics.median(dn), 4) if dn else "",
                     int(statistics.median(ncod)),
                     f"NG86 (Jukes-Cantor); saturated = dS undefined or > {ds_max}; median dS/dN "
                     f"over unsaturated pairs"))
    return rows


# --------------------------------------------------------------------------- #
def load_id_set(path):
    return set(utils.fasta_ids(path)) if path and os.path.exists(path) else set()


def load_ref_labels(path):
    labels = {}
    if not path or not os.path.exists(path):
        return labels
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1].strip() in ("Prometheus_ref", "Paralog_ref"):
                labels[parts[0].strip()] = parts[1].strip()
    return labels


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--arch-proteins", required=True, help="archaeal validated proteins FASTA")
    ap.add_argument("--bact-proteins", required=True, help="bacterial validated proteins FASTA")
    ap.add_argument("--ref-labels", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-tips", required=True)
    ap.add_argument("--out-clades", default="")
    ap.add_argument("--out-rooted-tree", default="")
    ap.add_argument("--ufboot", default="", help="IQ-TREE UFBoot replicate trees (--wbtl)")
    ap.add_argument("--max-boot-trees", type=int, default=1000)
    ap.add_argument("--min-ufboot", type=float, default=95.0)
    ap.add_argument("--min-alrt", type=float, default=80.0)
    ap.add_argument("--n-subsample", type=int, default=200)
    ap.add_argument("--subsample-min", type=int, default=5,
                    help="the larger domain is thinned to max(size of the smaller, this)")
    ap.add_argument("--robust-fraction", type=float, default=0.95)
    ap.add_argument("--ref-domain", choices=["Bacteria", "exclude"], default="Bacteria")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--codon", default="", help="cross-domain codon alignment (saturation check)")
    ap.add_argument("--out-saturation", default="")
    ap.add_argument("--n-pairs", type=int, default=200)
    ap.add_argument("--saturation-ds", type=float, default=3.0)
    args = ap.parse_args()
    try:
        run(args)
    except AnalysisError as exc:
        utils.log(f"FATAL: {exc}")
        sys.exit(1)


def write_rows(path, header, rows):
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v) for v in r) + "\n")


def run(args):
    for p in (args.out, args.out_tips, args.out_clades, args.out_rooted_tree, args.out_saturation):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    arch_ids, bact_ids = load_id_set(args.arch_proteins), load_id_set(args.bact_proteins)
    ref_labels = load_ref_labels(args.ref_labels)
    both = arch_ids & bact_ids
    if both:
        raise AnalysisError(f"ids present in both the archaeal and bacterial sets: {sorted(both)[:3]}")
    domain_of = {i: "Archaea" for i in arch_ids}
    domain_of.update({i: "Bacteria" for i in bact_ids})

    sat_rows = None
    if args.out_saturation:
        sat_rows = saturation(args.codon, domain_of, args.n_pairs, args.seed, args.saturation_ds)
        write_rows(args.out_saturation, SATURATION_HEADER, sat_rows)

    def bail(note):
        rows = [("status", "incomplete", note), ("direction_method", "topology-based", METHOD_NOTE)]
        rows += sat_summary(sat_rows)
        write_rows(args.out, SUMMARY_HEADER, rows)
        write_rows(args.out_tips, TIPS_HEADER, [])
        if args.out_clades:
            write_rows(args.out_clades, CLADES_HEADER, [])
        if args.out_rooted_tree:
            with open(args.out_rooted_tree, "w") as fh:
                fh.write("();\n")
        utils.log(f"crossdomain_analysis: {note}")

    txt = open(args.tree).read().strip() if os.path.exists(args.tree) else ""
    if not txt or txt in ("();", ";") or txt.count(",") < 2:
        bail("trivial/empty cross-domain tree (fewer than 4 sequences or one domain absent)")
        return
    tree = UTree(txt)

    # ---- tip classes ------------------------------------------------------------ #
    tips, states, out_labels, unknown = [], {}, [], []
    for v in tree.leaves:
        lb = tree.label[v]
        if ref_labels.get(lb) == "Paralog_ref":
            dom, src = "Paralog", "paralog_reference"
            out_labels.append(lb)
        elif ref_labels.get(lb) == "Prometheus_ref":
            src = "prometheus_reference"
            dom = "Bacteria" if args.ref_domain == "Bacteria" else "Reference"
            if dom == "Bacteria":
                states[lb] = "B"
        elif lb in domain_of:
            dom = domain_of[lb]
            src = "archaeal_validated" if dom == "Archaea" else "bacterial_validated"
            states[lb] = "A" if dom == "Archaea" else "B"
        else:
            dom, src = "Unknown", "unknown"
            unknown.append(lb)
        tips.append([lb, dom, src, states.get(lb, "")])
    if unknown:
        raise AnalysisError(f"{len(unknown)} tree tip(s) are in none of the input sets (ids "
                            f"renamed by the aligner/tree builder?): {unknown[:3]}")
    n_arch = sum(1 for s in states.values() if s == "A")
    n_bact = sum(1 for s in states.values() if s == "B")
    n_ref = sum(1 for t in tips if t[2] == "prometheus_reference")
    if n_arch == 0 or n_bact == 0:
        bail(f"need both domains in the tree (archaea={n_arch}, bacteria={n_bact})")
        return

    has_sup = bool(tree.support)

    def weak(u, v):
        s = tree.sup(u, v)
        if s is None:
            return True
        alrt, ufb = s
        ok = (alrt is None or alrt >= args.min_alrt) and (ufb is None or ufb >= args.min_ufboot)
        return not ok

    rows = [("status", "complete", "")]
    ml = analyse_tree(tree, states, out_labels)
    if ml is None:
        rooting = "midpoint(no_paralog_outgroup)"
        mid = midpoint_state(tree, states)
        rows += [("n_archaeal_tips", n_arch, "validated archaeal Prometheus on the tree"),
                 ("n_bacterial_tips", n_bact, "validated bacterial Prometheus (+ references "
                                              "unless --ref-domain exclude)"),
                 ("rooting", rooting, "no Paralog_ref tips: the direction cannot be polarised"),
                 ("midpoint_root_state", state_str(mid), "sensitivity only"),
                 ("ingroup_root_state", "NA", "no outgroup"),
                 ("inferred_direction", "ambiguous", "no paralog outgroup on the tree"),
                 ("direction_method", "topology-based", METHOD_NOTE)]
        rows += sat_summary(sat_rows)
        write_rows(args.out, SUMMARY_HEADER, rows)
        write_rows(args.out_tips, TIPS_HEADER, [t + [""] for t in tips])
        if args.out_clades:
            write_rows(args.out_clades, CLADES_HEADER, [])
        if args.out_rooted_tree:
            with open(args.out_rooted_tree, "w") as fh:
                fh.write(txt + ("\n" if txt.endswith(";") else ";\n"))
        return

    ing, aw, kids, order = ml["root"], ml["away"], ml["kids"], ml["order"]
    par = {ing: aw}
    for x in order:
        for y in kids[x]:
            par[y] = x
    nA, nB = counts(tree, kids, order, states)

    # ---- maximal archaeal clades on the ML tree ---------------------------------- #
    clades, clade_of_tip = [], {}
    for x in reversed(order):            # preorder
        if nA[x] > 0 and nB[x] == 0 and (x == ing or nB[par[x]] > 0):
            clades.append(x)
    clade_rows = []
    for k, x in enumerate(clades, 1):
        cid = f"AC{k}"
        stack, ex = [x], ""
        while stack:
            y = stack.pop()
            if not kids[y]:
                if states.get(tree.label[y]) == "A":
                    clade_of_tip[tree.label[y]] = cid
                    ex = ex or tree.label[y]
            stack.extend(kids[y])
        stem = tree.sup(x, par[x]) if kids[x] and x != ing else None
        if x == ing:
            row = [cid, nA[x], "NA", "NA", 0, 0, 0, 0, False, False, ex]
        else:
            p = par[x]
            sis = [y for y in kids[p] if y != x]
            sa, sb = sum(nA[y] for y in sis), sum(nB[y] for y in sis)
            if p == ing:
                oa = ob = 0
                psup = None
            else:
                q = par[p]
                outer = [y for y in kids[q] if y != p]
                oa, ob = sum(nA[y] for y in outer), sum(nB[y] for y in outer)
                psup = tree.sup(p, q)
            nested = sb > 0 and ob > 0
            nested_sup = bool(nested and p != ing and not weak(p, par[p]))
            fmt = (lambda s: "NA" if s is None else "/".join("" if v is None else f"{v:g}" for v in s))
            row = [cid, nA[x], fmt(stem) if kids[x] else "tip", fmt(psup), sa, sb, oa, ob,
                   nested, nested_sup, ex]
        clade_rows.append(row)
    n_clades = len(clades)
    n_nested_sup = sum(1 for r in clade_rows if r[9] is True)

    # ---- support-collapsed tree --------------------------------------------------- #
    if has_sup:
        _, ekids, eorder = oriented(tree, ing, aw, collapse_pred=weak)
        cres = fitch(tree, ekids, eorder, states)
        c_state, c_cost = cres[ing]
        cA, cB = counts(tree, ekids, eorder, states)
        n_groups = sum(1 for x in eorder if ekids[x] and cB[x] > 0 and
                       any(cA[y] > 0 and cB[y] == 0 for y in ekids[x]))
        if cB[ing] == 0:
            n_groups = 1
        out_sup = tree.sup(ing, aw)
    else:
        c_state, c_cost, n_groups, out_sup = frozenset(), None, None, None

    # ---- UFBoot replicate trees ---------------------------------------------------- #
    boot = {"n": 0, "A": 0, "B": 0, "AB": 0, "unrooted": 0}
    if args.ufboot and os.path.exists(args.ufboot) and os.path.getsize(args.ufboot) > 0:
        with open(args.ufboot) as fh:
            for line in fh:
                line = line.strip()
                if not line or boot["n"] >= args.max_boot_trees:
                    continue
                bt = UTree(line)
                if set(bt.by_label) != set(tree.by_label):
                    raise AnalysisError("UFBoot replicate tree tips differ from the ML tree")
                boot["n"] += 1
                r = analyse_tree(bt, states, out_labels)
                if r is None or not r["mono"]:
                    boot["unrooted"] += 1
                    continue
                boot["".join(sorted(r["state"])) or "AB"] += 1

    # ---- balanced subsampling ----------------------------------------------------- #
    rng = random.Random(args.seed)
    a_tips = sorted(lb for lb, s in states.items() if s == "A")
    b_tips = sorted(lb for lb, s in states.items() if s == "B")
    big, small = (b_tips, a_tips) if len(b_tips) >= len(a_tips) else (a_tips, b_tips)
    # thin the larger domain to the size of the smaller one, but never below
    # --subsample-min tips (a 1:1 sample of one archaeal tip cannot show nesting)
    k = min(len(big), max(len(small), args.subsample_min))
    sub = {"A": 0, "B": 0, "AB": 0}
    for _ in range(args.n_subsample):
        keep = set(small) | set(rng.sample(big, k))
        st = {lb: s for lb, s in states.items() if lb in keep}
        s = fitch(tree, kids, order, st)[ing][0]
        sub["".join(sorted(s)) or "AB"] += 1

    # ---- sensitivity: midpoint rooting, references excluded ------------------------ #
    mid = midpoint_state(tree, states)
    src_of = {t[0]: t[2] for t in tips}
    st_noref = {lb: s for lb, s in states.items() if src_of[lb] != "prometheus_reference"}
    noref_state = fitch(tree, kids, order, st_noref)[ing][0] if (
        any(s == "B" for s in st_noref.values())) else frozenset()

    # ---- direction call ----------------------------------------------------------- #
    def call(target):
        fails = []
        if not ml["mono"]:
            fails.append("outgroup not monophyletic")
        if ml["state"] != frozenset(target):
            fails.append(f"ML root state {state_str(ml['state'])}")
        if has_sup and c_state != frozenset(target):
            fails.append(f"collapsed-tree root state {state_str(c_state)}")
        if boot["n"]:
            if boot[target] / boot["n"] < args.robust_fraction:
                fails.append(f"UFBoot root {target} in {boot[target]}/{boot['n']}")
        elif has_sup and (out_sup is None or weak(ing, aw)):
            fails.append("outgroup split unsupported (no replicate trees)")
        elif not has_sup:
            fails.append("no support values on the tree")
        if args.n_subsample and sub[target] / args.n_subsample < args.robust_fraction:
            fails.append(f"balanced subsamples root {target} in {sub[target]}/{args.n_subsample}")
        return fails

    fails_b, fails_a = call("B"), call("A")
    if not fails_b:
        direction, dnote = "bacteria->archaea", "all robustness criteria met (bacterial root)"
    elif not fails_a:
        direction, dnote = "archaea->bacteria", "all robustness criteria met (archaeal root)"
    else:
        direction = "ambiguous"
        dnote = "not robust: " + "; ".join(fails_b if ml["state"] != frozenset("A") else fails_a)
    ml_dir = {frozenset("B"): "bacteria->archaea", frozenset("A"): "archaea->bacteria"}.get(
        ml["state"], "ambiguous")

    fmt_sup = "NA" if out_sup is None else "/".join("" if v is None else f"{v:g}" for v in out_sup)
    rooting = "outgroup" if ml["mono"] else "outgroup_nonmonophyletic(largest paralog clade)"
    rows += [
        ("n_archaeal_tips", n_arch, "validated archaeal Prometheus on the cross-domain tree"),
        ("n_bacterial_tips", n_bact, "validated bacterial Prometheus"
         + (f" + {n_ref} Prometheus reference(s)" if args.ref_domain == "Bacteria" else "")),
        ("n_reference_tips", n_ref, f"Prometheus_ref tips (ref-domain={args.ref_domain})"),
        ("n_paralog_tips", len(out_labels), "DNA2/NurA-HerA outgroup tips"),
        ("rooting", rooting, (f"intruding paralog tips: {ml['intruders'][:5]}"
                              if ml["intruders"] else "ingroup|paralog split")),
        ("outgroup_monophyletic", ml["mono"], "an edge separates every paralog from the ingroup"),
        ("outgroup_split_support", fmt_sup, "SH-aLRT/UFBoot of the ingroup|outgroup edge"),
        ("n_crossdomain_transitions", ml["cost"],
         "Fitch minimum Bacteria<->Archaea changes on the outgroup-rooted ingroup (ML tree)"),
        ("n_crossdomain_transitions_supported", "NA" if c_cost is None else c_cost,
         f"same with branches below SH-aLRT {args.min_alrt:g}/UFBoot {args.min_ufboot:g} "
         f"collapsed (lower bound)"),
        ("n_archaeal_clades", n_clades,
         "maximal archaeal-only clades = independent archaeal acquisitions if the polarity is "
         "bacteria->archaea without back-transfer"),
        ("n_archaeal_clades_supported", "NA" if n_groups is None else n_groups,
         "minimum archaeal groups compatible with the supported topology"),
        ("n_archaeal_clades_nested_supported", n_nested_sup,
         "archaeal clades whose placement inside a bacterial grade is supported"),
        ("archaeal_monophyletic", n_clades == 1,
         "True = a single archaeal clade on the ML tree"),
        ("ingroup_root_state", state_str(ml["state"]), "Fitch state set at the ingroup root (ML)"),
        ("collapsed_root_state", state_str(c_state) if has_sup else "NA",
         "Fitch root state after collapsing unsupported branches"),
        ("n_bootstrap_trees", boot["n"], "UFBoot replicate trees evaluated"),
        ("bootstrap_root_B_fraction", round(boot["B"] / boot["n"], 4) if boot["n"] else "NA",
         f"replicates with a bacterial ingroup root ({boot['unrooted']} without a clean "
         f"outgroup split count as unresolved)"),
        ("bootstrap_root_A_fraction", round(boot["A"] / boot["n"], 4) if boot["n"] else "NA",
         "replicates with an archaeal ingroup root"),
        ("n_subsample", args.n_subsample,
         f"subsamples: all {len(small)} tips of the smaller domain + {k} random tips of the "
         f"larger (seed {args.seed})"),
        ("subsample_root_B_fraction", round(sub["B"] / args.n_subsample, 4) if args.n_subsample else "NA",
         "controls the bacteria:archaea sampling imbalance"),
        ("subsample_root_A_fraction", round(sub["A"] / args.n_subsample, 4) if args.n_subsample else "NA",
         ""),
        ("midpoint_root_state", state_str(mid), "sensitivity: midpoint rooting of the character tips"),
        ("root_state_excluding_references", state_str(noref_state),
         "sensitivity: Prometheus_ref tips treated as missing"),
        ("inferred_direction_ml", ml_dir, "ML-tree root state alone (not robustness-checked)"),
        ("inferred_direction", direction, dnote),
        ("direction_method", "topology-based", METHOD_NOTE),
    ]
    rows += sat_summary(sat_rows)
    write_rows(args.out, SUMMARY_HEADER, rows)
    write_rows(args.out_tips, TIPS_HEADER, [t + [clade_of_tip.get(t[0], "")] for t in tips])
    if args.out_clades:
        write_rows(args.out_clades, CLADES_HEADER, clade_rows)
    if args.out_rooted_tree:
        with open(args.out_rooted_tree, "w") as fh:
            fh.write(rooted_newick(tree, ing, aw) + "\n")
    utils.log(f"crossdomain: archaea={n_arch} bacteria={n_bact} transitions={ml['cost']} "
              f"archaeal_clades={n_clades} root={state_str(ml['state'])} direction={direction}")


def sat_summary(sat_rows):
    if not sat_rows:
        return []
    cross = next((r for r in sat_rows if r[0] == "cross_domain"), None)
    if cross is None or not cross[1]:
        return [("dnds_saturation_warning", "NA", "saturation not assessable (no cross-domain "
                 "codon pairs)")]
    frac = cross[3]
    return [("cross_domain_dS_saturated_fraction", frac,
             "NG86 pairs with undefined or extreme dS (crossdomain_saturation.tsv)"),
            ("cross_domain_median_dS", cross[4] if cross[4] != "" else "NA",
             "median NG86 dS of unsaturated cross-domain pairs"),
            ("dnds_saturation_warning", bool(frac > 0.5),
             "True = synonymous sites saturated for most cross-domain pairs: treat the "
             "cross-domain FEL/MEME/BUSTED/aBSREL results as unreliable")]


def rooted_newick(tree, ing, aw):
    """Newick of the outgroup-rooted tree (root on the ingroup|outgroup edge)."""
    def fmt_sup(u, v):
        s = tree.sup(u, v)
        return "" if s is None else "/".join("" if x is None else f"{x:g}" for x in s)

    def write(root, away):
        par, kids, order = tree.orient(root, away)
        out = {}
        for x in order:
            ln = tree.adj[x][par[x]] if par[x] is not None else 0.0
            if par[x] == away:
                ln = tree.adj[root][away] / 2.0
            if not kids[x]:
                out[x] = f"{tree.label[x]}:{ln:.8g}"
            else:
                lab = fmt_sup(x, par[x]) if par[x] != away else ""
                out[x] = "(" + ",".join(out[y] for y in kids[x]) + f"){lab}:{ln:.8g}"
        return out[root]
    return "(" + write(ing, aw) + "," + write(aw, ing) + ");"


if __name__ == "__main__":
    main()
