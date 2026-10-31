"""Per-domain alignments + ML trees, compared with the whole-protein (family) tree.

Stage 5.2 (domain shuffling / fusion / independent acquisition). For each configured
domain category (utils.load_domain_spec; exact Pfam matching) every protein carrying it
contributes the sub-sequence spanning that domain's hmmscan ENVELOPE (utils.domain_intervals;
pieces of a split domain, e.g. AAA_11 + AAA_12 of the SF1 helicase, are spanned together
unless another category lies between them, in which case the longest same-category block
is used and the protein is counted as 'duplicated/rearranged').

Per domain (only proteins also present in the family tree; >= --min-taxa of them):
    align (MAFFT L-INS-i <= --linsi-max seqs, else FAMSA; or FAMSA throughout)
    -> ClipKIT (kpic-smart-gap) -> IQ-TREE 2 (--model, UFBoot2 -B 1000, --seed)
    -> AU test (iqtree2 -z: the domain ML tree vs the family topology pruned to the same
       taxa, -zb 10000 -au, model fixed to the one the ML run used)

Topology comparison (dendropy, both trees force-unrooted, restricted to shared taxa):
    rf_to_family     unweighted Robinson-Foulds distance (symmetric difference of splits)
    normalized_rf    RF / (2 (n - 3)), the maximum for two binary unrooted n-taxon trees
    n_supported_conflicts  domain-tree splits with UFBoot >= 95 that are INCOMPATIBLE with a
                     family-tree split with support >= 95 (well-supported discordance, as
                     opposed to RF noise from unresolved branches)
    au_pvalue_family_topology  AU-test p that the domain alignment fits the family topology
                     as well as its own ML tree (small p = the domain rejects forced congruence)

NO SILENT FAILURE: a missing or crashing mafft/famsa/clipkit/iqtree2 raises. Rows with a
note and no statistic are written only for genuinely insufficient data (domain absent,
fewer than --min-taxa carriers, trivial family tree).
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

COLS = ["domain", "n_tips", "rf_to_family", "normalized_rf", "note", "n_seqs", "n_splits_domain",
        "n_supported_conflicts", "au_pvalue_family_topology", "model", "n_rearranged"]
HEADER = "\t".join(COLS)
SUPPORT = 95.0


def fail(msg):
    raise SystemExit(f"[domain_trees] {msg}")


def tool(*names):
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    fail(f"required tool not on PATH: {' / '.join(names)}")


def tool_env():
    """Environment for external tools with an ABSOLUTE, existing TMPDIR: MAFFT cd's into its
    mktemp directory, so a relative TMPDIR (e.g. './.tmp') is not safe."""
    e = dict(os.environ)
    tmp = e.get("TMPDIR") or ""
    if not os.path.isabs(tmp):        # relative or unset: never create it relative to the cwd
        tmp = "/tmp"
    os.makedirs(tmp, exist_ok=True)
    e["TMPDIR"] = tmp
    return e


def run(cmd, tag, stdout_path=None):
    utils.log(f"[{tag}] " + " ".join(str(c) for c in cmd))
    if stdout_path:
        with open(stdout_path, "w") as out:
            p = subprocess.run([str(c) for c in cmd], stdout=out, stderr=subprocess.PIPE, text=True,
                               env=tool_env())
    else:
        p = subprocess.run([str(c) for c in cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                           env=tool_env())
    if p.returncode != 0:
        fail(f"{tag} failed (exit {p.returncode}): {(p.stderr or '')[-1500:]}")
    return p


# --------------------------------------------------------------------------- #
# domain extraction
# --------------------------------------------------------------------------- #
def domain_span(hits, spec, cat):
    """(start, end, rearranged) of category `cat` in one protein, or None.

    Envelope coordinates are used. Intervals of `cat` separated by another category's
    interval are treated as separate blocks (duplication/rearrangement); the block with
    the largest covered length is returned.
    """
    env = [dict(h, start=h.get("env_from", h["start"]), end=h.get("env_to", h["end"])) for h in hits]
    ivs = utils.domain_intervals(env, spec)
    mine = sorted(ivs.get(cat, []))
    if not mine:
        return None
    others = sorted(iv for c, v in ivs.items() if c != cat for iv in v)
    blocks = [[mine[0]]]
    for iv in mine[1:]:
        prev_end = blocks[-1][-1][1]
        if any(o[0] > prev_end and o[1] < iv[0] for o in others):
            blocks.append([iv])
        else:
            blocks[-1].append(iv)
    best = max(blocks, key=lambda b: sum(e - s + 1 for s, e in b))
    return best[0][0], max(e for _, e in best), len(blocks) > 1


def extract(proteins, hits, spec, cat):
    recs, n_rearr = [], 0
    for pid, seq in proteins.items():
        sp = domain_span(hits.get(pid, []), spec, cat)
        if sp is None:
            continue
        s, e, rearr = sp
        s, e = max(1, s), min(len(seq), e)
        if e - s + 1 < 10:
            continue
        recs.append((pid, seq[s - 1:e]))
        n_rearr += int(rearr)
    return recs, n_rearr


# --------------------------------------------------------------------------- #
# trees
# --------------------------------------------------------------------------- #
def read_tree(path, tns=None):
    import dendropy
    return dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True,
                             rooting="force-unrooted",
                             taxon_namespace=tns if tns is not None else dendropy.TaxonNamespace())


def tree_tips(path):
    return {lf.taxon.label for lf in read_tree(path).leaf_node_iter()}


def support_of(node):
    """UFBoot support from an IQ-TREE label ('95' or 'SH-aLRT/UFBoot' = '87.5/95')."""
    lab = node.label
    if lab is None:
        return None
    parts = str(lab).split("/")
    try:
        return float(parts[-1])
    except ValueError:
        return None


def splits_with_support(tree):
    """{frozenset(taxon labels on the side NOT containing the first taxon): support}."""
    present = {lf.taxon.label for lf in tree.leaf_node_iter()}
    ref = min(present)
    out = {}
    for nd in tree.postorder_internal_node_iter():
        if nd.parent_node is None:
            continue
        side = {lf.taxon.label for lf in nd.leaf_iter()}
        if ref in side:
            side = present - side
        if 1 < len(side) < len(present) - 1:
            out[frozenset(side)] = support_of(nd)
    return out


def compatible(a, b, universe):
    """Two splits are compatible iff one of the four intersections is empty."""
    ac, bc = universe - a, universe - b
    return not (a & b) or not (a & bc) or not (ac & b) or not (ac & bc)


def compare(domain_tree, family_tree, shared):
    """(n, rf, normalized_rf, n_splits_domain, n_supported_conflicts)."""
    import dendropy
    from dendropy.calculate import treecompare
    tns = dendropy.TaxonNamespace()
    dt = read_tree(domain_tree, tns)
    ft = read_tree(family_tree, tns)
    for t in (dt, ft):
        t.retain_taxa_with_labels(shared)
    n = len(shared)
    rf = treecompare.symmetric_difference(dt, ft)
    norm = rf / (2 * (n - 3)) if n > 3 else float("nan")
    ds, fs = splits_with_support(dt), splits_with_support(ft)
    uni = frozenset(shared)
    strong_f = [s for s, v in fs.items() if v is not None and v >= SUPPORT]
    n_conf = sum(1 for s, v in ds.items()
                 if v is not None and v >= SUPPORT and any(not compatible(s, f, uni) for f in strong_f))
    return n, rf, norm, len(ds), n_conf


def iqtree_model(prefix):
    """Model the IQ-TREE run actually used (from <prefix>.iqtree)."""
    rep = prefix + ".iqtree"
    if os.path.exists(rep):
        with open(rep) as fh:
            for line in fh:
                m = re.match(r"\s*Model of substitution:\s*(\S+)", line)
                if m:
                    return m.group(1)
    fail(f"could not read the substitution model from {rep}")


def au_pvalue(iqtree, aln, model, ml_tree, family_pruned, prefix, threads, seed):
    cand = prefix + ".candidates.nwk"
    with open(cand, "w") as out:
        for p in (ml_tree, family_pruned):
            with open(p) as fh:
                out.write(fh.read().strip() + "\n")
    run([iqtree, "-s", aln, "-m", model, "-z", cand, "-n", "0", "-zb", "10000", "-au",
         "-T", threads, "--seed", seed, "--prefix", prefix, "-redo", "-quiet"], "iqtree2 AU")
    rep = prefix + ".iqtree"
    rows, in_tab = {}, False
    with open(rep) as fh:
        for line in fh:
            if line.startswith("Tree      logL"):
                in_tab = True
                continue
            if in_tab:
                f = line.split()
                if not f:
                    if rows:
                        break
                    continue
                if f[0].isdigit():
                    # columns: Tree logL deltaL bp-RELL [+-] p-KH [+-] p-SH [+-] c-ELW [+-] p-AU [+-]
                    nums = [x for x in f[1:] if x not in ("+", "-")]
                    rows[int(f[0])] = float(nums[-1])
    if 2 not in rows:
        fail(f"could not parse the AU-test table in {rep}")
    return rows[2]


def write_pruned(family_tree, keep, path):
    t = read_tree(family_tree)
    t.retain_taxa_with_labels(set(keep))
    if len(t.seed_node.child_nodes()) == 2 and len(keep) > 2:
        t.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    s = t.as_string(schema="newick", suppress_rooting=True, suppress_internal_node_labels=True,
                    suppress_edge_lengths=False, unquoted_underscores=True)
    with open(path, "w") as out:
        out.write(s.strip() + "\n")


def align(records, out_aln, tool_name, linsi_max, threads):
    faa = out_aln + ".input.faa"
    utils.write_fasta(records, faa)
    if tool_name == "famsa" or (tool_name == "mafft" and len(records) > linsi_max):
        run([tool("famsa"), "-t", threads, faa, out_aln], "famsa")
        return "famsa"
    # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
    run([tool("mafft"), "--amino", "--localpair", "--maxiterate", "1000", "--thread", threads, "--threadit", "0", faa],
        "mafft L-INS-i", stdout_path=out_aln)
    return "mafft-linsi"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--family-tree", required=True)
    ap.add_argument("--domains", required=True, help="prometheus.domains JSON (utils.load_domain_spec)")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-taxa", type=int, default=6,
                    help="minimum carriers shared with the family tree (n-3 >= 3 splits)")
    ap.add_argument("--model", default="MFP", help="IQ-TREE -m for domain trees (default ModelFinder)")
    ap.add_argument("--ufboot", type=int, default=1000)
    ap.add_argument("--align-tool", default="mafft", choices=("mafft", "famsa"))
    ap.add_argument("--linsi-max", type=int, default=200)
    ap.add_argument("--trim", default="clipkit", choices=("clipkit", "none"))
    ap.add_argument("--clipkit-mode", default="kpic-smart-gap")
    ap.add_argument("--au-test", default="True")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    spec = utils.load_domain_spec(args.domains)
    if not spec:
        fail("--domains is empty")
    proteins = {n: s.rstrip("*") for n, s in utils.read_fasta(args.proteins)}
    hits = utils.parse_domtbl(args.domtbl)

    rows = []

    def row(domain, **kw):
        r = {c: "" for c in COLS}
        r["domain"] = domain
        r.update(kw)
        rows.append(r)

    fam_txt = open(args.family_tree).read().strip() if os.path.exists(args.family_tree) else ""
    family_ok = fam_txt not in ("", "();", ";")
    fam_tips = tree_tips(args.family_tree) if family_ok else set()

    for cat in utils.canonical_order(spec):
        recs, n_rearr = extract(proteins, hits, spec, cat)
        n_seqs = len(recs)
        if not family_ok:
            row(cat, n_tips=0, n_seqs=n_seqs, note="not computed: family tree trivial/empty")
            continue
        recs = [(n, s) for n, s in recs if n in fam_tips]
        n_tips = len(recs)
        if n_tips < args.min_taxa:
            row(cat, n_tips=n_tips, n_seqs=n_seqs, n_rearranged=n_rearr,
                note=f"not computed: {n_tips} carriers in the family tree (< {args.min_taxa})")
            continue
        aln = os.path.join(args.workdir, f"{cat}.aln")
        method = align(recs, aln, args.align_tool, args.linsi_max, args.threads)
        tree_in = aln
        if args.trim == "clipkit":
            trim = os.path.join(args.workdir, f"{cat}.trim.aln")
            run([tool("clipkit"), aln, "-m", args.clipkit_mode, "-o", trim], "clipkit")
            tree_in = trim
        prefix = os.path.join(args.workdir, cat)
        cmd = [tool("iqtree2", "iqtree"), "-s", tree_in, "-m", args.model, "-T", args.threads,
               "--seed", args.seed, "--prefix", prefix, "-redo", "-quiet"]
        if args.ufboot:
            cmd += ["-B", args.ufboot]
        run(cmd, f"iqtree2 {cat}")
        treefile = prefix + ".treefile"
        if not os.path.exists(treefile) or os.path.getsize(treefile) == 0:
            fail(f"iqtree2 wrote no tree for {cat}")
        model = iqtree_model(prefix)
        shared = tree_tips(treefile) & fam_tips
        n, rf, norm, n_spl, n_conf = compare(treefile, args.family_tree, shared)
        au = ""
        if str(args.au_test).lower() in ("true", "1", "yes"):
            fam_pruned = prefix + ".family_pruned.nwk"
            write_pruned(args.family_tree, shared, fam_pruned)
            au = round(au_pvalue(tool("iqtree2", "iqtree"), tree_in, model, treefile, fam_pruned,
                                 prefix + ".au", args.threads, args.seed), 6)
        if not args.ufboot:
            n_conf = ""       # no UFBoot support on the domain tree: not computable
        row(cat, n_tips=n, rf_to_family=rf, normalized_rf=round(norm, 6), n_seqs=n_seqs,
            n_splits_domain=n_spl, n_supported_conflicts=n_conf, au_pvalue_family_topology=au,
            model=model, n_rearranged=n_rearr,
            note=f"ok ({method}; {args.trim}; UFBoot>={SUPPORT:g} for conflicts)")

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in COLS) + "\n")
    utils.log(f"wrote {len(rows)} domain row(s) to {args.out}")


if __name__ == "__main__":
    main()
