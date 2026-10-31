"""Selection analysis on the Prometheus codon alignment (HyPhy FEL/MEME/BUSTED/aBSREL/RELAX).

Stage 5.6. Diversifying selection concentrated on the nuclease/REase domains would be the
molecular signature of host-virus conflict; RELAX asks whether selection intensity differs
on methanogen lineages.

Verified against HyPhy 2.5.62 (the pinned version) by running every method on a small
codon alignment (CLI ``hyphy ENV='SetParameter(RANDOM_SEED,<s>,0);' CPU=<n> <method>
--alignment A --tree T --output J``):
  * FEL   MLE.headers: alpha, beta, alpha=beta, LRT, p-value, Total branch length
  * MEME  MLE.headers: &alpha;, &beta;<sup>1</sup>, p<sup>1</sup>, &beta;<sup>+</sup>,
          p<sup>+</sup>, LRT, p-value, # branches under selection, ... (2.5.103 inserts
          'LRT MEME vs FEL'), so columns are matched by NORMALISED NAME, never by index
  * BUSTED 'test results': {LRT, p-value};  aBSREL 'test results': {P-value threshold,
          positive test results, tested};  RELAX 'test results': {LRT, p-value,
          relaxation or intensification parameter}
  * 'data partitions'[i].coverage lists the 0-based codon sites of each partition, used
    to map MLE rows back to alignment columns (also for GARD-partitioned input)
  * RELAX classic mode: branches labelled {Test} / {Reference}; every other branch is put
    in an 'Unclassified' nuisance set with its own omega distribution.
  * HyPhy seeds its RNG from time+pid (global_things.cpp GlobalStartup); BUSTED's random
    initial grid made two unseeded 2.5.62 runs on the same data end at different
    log-likelihoods. The ENV= statement SetParameter(RANDOM_SEED, s, 0) re-seeds it
    (batchlanruntime.cpp), and seeded runs are bit-identical, so --seed (default 12345) is
    always passed.

Behaviour:
  * the analysed set = codon records  ∩  gene-tree tips  [∩ catalytic_status == intact];
    the gene tree is ALWAYS pruned to that set and written without internal node labels;
  * RELAX labelling (--relax-labelling): 'clade' (default) labels every branch whose
    descendant tips are all methanogens Test and all non-methanogens Reference (mixed
    clades and ANME/unknown hosts Unclassified); 'terminal' labels tip branches only;
  * site table: one row per codon column; FEL (beta/alpha, p, class) plus MEME columns;
    domain per column = majority domain over the sequences with a residue there (each
    sequence's own hmmscan envelope, through utils.domain_intervals);
  * recombination: every HyPhy test assumes one tree. --gard True runs a GARD screen; when
    it finds breakpoints FEL/MEME are re-run on the GARD-partitioned alignment and the
    branch-level tests carry a caveat. Without GARD the limitation is recorded in a
    'GARD' summary row;
  * NO SILENT FAILURE: a missing hyphy, a crashed method, an unparsable JSON or a malformed
    input raises (non-zero exit). Header-only/'not run' outputs only for genuinely empty
    or too-small inputs, or when run != True.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

SUMMARY_HEADER = "method\tstatistic\tvalue\tp_value\tn_sites_or_branches\tnote"
SITE_COLS = ["site", "alignment_col", "domain", "beta_alpha_or_dnds", "p_value", "classification",
             "fel_alpha", "fel_beta", "meme_p_value", "meme_classification", "n_seqs_with_codon",
             "domain_support"]
SITES_HEADER = "\t".join(SITE_COLS)
MIN_TIPS = 4
SITE_PVAL = 0.1
BRANCH_PVAL = 0.05
DEFAULT_SEED = 12345
ALL_METHODS = ("FEL", "MEME", "BUSTED", "aBSREL", "RELAX")
STATS = {"FEL": "sites_positive", "MEME": "sites_episodic", "BUSTED": "episodic_diversifying_p",
         "aBSREL": "branches_selected", "RELAX": "K"}
STOPS = {"TAA", "TAG", "TGA"}


class SelectionError(SystemExit):
    pass


def fail(msg):
    raise SelectionError(f"[run_selection] {msg}")


def truthy(s):
    return str(s).strip().lower() in ("true", "1", "yes", "t")


# =========================================================================
# inputs
# =========================================================================
def read_codon(path):
    recs = list(utils.read_fasta(path))
    if not recs:
        return recs
    names = [n for n, _ in recs]
    if len(set(names)) != len(names):
        fail(f"duplicated ids in {path}")
    lens = {len(s) for _, s in recs}
    if len(lens) != 1:
        fail(f"{path} is not an alignment (record lengths {sorted(lens)[:5]})")
    L = lens.pop()
    if L % 3:
        fail(f"{path}: alignment length {L} is not a multiple of 3")
    for n, s in recs:
        s = s.upper()
        for i in range(0, L, 3):
            if s[i:i + 3] in STOPS:
                fail(f"{path}: stop codon in {n} at codon {i // 3 + 1} (HyPhy rejects stop codons)")
    return [(n, s.upper()) for n, s in recs]


def tree_is_trivial(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return True
    with open(path) as fh:
        return fh.read().strip() in ("", "();", ";")


def load_tree(path):
    import dendropy
    return dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True,
                             rooting="force-unrooted")


def clean_newick(tree, tip_label=None, node_label=None):
    """Newick with tip labels, optional {annotations}, branch lengths; no support values."""
    def lab(nd):
        if nd.is_leaf():
            s = nd.taxon.label if tip_label is None else tip_label(nd)
        else:
            s = "" if node_label is None else node_label(nd)
        bl = nd.edge.length
        return s + (f":{bl:.10g}" if bl is not None and nd.parent_node is not None else "")

    def rec(nd):
        if nd.is_leaf():
            return lab(nd)
        return "(" + ",".join(rec(c) for c in nd.child_nodes()) + ")" + lab(nd)
    return rec(tree.seed_node) + ";"


def prune_tree(tree, keep):
    """Prune to `keep`; return an UNROOTED tree (basal bifurcation collapsed), so every
    branch -- and every RELAX label -- refers to one branch of the unrooted tree HyPhy fits."""
    keep = set(keep)
    tips = {lf.taxon.label for lf in tree.leaf_node_iter()}
    drop = tips - keep
    if drop:
        tree.prune_taxa_with_labels(drop, suppress_unifurcations=True)
    if len(tree.seed_node.child_nodes()) == 2 and len(keep) > 2:
        tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    return tree


def intact_loci(catalog_path):
    import pandas as pd
    df = pd.read_csv(catalog_path, sep="\t", dtype=str, keep_default_na=False)
    for c in ("locus_id", "catalytic_status"):
        if c not in df.columns:
            fail(f"catalog {catalog_path} lacks column '{c}' (SCHEMA prometheus_catalog.tsv)")
    return set(df.loc[df["catalytic_status"].str.strip().str.lower() == "intact", "locus_id"].str.strip())


# mcr_group values (scaffolds/methanogen_calls) whose host is neither a confident methanogen
# nor a confident non-methanogen: kept out of both RELAX sets (Unclassified nuisance set).
UNCERTAIN_MCR = {"anme", "ambiguous", "unplaced", "not_placed"}


def host_states(master_path):
    """genome_id -> 'M' (methanogen), 'N' (non-methanogen) or 'U' (unclassified).

    M: methanogen == True (mcrA-defined). U: ANME (methanotrophs with MCR), mcr hosts whose
    McrA placement is ambiguous/unplaced/not_placed, or an unparsable methanogen value.
    N: everything else with methanogen == False -- no mcrA, or mcr_group 'other' (MCR-like
    alkane oxidisers, not methanogens).
    """
    import pandas as pd
    df = pd.read_csv(master_path, sep="\t", dtype=str, keep_default_na=False)
    for c in ("genome_id", "methanogen"):
        if c not in df.columns:
            fail(f"master table {master_path} lacks column '{c}'")
    out = {}
    for _, r in df.iterrows():
        g = str(r["genome_id"]).strip()
        m = str(r["methanogen"]).strip().lower()
        grp = str(r.get("mcr_group", "")).strip().lower()
        anme = str(r.get("is_anme", "")).strip().lower() in ("true", "1", "yes", "t")
        if m in ("true", "1", "yes", "t"):
            out[g] = "M"
        elif anme or grp in UNCERTAIN_MCR or m not in ("false", "0", "no", "f"):
            out[g] = "U"
        else:
            out[g] = "N"
    return out


# =========================================================================
# HyPhy
# =========================================================================
def hyphy_bin():
    h = shutil.which("hyphy") or shutil.which("HYPHYMP")
    if h is None:
        fail("hyphy not found on PATH (selection.run is True)")
    return h


def run_hyphy(method, aln, tree, out_json, threads, extra=None, log_path=None, seed=None):
    cmd = [hyphy_bin()]
    if seed is not None:
        cmd.append(f"ENV=SetParameter(RANDOM_SEED,{int(seed)},0);")
    cmd += [f"CPU={max(1, int(threads))}", method, "--alignment", os.path.abspath(aln)]
    if tree:
        cmd += ["--tree", os.path.abspath(tree)]
    cmd += ["--output", os.path.abspath(out_json)] + list(extra or [])
    utils.log(f"[{method}] " + " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if log_path:
        with open(log_path, "w") as fh:
            fh.write(proc.stdout)
    if proc.returncode != 0:
        fail(f"hyphy {method} exited {proc.returncode}: {proc.stdout[-1500:]}")
    if not (os.path.exists(out_json) and os.path.getsize(out_json) > 0):
        fail(f"hyphy {method} wrote no JSON ({out_json}): {proc.stdout[-1500:]}")
    try:
        with open(out_json) as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        fail(f"hyphy {method} JSON {out_json} is not valid JSON: {e}")


def _norm(h):
    """'&beta;<sup>+</sup>' -> 'beta+', 'p-value' -> 'p-value', 'alpha=beta' kept."""
    h = re.sub(r"</?sup>", "", str(h))
    h = re.sub(r"&([a-zA-Z]+);", r"\1", h)
    return h.strip().lower()


def mle_table(js, method):
    """[(alignment_site_0based, {normalised_header: value})] from an MLE block."""
    mle = js.get("MLE")
    if not isinstance(mle, dict) or "headers" not in mle or "content" not in mle:
        fail(f"{method} JSON has no MLE headers/content block")
    hdr = [_norm(h[0] if isinstance(h, (list, tuple)) else h) for h in mle["headers"]]
    parts = js.get("data partitions", {}) or {}
    out = []
    for key in sorted(mle["content"], key=lambda k: int(k)):
        rows = mle["content"][key]
        cov = parts.get(str(key), {}).get("coverage")
        sites = cov[0] if (isinstance(cov, list) and cov and isinstance(cov[0], list)) else None
        if sites is not None and len(sites) != len(rows):
            fail(f"{method}: partition {key} has {len(rows)} MLE rows but {len(sites)} covered sites")
        if sites is None:
            if len(mle["content"]) > 1:
                fail(f"{method}: multi-partition MLE without 'data partitions' coverage")
            sites = list(range(len(rows)))
        for s, r in zip(sites, rows):
            if len(r) != len(hdr):
                fail(f"{method}: MLE row width {len(r)} != {len(hdr)} headers")
            out.append((int(s), dict(zip(hdr, r))))
    return out


def fnum(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return x


def need(d, key, method):
    if key not in d:
        fail(f"{method}: expected column '{key}' not in MLE headers {sorted(d)}")
    return fnum(d[key])


def parse_fel(js, site_p):
    rows = {}
    for s, d in mle_table(js, "FEL"):
        a, b, p = need(d, "alpha", "FEL"), need(d, "beta", "FEL"), need(d, "p-value", "FEL")
        if p <= site_p and b > a:
            cls = "positive"
        elif p <= site_p and b < a:
            cls = "negative"
        else:
            cls = "neutral"
        if a > 0:
            ratio = b / a
        elif b > 0:
            ratio = float("inf")     # beta > 0 with alpha = 0: dN/dS unbounded
        else:
            ratio = float("nan")     # invariant / no information
        rows[s] = {"alpha": a, "beta": b, "p": p, "ratio": ratio, "class": cls}
    return rows


def parse_meme(js, site_p):
    rows = {}
    for s, d in mle_table(js, "MEME"):
        p = need(d, "p-value", "MEME")
        rows[s] = {"p": p, "class": "episodic" if p <= site_p else "neutral",
                   "beta_plus": fnum(d.get("beta+")), "alpha": fnum(d.get("alpha"))}
    return rows


def test_results(js, method, keys):
    tr = js.get("test results")
    if not isinstance(tr, dict):
        fail(f"{method} JSON has no 'test results' block")
    out = {}
    for k in keys:
        if k not in tr:
            fail(f"{method} 'test results' lacks '{k}' (has {sorted(tr)})")
        out[k] = fnum(tr[k])
    return out


def n_branches(js):
    try:
        return sum(1 for _ in js["branch attributes"]["0"])
    except (KeyError, TypeError):
        return ""


# =========================================================================
# RELAX labelling
# =========================================================================
def relax_labels(tree, states, mode):
    """Assign 'Test' / 'Reference' to branches (keyed by the node BELOW the branch).

    states: tip label -> 'M' | 'N' | 'U'. The tree is unrooted (basal bifurcation collapsed),
    so a branch is a bipartition of the tips and the rule must not depend on where the
    arbitrary newick root sits:
      'clade'    Test if one side of the branch is all-M (the branch lies inside, or leads
                 to, a methanogen-only group), Reference if one side is all-N; a branch with
                 an all-M side AND an all-N side (the single transition branch) and every
                 other branch stay unlabelled (RELAX 'Unclassified' nuisance set);
      'terminal' tip branches only, by the host state of the tip.
    """
    all_states = [states.get(lf.taxon.label, "U") for lf in tree.leaf_node_iter()]
    total = {}
    for st in all_states:
        total[st] = total.get(st, 0) + 1
    desc = {}
    lab = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            cnt = {states.get(nd.taxon.label, "U"): 1}
        else:
            cnt = {}
            for ch in nd.child_nodes():
                for k, v in desc[ch].items():
                    cnt[k] = cnt.get(k, 0) + v
        desc[nd] = cnt
        if nd.parent_node is None or (mode == "terminal" and not nd.is_leaf()):
            continue
        other = {k: total[k] - cnt.get(k, 0) for k in total if total[k] - cnt.get(k, 0) > 0}
        sides = (set(cnt), set(other))
        is_m = {"M"} in sides
        is_n = {"N"} in sides
        if is_m and not is_n:
            lab[nd] = "Test"
        elif is_n and not is_m:
            lab[nd] = "Reference"
    return lab


# =========================================================================
# site -> domain
# =========================================================================
def column_domains(records, hits, spec):
    """codon column (0-based) -> (majority domain or 'NA', support fraction, n residues)."""
    from collections import Counter
    L = len(records[0][1]) // 3 if records else 0
    votes = [Counter() for _ in range(L)]
    cover = [0] * L
    for name, seq in records:
        env = [dict(h, start=h.get("env_from", h["start"]), end=h.get("env_to", h["end"]))
               for h in hits.get(name, [])]
        ivs = utils.domain_intervals(env, spec)
        res_dom = {}
        for cat, spans in ivs.items():
            for a, b in spans:
                for r in range(a, b + 1):
                    res_dom.setdefault(r, cat)
        aa = 0
        for c in range(L):
            cod = seq[3 * c:3 * c + 3]
            if set(cod) <= set("-."):
                continue          # gap codon (an NNN codon is still a residue)
            aa += 1
            cover[c] += 1
            votes[c][res_dom.get(aa, "NA")] += 1
    out = []
    for c in range(L):
        if not cover[c]:
            out.append(("NA", "", 0))
            continue
        dom, n = votes[c].most_common(1)[0]
        frac = n / cover[c]
        out.append((dom if (dom != "NA" and frac >= 0.5) else "NA", round(frac, 4), cover[c]))
    return out


# =========================================================================
# GARD
# =========================================================================
def run_gard(codon_path, workdir, threads, mode, seed=None):
    out_json = os.path.join(workdir, "GARD.json")
    best = os.path.join(workdir, "GARD.best-gard")
    js = run_hyphy("gard", codon_path, None, out_json, threads,
                   extra=["--type", "codon", "--mode", mode, "--output-lf", best],
                   log_path=os.path.join(workdir, "GARD.log"), seed=seed)
    bp = js.get("breakpointData")
    if not isinstance(bp, dict) or not bp:
        fail("GARD JSON has no breakpointData")
    n_bp = len(bp) - 1
    single, bestaic = fnum(js.get("singleTreeAICc")), fnum(js.get("bestModelAICc"))
    return n_bp, single - bestaic, (best if n_bp > 0 and os.path.exists(best) else None)


# =========================================================================
# main
# =========================================================================
def write_outputs(summary_rows, site_rows, args):
    with open(args.out_summary, "w") as out:
        out.write(SUMMARY_HEADER + "\n")
        for r in summary_rows:
            out.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    with open(args.out_sites, "w") as out:
        out.write(SITES_HEADER + "\n")
        for r in site_rows:
            out.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in SITE_COLS) + "\n")


def fmt(x, nd=6):
    if x is None or (isinstance(x, float) and (math.isnan(x))):
        return ""
    if isinstance(x, float) and math.isinf(x):
        return "inf"
    return round(x, nd) if isinstance(x, float) else x


def not_run(methods, reason, gard_row=True):
    rows = [(m, STATS[m], "", "", "", f"not run: {reason}") for m in methods]
    if gard_row:
        rows.append(("GARD", "breakpoints", "", "", "", f"not run: {reason}"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--codon", required=True)
    ap.add_argument("--gene-tree", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--methods", required=True, help="comma list: FEL,MEME,BUSTED,aBSREL,RELAX")
    ap.add_argument("--relax-partition", default="methanogen", help="methanogen | none")
    ap.add_argument("--relax-labelling", default="clade", choices=("clade", "terminal"))
    ap.add_argument("--intact-only", default="True")
    ap.add_argument("--site-pvalue", type=float, default=SITE_PVAL)
    ap.add_argument("--branch-pvalue", type=float, default=BRANCH_PVAL,
                    help="aBSREL corrected-p threshold (HyPhy default 0.05)")
    ap.add_argument("--gard", default="False", help="True: GARD recombination screen first")
    ap.add_argument("--gard-mode", default="Faster", choices=("Faster", "Normal"))
    ap.add_argument("--run", default="True")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED,
                    help="HyPhy RNG seed (SetParameter(RANDOM_SEED)); default 12345")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-sites", required=True)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    for p in (args.out_summary, args.out_sites):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    bad = [m for m in methods if m not in ALL_METHODS]
    if bad:
        fail(f"unknown method(s) {bad}; allowed {ALL_METHODS}")

    if not truthy(args.run):
        write_outputs(not_run(methods, "evolution.selection.run disabled"), [], args)
        utils.log("selection disabled; 'not run' rows written")
        return
    records = read_codon(args.codon)
    if len(records) < MIN_TIPS or tree_is_trivial(args.gene_tree):
        write_outputs(not_run(methods, f"too few codon sequences ({len(records)}) or trivial gene tree"),
                      [], args)
        utils.log("too few sequences / trivial tree; 'not run' rows written")
        return

    spec = utils.load_domain_spec(args.domains)
    hits = utils.parse_domtbl(args.domtbl)
    tree = load_tree(args.gene_tree)
    tips = {lf.taxon.label for lf in tree.leaf_node_iter()}
    keep = [n for n, _ in records if n in tips]
    notes = []
    extra_tips = tips - {n for n, _ in records}
    if extra_tips:                       # e.g. reference/paralog tips, or loci without a CDS
        notes.append(f"{len(extra_tips)} gene-tree tip(s) without a codon record pruned")
    if len(keep) < len(records):
        notes.append(f"{len(records) - len(keep)} codon record(s) absent from the gene tree dropped")
    if truthy(args.intact_only):
        intact = intact_loci(args.catalog)
        n0 = len(keep)
        keep = [n for n in keep if n in intact]
        notes.append(f"intact-only: {len(keep)}/{n0} loci catalytically intact")
    if len(keep) < MIN_TIPS:
        write_outputs(not_run(methods, f"only {len(keep)} analysable loci (<{MIN_TIPS}); "
                              + "; ".join(notes)), [], args)
        utils.log(f"only {len(keep)} analysable loci; 'not run' rows written")
        return
    keep_set = set(keep)
    records = [(n, s) for n, s in records if n in keep_set]
    codon_path = os.path.join(args.workdir, "codon.analysed.fasta")
    utils.write_fasta(records, codon_path)
    tree = prune_tree(tree, keep_set)
    tree_path = os.path.join(args.workdir, "gene.analysed.nwk")
    with open(tree_path, "w") as fh:
        fh.write(clean_newick(tree) + "\n")
    base_note = f"n={len(records)} loci x {len(records[0][1]) // 3} codons; " + "; ".join(notes)
    site_p = float(args.site_pvalue)

    # ---- recombination screen -------------------------------------------------
    summary, site_aln = [], codon_path
    rec_caveat = "assumes no recombination (single tree; GARD not run)"
    if truthy(args.gard):
        n_bp, d_aic, best = run_gard(codon_path, args.workdir, args.threads, args.gard_mode, seed=args.seed)
        summary.append(("GARD", "breakpoints", n_bp, "", len(records[0][1]) // 3,
                        f"GARD {args.gard_mode} codon screen; dAICc(single tree - best)={fmt(d_aic, 3)}"))
        if n_bp > 0:
            site_aln = best
            rec_caveat = (f"GARD found {n_bp} breakpoint(s): FEL/MEME use the GARD partitions; "
                          "this branch-level test uses one tree and may be biased")
        else:
            rec_caveat = "GARD found no breakpoint"
    else:
        summary.append(("GARD", "breakpoints", "", "", "",
                        "not run (evolution.selection.gard false): all HyPhy tests assume a single "
                        "tree, which domain shuffling/recombination would violate"))

    fel, meme = {}, {}
    site_tree = None if site_aln != codon_path else tree_path
    if "FEL" in methods:
        js = run_hyphy("fel", site_aln, site_tree, os.path.join(args.workdir, "FEL.json"), args.threads,
                       extra=["--pvalue", str(site_p)], log_path=os.path.join(args.workdir, "FEL.log"),
                       seed=args.seed)
        fel = parse_fel(js, site_p)
        n_pos = sum(1 for r in fel.values() if r["class"] == "positive")
        n_neg = sum(1 for r in fel.values() if r["class"] == "negative")
        summary.append(("FEL", "sites_positive", n_pos, "", len(fel),
                        f"p<={site_p} and beta>alpha; {base_note}"))
        summary.append(("FEL", "sites_negative", n_neg, "", len(fel), f"p<={site_p} and beta<alpha"))
    if "MEME" in methods:
        js = run_hyphy("meme", site_aln, site_tree, os.path.join(args.workdir, "MEME.json"), args.threads,
                       extra=["--pvalue", str(site_p)], log_path=os.path.join(args.workdir, "MEME.log"),
                       seed=args.seed)
        meme = parse_meme(js, site_p)
        summary.append(("MEME", "sites_episodic", sum(1 for r in meme.values() if r["class"] == "episodic"),
                        "", len(meme), f"p<={site_p}; {base_note}"))
    if "BUSTED" in methods:
        js = run_hyphy("busted", codon_path, tree_path, os.path.join(args.workdir, "BUSTED.json"),
                       args.threads, log_path=os.path.join(args.workdir, "BUSTED.log"), seed=args.seed)
        tr = test_results(js, "BUSTED", ("LRT", "p-value"))
        summary.append(("BUSTED", "episodic_diversifying_p", fmt(tr["p-value"]), fmt(tr["p-value"]),
                        n_branches(js), f"gene-wide LRT={fmt(tr['LRT'], 4)}; {rec_caveat}"))
    if "aBSREL" in methods:
        js = run_hyphy("absrel", codon_path, tree_path, os.path.join(args.workdir, "aBSREL.json"),
                       args.threads, log_path=os.path.join(args.workdir, "aBSREL.log"), seed=args.seed)
        tr = test_results(js, "aBSREL", ("positive test results", "tested"))
        n_sel = tr["positive test results"]
        if abs(args.branch_pvalue - fnum(js.get("test results", {}).get("P-value threshold", BRANCH_PVAL))) > 1e-12:
            n_sel = 0
            for _b, d in (js.get("branch attributes", {}).get("0", {}) or {}).items():
                if isinstance(d, dict) and fnum(d.get("Corrected P-value")) <= args.branch_pvalue:
                    n_sel += 1
        summary.append(("aBSREL", "branches_selected", int(n_sel), "", int(tr["tested"]),
                        f"Holm-corrected p<={args.branch_pvalue}; {rec_caveat}"))
    if "RELAX" in methods:
        if args.relax_partition != "methanogen":
            summary.append(("RELAX", "K", "", "", "", f"not run: relax-partition '{args.relax_partition}'"))
        else:
            states = host_states(args.master)
            tip_state = {n: states.get(utils.genome_of(n), "U") for n in keep}
            labels = relax_labels(tree, tip_state, args.relax_labelling)
            n_test = sum(1 for v in labels.values() if v == "Test")
            n_ref = sum(1 for v in labels.values() if v == "Reference")
            n_tip_t = sum(1 for nd, v in labels.items() if v == "Test" and nd.is_leaf())
            if n_test < 1 or n_ref < 1:
                summary.append(("RELAX", "K", "", "", "", f"not run: no Test/Reference split "
                                f"(Test={n_test}, Reference={n_ref} branches)"))
            else:
                lab_path = os.path.join(args.workdir, "relax.labelled.nwk")
                with open(lab_path, "w") as fh:
                    fh.write(clean_newick(
                        tree,
                        tip_label=lambda nd: nd.taxon.label + (f"{{{labels[nd]}}}" if nd in labels else ""),
                        node_label=lambda nd: f"{{{labels[nd]}}}" if nd in labels else "") + "\n")
                js = run_hyphy("relax", codon_path, lab_path, os.path.join(args.workdir, "RELAX.json"),
                               args.threads, extra=["--test", "Test", "--reference", "Reference"],
                               log_path=os.path.join(args.workdir, "RELAX.log"), seed=args.seed)
                tr = test_results(js, "RELAX", ("LRT", "p-value", "relaxation or intensification parameter"))
                n_unc = sum(1 for nd in tree.preorder_node_iter() if nd.parent_node is not None) - n_test - n_ref
                summary.append(("RELAX", "K", fmt(tr["relaxation or intensification parameter"], 4),
                                fmt(tr["p-value"]), n_test + n_ref,
                                f"{args.relax_labelling} labelling: Test={n_test} methanogen branches "
                                f"({n_tip_t} terminal), Reference={n_ref}, Unclassified={n_unc} (nuisance); "
                                f"LRT={fmt(tr['LRT'], 4)}; K>1 intensified, K<1 relaxed; {rec_caveat}"))

    # ---- site table -------------------------------------------------------------
    site_rows = []
    if fel or meme:
        doms = column_domains(records, hits, spec)
        L = len(records[0][1]) // 3
        for c in range(L):
            f, m = fel.get(c), meme.get(c)
            if f is None and m is None:
                continue
            dom, sup, cov = doms[c]
            src = f if f is not None else None
            site_rows.append({
                "site": c + 1, "alignment_col": c + 1, "domain": dom,
                "beta_alpha_or_dnds": fmt(src["ratio"], 4) if src else "",
                "p_value": fmt(src["p"]) if src else (fmt(m["p"]) if m else ""),
                "classification": src["class"] if src else ("positive" if m and m["class"] == "episodic" else "neutral"),
                "fel_alpha": fmt(f["alpha"]) if f else "", "fel_beta": fmt(f["beta"]) if f else "",
                "meme_p_value": fmt(m["p"]) if m else "", "meme_classification": m["class"] if m else "",
                "n_seqs_with_codon": cov, "domain_support": sup,
            })
    write_outputs(summary, site_rows, args)
    utils.log(f"selection: {len(summary)} summary row(s), {len(site_rows)} site row(s)")


if __name__ == "__main__":
    main()
