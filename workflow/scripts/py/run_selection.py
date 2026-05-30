"""Selection analysis on the Prometheus codon alignment (HyPhy FEL/MEME/BUSTED/aBSREL/RELAX).

Estimates the selection regime on catalytically-intact Prometheus loci (Stage 5.6).
Diversifying selection concentrated on the nuclease/REase domains would be the
molecular signature of host-virus conflict; RELAX asks whether that pressure differs
on methanogen branches.

Each HyPhy method is run independently and guarded, so one failure never blocks the
others or the DAG. RELAX additionally writes a HyPhy-labelled newick that tags
methanogen tips ``{Test}`` and the rest ``{Reference}`` (methanogen genome set read
from the master table).

Outputs:
  selection_summary.tsv  one headline row per method (statistic, value, p, n, note)
  site_table.tsv         per-site FEL/MEME results mapped to a domain

Guards: run != 'True', <4 codon sequences, or a trivial gene tree -> header-only
summary + sites. intact-only restricts to catalog catalytic_status=='intact'; if that
leaves <4 tips, header-only with a note.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

SUMMARY_HEADER = "method\tstatistic\tvalue\tp_value\tn_sites_or_branches\tnote"
SITES_HEADER = "site\talignment_col\tdomain\tbeta_alpha_or_dnds\tp_value\tclassification"
MIN_TIPS = 4
SITE_PVAL = 0.1  # FEL/MEME positive-selection significance threshold


# =========================================================================
# tree / alignment helpers
# =========================================================================
def tree_text(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return ""


def is_trivial(path):
    return tree_text(path) in ("", "();", ";")


def load_tree(path):
    import dendropy
    return dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)


def prune_tree_to(tree, keep_labels, out_path):
    """Retain only `keep_labels` tips; write newick. Returns n retained tips."""
    keep = set(keep_labels)
    try:
        tree.retain_taxa_with_labels(keep)
    except Exception as e:  # noqa: BLE001
        utils.log(f"tree prune failed: {e}")
        return 0
    n = sum(1 for _ in tree.leaf_node_iter())
    tree.write(path=out_path, schema="newick", unquoted_underscores=True,
               suppress_rooting=True)
    return n


# =========================================================================
# catalog / master-table inputs
# =========================================================================
def intact_loci(catalog_path):
    """Set of locus_ids with catalytic_status == 'intact'."""
    intact = set()
    if not catalog_path or not os.path.exists(catalog_path):
        return intact
    try:
        import pandas as pd
        df = pd.read_csv(catalog_path, sep="\t", dtype=str)
        if "catalytic_status" in df.columns and "locus_id" in df.columns:
            sub = df[df["catalytic_status"].astype(str).str.lower() == "intact"]
            intact = set(sub["locus_id"].astype(str))
    except Exception as e:  # noqa: BLE001
        utils.log(f"could not read catalog: {e}")
    return intact


def methanogen_genomes(master_path):
    """Set of genome_ids flagged methanogen==True in the master table."""
    mg = set()
    if not master_path or not os.path.exists(master_path):
        return mg
    try:
        import pandas as pd
        df = pd.read_csv(master_path, sep="\t", dtype=str)
        if "methanogen" in df.columns and "genome_id" in df.columns:
            sub = df[df["methanogen"].astype(str).str.lower().isin(("true", "1", "yes"))]
            mg = set(sub["genome_id"].astype(str))
    except Exception as e:  # noqa: BLE001
        utils.log(f"could not read master table: {e}")
    return mg


# =========================================================================
# HyPhy runners
# =========================================================================
def hyphy_path():
    return shutil.which("hyphy") or shutil.which("HYPHYMP")


def run_hyphy(method, codon, tree, out_json, extra=None):
    """Run one HyPhy analysis; return True on success (output JSON present)."""
    hyphy = hyphy_path()
    if hyphy is None:
        utils.log(f"[{method}] hyphy not on PATH; skipping")
        return False
    cmd = [hyphy, method, "--alignment", os.path.abspath(codon),
           "--tree", os.path.abspath(tree), "--output", os.path.abspath(out_json)]
    if extra:
        cmd += extra
    utils.log(f"[{method}] running: " + " ".join(cmd))
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if proc.returncode != 0:
            err = proc.stderr.decode(errors="replace").strip()
            utils.log(f"[{method}] exit {proc.returncode}: {err[:500]}")
            return False
    except Exception as e:  # noqa: BLE001
        utils.log(f"[{method}] failed: {e}")
        return False
    if not (os.path.exists(out_json) and os.path.getsize(out_json) > 0):
        utils.log(f"[{method}] produced no JSON output")
        return False
    return True


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except Exception as e:  # noqa: BLE001
        utils.log(f"could not parse JSON {path}: {e}")
        return None


# ---- per-method JSON parsing ------------------------------------------------
def _mle_rows(js):
    """Return (headers, rows) from a HyPhy MLE content block (FEL/MEME)."""
    content = js.get("MLE", {}).get("content", {})
    headers = js.get("MLE", {}).get("headers", [])
    # content is keyed by partition index ('0', '1', ...); concatenate partitions
    rows = []
    for _k, part in sorted(content.items()):
        if isinstance(part, list):
            rows.extend(part)
    # headers are [[name, description], ...]
    hdr_names = [h[0] if isinstance(h, (list, tuple)) else str(h) for h in headers]
    return hdr_names, rows


def _col_idx(hdr_names, *candidates):
    low = [h.lower() for h in hdr_names]
    for cand in candidates:
        c = cand.lower()
        for i, h in enumerate(low):
            if h == c:
                return i
        for i, h in enumerate(low):
            if c in h:
                return i
    return None


def parse_fel(path):
    """FEL: per-site alpha/beta + p-value. Return (n_pos_sites, site_records)."""
    js = load_json(path)
    if not js:
        return None, []
    hdr, rows = _mle_rows(js)
    i_alpha = _col_idx(hdr, "alpha")
    i_beta = _col_idx(hdr, "beta")
    i_p = _col_idx(hdr, "p-value", "p_value", "pvalue")
    sites = []
    n_pos = 0
    for idx, r in enumerate(rows, start=1):
        try:
            alpha = float(r[i_alpha]) if i_alpha is not None else float("nan")
            beta = float(r[i_beta]) if i_beta is not None else float("nan")
            pval = float(r[i_p]) if i_p is not None else float("nan")
        except (TypeError, ValueError, IndexError):
            continue
        # beta>alpha => positive (diversifying); beta<alpha => negative (purifying)
        if pval <= SITE_PVAL and beta > alpha:
            cls = "positive"
        elif pval <= SITE_PVAL and beta < alpha:
            cls = "negative"
        else:
            cls = "neutral"
        if cls == "positive":
            n_pos += 1
        ratio = (beta / alpha) if (alpha and alpha > 0) else ""
        sites.append({"site": idx, "stat": (round(ratio, 4) if ratio != "" else ""),
                      "p": pval, "class": cls})
    return n_pos, sites


def parse_meme(path):
    """MEME: per-site episodic positive selection p-value. Return (n_pos, sites)."""
    js = load_json(path)
    if not js:
        return None, []
    hdr, rows = _mle_rows(js)
    i_p = _col_idx(hdr, "p-value", "p_value", "pvalue")
    i_beta_plus = _col_idx(hdr, "beta+", "&beta;+", "beta (positive)")
    i_alpha = _col_idx(hdr, "alpha")
    sites = []
    n_pos = 0
    for idx, r in enumerate(rows, start=1):
        try:
            pval = float(r[i_p]) if i_p is not None else float("nan")
        except (TypeError, ValueError, IndexError):
            continue
        cls = "positive" if pval <= SITE_PVAL else "neutral"
        if cls == "positive":
            n_pos += 1
        stat = ""
        try:
            if i_beta_plus is not None and i_alpha is not None:
                bp = float(r[i_beta_plus]); a = float(r[i_alpha])
                stat = round(bp / a, 4) if a > 0 else ""
        except (TypeError, ValueError, IndexError):
            stat = ""
        sites.append({"site": idx, "stat": stat, "p": pval, "class": cls})
    return n_pos, sites


def parse_busted(path):
    """BUSTED: gene-wide test p-value for episodic diversifying selection."""
    js = load_json(path)
    if not js:
        return None
    tr = js.get("test results", {}) or js.get("test_results", {})
    p = tr.get("p-value", tr.get("p_value"))
    try:
        return float(p)
    except (TypeError, ValueError):
        return None


def parse_absrel(path):
    """aBSREL: number of branches under selection (corrected p <= 0.05)."""
    js = load_json(path)
    if not js:
        return None
    n = js.get("test results", {}).get("positive test results")
    if isinstance(n, (int, float)):
        return int(n)
    # fall back: count branch attributes with corrected p <= 0.05
    count = 0
    found = False
    battr = js.get("branch attributes", {}).get("0", {})
    for _b, d in battr.items():
        if not isinstance(d, dict):
            continue
        p = d.get("Corrected P-value", d.get("Corrected P-value "))
        if p is None:
            continue
        found = True
        try:
            if float(p) <= 0.05:
                count += 1
        except (TypeError, ValueError):
            pass
    return count if found else None


def parse_relax(path):
    """RELAX: selection-intensity K and the LRT p-value."""
    js = load_json(path)
    if not js:
        return None, None
    tr = js.get("test results", {})
    p = tr.get("p-value", tr.get("p_value"))
    k = tr.get("relaxation or intensification parameter", tr.get("k"))
    try:
        p = float(p)
    except (TypeError, ValueError):
        p = None
    try:
        k = float(k)
    except (TypeError, ValueError):
        k = None
    return k, p


# =========================================================================
# RELAX branch labelling
# =========================================================================
def write_relax_labelled(gene_tree_path, methanogens, out_path):
    """Write a newick where methanogen tips are name{Test}, others name{Reference}.

    Returns (n_test, n_reference). Uses HyPhy node-label syntax (label in braces
    immediately after the taxon name). Falls back to no labelling if dendropy fails.
    """
    try:
        t = load_tree(gene_tree_path)
    except Exception as e:  # noqa: BLE001
        utils.log(f"RELAX: could not load gene tree: {e}")
        return 0, 0
    n_test = n_ref = 0
    for leaf in t.leaf_node_iter():
        if leaf.taxon is None:
            continue
        label = leaf.taxon.label
        genome = utils.genome_of(label)
        tag = "Test" if genome in methanogens else "Reference"
        if tag == "Test":
            n_test += 1
        else:
            n_ref += 1
        leaf.taxon.label = f"{label}{{{tag}}}"
    # write with underscores preserved; braces must survive, so don't quote
    newick = t.as_string(schema="newick", unquoted_underscores=True,
                         suppress_rooting=True)
    # dendropy may quote labels containing '{'; strip the surrounding quotes.
    newick = newick.replace("'", "")
    with open(out_path, "w") as out:
        out.write(newick.strip() + "\n")
    return n_test, n_ref


# =========================================================================
# site -> domain mapping
# =========================================================================
def _matches_domain(hit, names, accs):
    dom = str(hit.get("domain", "")).lower()
    acc = str(hit.get("accession", ""))
    acc_base = acc.split(".")[0]
    if accs and (acc in accs or acc_base in accs):
        return True
    for nm in names:
        if nm and nm.lower() in dom:
            return True
    return False


def build_codon_domain_map(codon_records, domains, hits):
    """Approximate map: codon alignment column -> domain.

    Approximation (documented): take the reference (longest un-gapped) codon-aligned
    sequence; walk its alignment columns; for each ungapped codon, advance a residue
    counter and look up which domain (if any) covers that residue in the reference
    protein's domtbl ranges. Project that residue->domain assignment onto the
    alignment column. Columns where the reference is gapped, or no domain covers the
    residue, are 'NA'. This uses a single reference sequence as the coordinate frame
    rather than a per-column majority vote, which is a coarse but transparent proxy.

    Returns dict {alignment_col (1-based) -> domain_name}.
    """
    col_domain = {}
    if not codon_records:
        return col_domain
    # choose reference = record with the most non-gap characters
    ref_id, ref_seq = max(codon_records, key=lambda kv: sum(c not in "-." for c in kv[1]))
    ref_hits = hits.get(ref_id, [])
    if not ref_hits:
        utils.log(f"site->domain: reference {ref_id} has no domtbl hits; all NA")
        return col_domain

    # build residue(1-based aa) -> domain for the reference protein
    res_domain = {}
    for dname, spec in domains.items():
        spec = spec or {}
        names = spec.get("names", [])
        accs = set(spec.get("accessions", []))
        for h in ref_hits:
            if _matches_domain(h, names, accs):
                for r in range(int(h["start"]), int(h["end"]) + 1):
                    res_domain.setdefault(r, dname)

    if not res_domain:
        return col_domain

    # walk codon-alignment columns; each codon = 3 nt; ungapped codons advance aa idx
    aa_idx = 0
    n = len(ref_seq)
    col = 0
    i = 0
    while i < n:
        codon = ref_seq[i:i + 3]
        col += 1
        is_gap = all(c in "-." for c in codon)
        if not is_gap:
            aa_idx += 1
            dom = res_domain.get(aa_idx)
            if dom:
                col_domain[col] = dom
        i += 3
    return col_domain


# =========================================================================
# main
# =========================================================================
def write_headers(summary_path, sites_path):
    utils.write_empty(summary_path, SUMMARY_HEADER)
    utils.write_empty(sites_path, SITES_HEADER)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--codon", required=True)
    ap.add_argument("--gene-tree", required=True)
    ap.add_argument("--catalog", required=True)
    ap.add_argument("--master", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--domains", required=True)
    ap.add_argument("--methods", required=True, help="comma list: FEL,MEME,BUSTED,aBSREL,RELAX")
    ap.add_argument("--relax-partition", default="methanogen")
    ap.add_argument("--intact-only", default="True")
    ap.add_argument("--run", default="True")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-sites", required=True)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    for p in (args.out_summary, args.out_sites):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    # ---- top-level guards ----
    if str(args.run) != "True":
        write_headers(args.out_summary, args.out_sites)
        utils.log("run != True; header-only outputs")
        return

    codon_records = list(utils.read_fasta(args.codon))
    if len(codon_records) < MIN_TIPS or is_trivial(args.gene_tree):
        write_headers(args.out_summary, args.out_sites)
        utils.log(f"too few codon seqs ({len(codon_records)}) or trivial gene tree; "
                  "header-only outputs")
        return

    try:
        domains = json.loads(args.domains) if args.domains else {}
    except json.JSONDecodeError:
        domains = {}
    hits = utils.parse_domtbl(args.domtbl)

    # ---- intact-only restriction ----
    codon_path = os.path.abspath(args.codon)
    tree_path = os.path.abspath(args.gene_tree)
    if str(args.intact_only) == "True":
        intact = intact_loci(args.catalog)
        if intact:
            kept = [(n, s) for n, s in codon_records if n in intact]
            if len(kept) < MIN_TIPS:
                write_headers(args.out_summary, args.out_sites)
                utils.log(f"intact-only left {len(kept)} loci (<{MIN_TIPS}); "
                          "header-only outputs")
                return
            codon_path = os.path.join(args.workdir, "codon.intact.fasta")
            utils.write_fasta(kept, codon_path)
            codon_records = kept
            # prune the gene tree to the intact tips
            try:
                tr = load_tree(args.gene_tree)
                tree_path = os.path.join(args.workdir, "gene.intact.nwk")
                n_kept = prune_tree_to(tr, [n for n, _ in kept], tree_path)
                if n_kept < MIN_TIPS:
                    write_headers(args.out_summary, args.out_sites)
                    utils.log(f"intact-only tree has {n_kept} tips (<{MIN_TIPS}); "
                              "header-only outputs")
                    return
            except Exception as e:  # noqa: BLE001
                utils.log(f"could not prune gene tree to intact tips: {e}; "
                          "using full gene tree")
                tree_path = os.path.abspath(args.gene_tree)
        else:
            utils.log("intact-only requested but no intact loci in catalog; "
                      "using all loci")

    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    summary_rows = []
    fel_sites, meme_sites = [], []

    # ---- FEL ----
    if "FEL" in methods:
        out = os.path.join(args.workdir, "FEL.json")
        if run_hyphy("fel", codon_path, tree_path, out):
            n_pos, sites = parse_fel(out)
            fel_sites = sites
            if n_pos is None:
                summary_rows.append(("FEL", "sites_positive", "", "", "", "ran; unparsed JSON"))
            else:
                summary_rows.append(("FEL", "sites_positive", n_pos, "", n_pos,
                                     f"p<={SITE_PVAL}, beta>alpha"))
        else:
            summary_rows.append(("FEL", "sites_positive", "", "", "", "skipped/failed"))

    # ---- MEME ----
    if "MEME" in methods:
        out = os.path.join(args.workdir, "MEME.json")
        if run_hyphy("meme", codon_path, tree_path, out):
            n_pos, sites = parse_meme(out)
            meme_sites = sites
            if n_pos is None:
                summary_rows.append(("MEME", "sites_episodic", "", "", "", "ran; unparsed JSON"))
            else:
                summary_rows.append(("MEME", "sites_episodic", n_pos, "", n_pos,
                                     f"p<={SITE_PVAL}"))
        else:
            summary_rows.append(("MEME", "sites_episodic", "", "", "", "skipped/failed"))

    # ---- BUSTED ----
    if "BUSTED" in methods:
        out = os.path.join(args.workdir, "BUSTED.json")
        if run_hyphy("busted", codon_path, tree_path, out):
            p = parse_busted(out)
            summary_rows.append(("BUSTED", "episodic_diversifying_p",
                                 (round(p, 6) if p is not None else ""),
                                 (round(p, 6) if p is not None else ""), "",
                                 "gene-wide test"))
        else:
            summary_rows.append(("BUSTED", "episodic_diversifying_p", "", "", "",
                                 "skipped/failed"))

    # ---- aBSREL ----
    if "aBSREL" in methods:
        out = os.path.join(args.workdir, "aBSREL.json")
        if run_hyphy("absrel", codon_path, tree_path, out):
            n = parse_absrel(out)
            summary_rows.append(("aBSREL", "branches_selected",
                                 (n if n is not None else ""), "",
                                 (n if n is not None else ""),
                                 "corrected p<=0.05"))
        else:
            summary_rows.append(("aBSREL", "branches_selected", "", "", "",
                                 "skipped/failed"))

    # ---- RELAX (methanogen Test vs Reference) ----
    if "RELAX" in methods:
        if args.relax_partition == "methanogen":
            mg = methanogen_genomes(args.master)
            labelled = os.path.join(args.workdir, "relax.labelled.nwk")
            n_test, n_ref = write_relax_labelled(tree_path, mg, labelled)
            if n_test < 1 or n_ref < 1:
                summary_rows.append(("RELAX", "K", "", "", "",
                                     f"no valid Test/Reference split "
                                     f"(Test={n_test}, Reference={n_ref})"))
            else:
                out = os.path.join(args.workdir, "RELAX.json")
                if run_hyphy("relax", codon_path, labelled, out,
                             extra=["--test", "Test", "--reference", "Reference"]):
                    k, p = parse_relax(out)
                    summary_rows.append(("RELAX", "K",
                                         (round(k, 4) if k is not None else ""),
                                         (round(p, 6) if p is not None else ""), "",
                                         f"Test={n_test} methanogen, Reference={n_ref}"))
                else:
                    summary_rows.append(("RELAX", "K", "", "", "", "skipped/failed"))
        else:
            summary_rows.append(("RELAX", "K", "", "", "",
                                 f"unsupported relax-partition '{args.relax_partition}'"))

    # ---- write summary ----
    with open(args.out_summary, "w") as out:
        out.write(SUMMARY_HEADER + "\n")
        for method, stat, val, p, n, note in summary_rows:
            out.write(f"{method}\t{stat}\t{val}\t{p}\t{n}\t{note}\n")

    # ---- site table (FEL preferred; fall back to MEME) ----
    site_source = None
    site_records = []
    if fel_sites:
        site_source, site_records = "FEL", fel_sites
    elif meme_sites:
        site_source, site_records = "MEME", meme_sites

    if not site_records:
        utils.write_empty(args.out_sites, SITES_HEADER)
        utils.log("no site-level method produced results; header-only site table")
    else:
        col_domain = build_codon_domain_map(codon_records, domains, hits)
        with open(args.out_sites, "w") as out:
            out.write(SITES_HEADER + "\n")
            for rec in site_records:
                col = rec["site"]
                domain = col_domain.get(col, "NA")
                p = rec["p"]
                p_str = (round(p, 6) if isinstance(p, float) and p == p else "")
                out.write(
                    f"{rec['site']}\t{col}\t{domain}\t{rec['stat']}\t{p_str}\t{rec['class']}\n"
                )
        utils.log(f"site table from {site_source}: {len(site_records)} sites")

    utils.log(f"selection: {len(summary_rows)} method row(s) written")


if __name__ == "__main__":
    main()
