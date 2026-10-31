"""Prune a large reference tree (e.g. GTDB bac120) to a focal genome set.

For the bacterial co-occurrence / niche / viral tests we need a species tree; the
defensible, tractable choice is to PRUNE the published GTDB bac120 reference tree to
the focal genomes. GTDB tips are species-REPRESENTATIVE accessions with a source prefix
(`RS_GCF_000005845.2`, `GB_GCA_...`), internal labels are quoted support+taxonomy
strings ('100.0:p__Pseudomonadota'), and user genome ids often lack the prefix, carry
another accession version, the paired GenBank/RefSeq accession, or are not species
representatives at all. Each focal genome is therefore matched to ONE tip by the first
tier that resolves it uniquely:

  exact               label == genome_id
  prefix_stripped     after removing RS_/GB_ from both
  version_stripped    ... and a trailing accession version (.1, .2 ...)
  paired_accession    same NCBI assembly number, GCA_ <-> GCF_ (paired assemblies)
  species_representative
                      (needs --taxonomy and --ref-taxonomy) the tip that represents the
                      genome's GTDB species (s__ name)

A key that maps to >1 tip is 'ambiguous' and not used. When several focal genomes match
the same tip (e.g. conspecific non-representatives), the best match tier wins, then the
highest quality (completeness - 5 x contamination, from --quality), then the smallest
id; the others are 'collapsed' (one genome per species-tree tip, as in the archaeal
species-representative design). Kept tips are RELABELLED to the focal genome_id, so the
pruned tree's tips equal the master table's tree_tip, internal labels are dropped
(--keep-internal-labels keeps them) and branch lengths/rooting are preserved.

Every genome's outcome is written to --out-report. Fewer matched genomes than
--min-match-fraction of the focal set is an ERROR (exit 1), never a silently tiny tree.
An empty focal set yields a trivial '();' tree.
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

REPORT_HEADER = "genome_id\ttree_tip\tmatch\tnote"
TIERS = ["exact", "prefix_stripped", "version_stripped", "paired_accession",
         "species_representative"]
SAFE_ID = re.compile(r"^[A-Za-z0-9_.\-]+$")


def key_exact(s):
    return s


def key_prefix(s):
    return re.sub(r"^(RS_|GB_)", "", s)


def key_version(s):
    return re.sub(r"\.\d+$", "", key_prefix(s))


def key_paired(s):
    m = re.match(r"^(?:RS_|GB_)?GC[AF]_(\d{6,})(?:\.\d+)?$", s)
    return f"GCX_{m.group(1)}" if m else None


KEYS = {"exact": key_exact, "prefix_stripped": key_prefix,
        "version_stripped": key_version, "paired_accession": key_paired}


def _id_col(df, preferred=""):
    for c in (preferred, "genome_id", "user_genome"):
        if c and c in df.columns:
            return c
    return df.columns[0]


def load_ids(path, col=""):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    if df.shape[1] == 0:
        return []
    c = _id_col(df, col)
    return list(dict.fromkeys(str(g).strip() for g in df[c].dropna() if str(g).strip()))


def load_quality(path):
    """genome_id -> completeness - 5*contamination (missing -> -inf)."""
    import pandas as pd
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    low = {c.lower(): c for c in df.columns}
    gcol = _id_col(df)
    comp, cont = low.get("completeness"), low.get("contamination")
    out = {}
    for _, r in df.iterrows():
        g = str(r[gcol]).strip()
        try:
            q = float(r[comp]) if comp else float("-inf")
            if cont and str(r[cont]).strip() not in ("", "nan"):
                q -= 5.0 * float(r[cont])
        except (TypeError, ValueError):
            q = float("-inf")
        out[g] = q
    return out


def species_of(classification):
    for tok in str(classification or "").split(";"):
        tok = tok.strip()
        if tok.startswith("s__") and tok[3:].strip():
            return tok[3:].strip()
    return ""


def load_user_species(path):
    """genome_id -> GTDB species name from a classification or gtdb_species column."""
    import pandas as pd
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    gcol = _id_col(df)
    ccol = next((c for c in df.columns if c.lower() == "classification"), None)
    scol = "gtdb_species" if "gtdb_species" in df.columns else None
    out = {}
    for _, r in df.iterrows():
        g = str(r[gcol]).strip()
        sp = species_of(r[ccol]) if ccol else ""
        if not sp and scol and str(r[scol]).strip() not in ("", "nan"):
            sp = re.sub(r"^s__", "", str(r[scol]).strip())
        if g and sp:
            out[g] = sp
    return out


def load_ref_taxonomy(path):
    """GTDB {bac120,ar53}_taxonomy.tsv (accession<TAB>classification, no header)."""
    out = {}
    if not path:
        return out
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] and not parts[0].startswith("#"):
                out[parts[0].strip()] = species_of(parts[1])
    return out


def match_genomes(keep, tips, user_species=None, ref_tax=None, tiers=None, quality=None):
    """Return {genome: (tip, tier)} plus a report {genome: (tip, status, note)}."""
    tiers = tiers or TIERS
    quality = quality or {}
    assigned, report = {}, {}
    pending = list(keep)
    for tier in tiers:
        if tier == "species_representative":
            if not (user_species is not None and ref_tax):
                continue
            tip_sp = {}
            for t in tips:
                sp = ref_tax.get(t) or ref_tax.get(key_prefix(t), "")
                if sp:
                    tip_sp.setdefault(sp, []).append(t)
            ref_by_key = {}
            for acc, sp in ref_tax.items():
                for k in (key_prefix(acc), key_version(acc)):
                    ref_by_key.setdefault(k, set()).add(sp)
            for g in pending:
                sp = user_species.get(g, "")
                if not sp:
                    cand = ref_by_key.get(key_prefix(g)) or ref_by_key.get(key_version(g)) or set()
                    sp = next(iter(cand)) if len(cand) == 1 else ""
                hits = tip_sp.get(sp, []) if sp else []
                if len(hits) == 1:
                    assigned[g] = (hits[0], tier)
                elif len(hits) > 1:
                    report[g] = ("", "ambiguous", f"species {sp} has {len(hits)} tips")
        else:
            fn = KEYS[tier]
            idx = {}
            for t in tips:
                k = fn(t)
                if k:
                    idx.setdefault(k, []).append(t)
            for g in pending:
                k = fn(g)
                hits = idx.get(k, []) if k else []
                if len(hits) == 1:
                    assigned[g] = (hits[0], tier)
                elif len(hits) > 1 and g not in report:
                    report[g] = ("", "ambiguous", f"{tier} key '{k}' matches {len(hits)} tips: "
                                                  f"{hits[:3]}")
        pending = [g for g in pending if g not in assigned]
    # one genome per tip: best tier, then quality, then id
    by_tip = {}
    for g, (t, tier) in assigned.items():
        by_tip.setdefault(t, []).append(g)
    for t, gs in by_tip.items():
        gs.sort(key=lambda g: (TIERS.index(assigned[g][1]), -quality.get(g, float("-inf")), g))
        for loser in gs[1:]:
            report[loser] = (t, "collapsed", f"tip already represents {gs[0]}")
            del assigned[loser]
    for g, (t, tier) in assigned.items():
        report[g] = (t, tier, "")
    for g in keep:
        report.setdefault(g, ("", "unmatched", "no tree tip under any tier"))
    return assigned, report


def write_report(path, keep, report):
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(REPORT_HEADER + "\n")
        for g in keep:
            t, status, note = report[g]
            fh.write(f"{g}\t{t}\t{status}\t{note}\n")


def write_trivial(path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("();\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--in-tree", required=True)
    ap.add_argument("--keep", default="", help="TSV with a genome_id column (or 1st column)")
    ap.add_argument("--keep-col", default="", help="id column of --keep (default genome_id)")
    ap.add_argument("--keep-fasta", default="", help="alternatively: keep the ids of this FASTA")
    ap.add_argument("--out-tree", required=True)
    ap.add_argument("--out-report", default="")
    ap.add_argument("--match", default="all", choices=["all", "exact"],
                    help="'exact' = label equality only (gene trees)")
    ap.add_argument("--taxonomy", default="", help="focal genome taxonomy (genome_id + "
                    "classification or gtdb_species) for species-representative matching")
    ap.add_argument("--ref-taxonomy", default="", help="GTDB bac120_taxonomy.tsv of the tree release")
    ap.add_argument("--quality", default="", help="QC TSV (completeness/contamination) to "
                    "choose among conspecific genomes")
    ap.add_argument("--min-match-fraction", type=float, default=0.9)
    ap.add_argument("--no-relabel", action="store_true", help="keep the reference tip labels")
    ap.add_argument("--keep-internal-labels", action="store_true")
    args = ap.parse_args()

    import dendropy
    if args.keep_fasta:
        keep = list(dict.fromkeys(utils.fasta_ids(args.keep_fasta)))
    elif args.keep:
        keep = load_ids(args.keep, args.keep_col)
    else:
        ap.error("one of --keep / --keep-fasta is required")
    if not keep:
        write_trivial(args.out_tree)
        write_report(args.out_report, [], {})
        utils.log("prune_tree: empty focal set; wrote trivial tree.")
        return
    unsafe = [g for g in keep if not SAFE_ID.match(g)]
    if unsafe:
        utils.log(f"FATAL: focal id(s) not newick-safe ([A-Za-z0-9_.-] only): {unsafe[:3]}")
        sys.exit(1)
    if not os.path.exists(args.in_tree) or os.path.getsize(args.in_tree) == 0:
        utils.log(f"FATAL: reference tree {args.in_tree} is missing or empty")
        sys.exit(1)

    tree = dendropy.Tree.get(path=args.in_tree, schema="newick", preserve_underscores=True,
                             rooting="force-rooted")
    tips = [lf.taxon.label for lf in tree.leaf_node_iter() if lf.taxon is not None]
    dup = sorted(t for t, n in Counter(tips).items() if n > 1)
    if dup:
        utils.log(f"FATAL: duplicated tip labels in the reference tree: {dup[:3]}")
        sys.exit(1)
    tiers = ["exact"] if args.match == "exact" else TIERS
    user_species = load_user_species(args.taxonomy) if args.taxonomy else None
    ref_tax = load_ref_taxonomy(args.ref_taxonomy)
    assigned, report = match_genomes(keep, tips, user_species, ref_tax, tiers,
                                     load_quality(args.quality))
    write_report(args.out_report, keep, report)

    counts = {}
    for _, status, _n in report.values():
        counts[status] = counts.get(status, 0) + 1
    frac = len(assigned) / len(keep)
    utils.log(f"prune_tree: {len(keep)} focal genomes vs {len(tips)} reference tips -> "
              f"{len(assigned)} on the tree ({frac:.1%}); outcome counts {counts}")
    if frac < args.min_match_fraction:
        ex_g = [g for g in keep if report[g][1] in ("unmatched", "ambiguous")][:3]
        utils.log(f"FATAL: only {len(assigned)}/{len(keep)} focal genomes ({frac:.1%}) map to "
                  f"the reference tree (minimum {args.min_match_fraction:.0%}). Unmatched e.g. "
                  f"{ex_g}; tree tips look like {tips[:3]}. Use GTDB accessions as genome_id, "
                  f"supply bacteria.gtdb_bac_taxonomy for species-representative matching, or "
                  f"lower bacteria.min_tree_match_fraction deliberately.")
        sys.exit(1)
    if not assigned:
        utils.log("FATAL: no focal genome is on the reference tree")
        sys.exit(1)

    tip2genome = {t: g for g, (t, _tier) in assigned.items()}
    tree.retain_taxa_with_labels(set(tip2genome))
    if not args.no_relabel:
        for lf in tree.leaf_node_iter():
            lf.taxon.label = tip2genome[lf.taxon.label]
    tree.write(path=args.out_tree, schema="newick", unquoted_underscores=True,
               suppress_rooting=True,
               suppress_internal_node_labels=not args.keep_internal_labels)
    n = sum(1 for _ in tree.leaf_node_iter())
    utils.log(f"prune_tree: kept {n} tips -> {args.out_tree}")


if __name__ == "__main__":
    main()
