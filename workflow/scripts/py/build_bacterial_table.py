"""Build a slim bacterial master table for the cross-domain comparative analyses.

One row per bacterial genome of the analysis universe (--samples; else qc ∪ taxonomy),
with the columns the reused archaeal scripts need:
  genome_id, tree_tip, gtdb_domain..gtdb_species, completeness, contamination,
  genome_size, pass_primary, prometheus_present, prometheus_tier, prometheus_n_loci,
  n_defense_systems, n_defense_systems_excl_prometheus,
  + additive: detect_defensefinder, prometheus_status_primary, prometheus_status_sensitivity

`compute_prevalence.py`, `cooccurrence.R` and `ecology_stats.R` are reused unchanged on
this table. Every definition is the archaeal master's (build_master_table.py, Contract
v2), re-implemented here so the bacterial branch does not depend on that script's
internal helpers:
  * prometheus_tier = best NON-REJECTED catalog tier (A > B > C; '' when none),
    prometheus_present = has >= 1 non-rejected locus, prometheus_n_loci = their count;
    a tier outside {A, B, C, Rejected} is an error;
  * prometheus_status_{primary,sensitivity} = Contract-v2 presence rule for the tier
    sets --primary-tiers / --sensitivity-tiers: present (tier in T) | absent (tier '')
    | excluded (validated, best tier outside T);
  * pass_primary = completeness >= --min-completeness AND contamination <=
    --max-contamination (the archaeal analysis set's quality bar: preprocess keeps
    completeness >= 90, contamination <= 5), so prevalence is compared at matched
    quality. An unknown completeness (or, when the QC table has a contamination column,
    an unknown contamination) does not pass;
  * n_defense_systems = DefenseFinder systems with activity != Antidefense;
    n_defense_systems_excl_prometheus = the same minus systems whose type EQUALS
    --system (case-insensitive) -- the focal system is never inside its own covariate;
  * tree_tip = genome_id (prune_tree_to_taxa.py relabels the GTDB tips to genome_id).
Taxonomy may be a GTDB-Tk summary (user_genome + classification), a genome_id +
classification table, or genome_id + gtdb_* columns. Duplicated ids are an error, and so
is a missing per-genome DefenseFinder table (extract_bacterial_prometheus.py writes one,
header-only when empty, for every genome of the universe).
"""
from __future__ import annotations
import argparse
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

GTDB_RANKS = [("d__", "gtdb_domain"), ("p__", "gtdb_phylum"), ("c__", "gtdb_class"),
              ("o__", "gtdb_order"), ("f__", "gtdb_family"), ("g__", "gtdb_genus"),
              ("s__", "gtdb_species")]
RANK_COLS = [col for _, col in GTDB_RANKS]
HEADER = ["genome_id", "tree_tip", *RANK_COLS,
          "completeness", "contamination", "genome_size", "pass_primary",
          "prometheus_present", "prometheus_tier", "prometheus_n_loci",
          "n_defense_systems", "n_defense_systems_excl_prometheus",
          # additive
          "detect_defensefinder", "prometheus_status_primary", "prometheus_status_sensitivity"]
SIZE_COLS = ("genome_size", "length", "total_length", "genome_length", "size", "contig_bp")
TIER_RANK = {"A": 3, "B": 2, "C": 1}
VALID_TIERS = {"A", "B", "C", "Rejected"}


class InputError(Exception):
    pass


def _blank(v) -> str:
    s = "" if v is None else str(v).strip()
    return "" if s.lower() in ("nan", "none", "na", "<na>") else s


def _read(path, what):
    import pandas as pd
    if not path or not os.path.exists(path):
        raise InputError(f"{what} not found: {path!r}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    df.columns = [c.strip() for c in df.columns]
    return df


def _id_col(df):
    for c in ("genome_id", "user_genome"):
        if c in df.columns:
            return c
    return df.columns[0]


def _ids(df, col, what):
    ids = [str(g).strip() for g in df[col]]
    nonblank = [g for g in ids if g]
    seen, dup = set(), set()
    for g in nonblank:
        if g in seen:
            dup.add(g)
        seen.add(g)
    if dup:
        raise InputError(f"{len(dup)} duplicated genome id(s) in {what}: {sorted(dup)[:5]}")
    return ids


def parse_tiers(spec):
    out = {t.strip() for t in str(spec or "").split(",") if t.strip()}
    if not out or out - set(TIER_RANK):
        raise InputError(f"tier set {spec!r} must be a non-empty subset of A,B,C")
    return out


def tier_status(tier, tiers):
    """Contract v2 item 6: present | absent | excluded."""
    if tier == "":
        return "absent"
    return "present" if tier in tiers else "excluded"


def load_qc(path):
    """genome_id -> {completeness, contamination, genome_size}; plus has_contamination."""
    df = _read(path, "QC table")
    if df.shape[1] == 0:
        return {}, False
    low = {c.lower(): c for c in df.columns}
    comp = next((low[c] for c in ("completeness", "comp") if c in low), None)
    cont = next((low[c] for c in ("contamination", "contam") if c in low), None)
    size = next((low[c] for c in SIZE_COLS if c in low), None)
    if comp is None:
        raise InputError(f"QC table {path} has no completeness column")
    out = {}
    for gid, (_, r) in zip(_ids(df, _id_col(df), f"QC table {path}"), df.iterrows()):
        if gid:
            out[gid] = {"completeness": _blank(r[comp]),
                        "contamination": _blank(r[cont]) if cont else "",
                        "genome_size": _blank(r[size]) if size else ""}
    return out, cont is not None


def load_taxonomy(path):
    """genome_id -> {gtdb_domain..gtdb_species} from any of the supported layouts."""
    df = _read(path, "taxonomy")
    if df.shape[1] == 0:
        return {}
    ccol = next((c for c in df.columns if c.lower() == "classification"), None)
    have_ranks = [c for c in RANK_COLS if c in df.columns]
    if ccol is None and not have_ranks:
        raise InputError(f"taxonomy {path} has neither a 'classification' column nor gtdb_* "
                         f"rank columns")
    out = {}
    for gid, (_, r) in zip(_ids(df, _id_col(df), f"taxonomy {path}"), df.iterrows()):
        if not gid:
            continue
        ranks = {col: "" for col in RANK_COLS}
        cstr = _blank(r[ccol]) if ccol else ""
        if cstr and not cstr.lower().startswith("unclassified") and cstr.lower() not in ("n/a", "na"):
            for token in (p.strip() for p in cstr.split(";")):
                for pre, col in GTDB_RANKS:
                    if token.startswith(pre):
                        ranks[col] = token[len(pre):].strip()
        for col in have_ranks:                      # explicit gtdb_* columns win
            v = _blank(r[col])
            pre = dict((c, p) for p, c in GTDB_RANKS)[col]
            if v.startswith(pre):
                v = v[len(pre):].strip()
            if v:
                ranks[col] = v
        out[gid] = ranks
    return out


def load_samples(path):
    df = _read(path, "samples")
    if df.shape[1] == 0:
        return []
    return [g for g in _ids(df, _id_col(df), f"samples {path}") if g]


def load_catalog(path):
    """-> ({genome: {present, n_loci, best}} over non-rejected loci, detected genomes)."""
    df = _read(path, "bacterial Prometheus catalog")
    out, detected = {}, set()
    if df.empty:
        return out, detected
    for col in ("genome_id", "confidence_tier"):
        if col not in df.columns:
            raise InputError(f"{path} has no {col} column")
    for gid, tier in zip(df["genome_id"], df["confidence_tier"]):
        gid, tier = gid.strip(), tier.strip()
        if not gid:
            continue
        if tier not in VALID_TIERS:
            raise InputError(f"catalog tier {tier!r} for {gid} is not one of {sorted(VALID_TIERS)}")
        detected.add(gid)
        rec = out.setdefault(gid, {"present": False, "n_loci": 0, "best": ""})
        if tier == "Rejected":
            continue
        rec["present"] = True
        rec["n_loci"] += 1
        if TIER_RANK[tier] > TIER_RANK.get(rec["best"], 0):
            rec["best"] = tier
    return out, detected


def count_defense(df_dir, genome, systems):
    """(n_defense, n_defense_excl_focal) from {df_dir}/{g}/{g}_defense_finder_systems.tsv."""
    f = os.path.join(df_dir, genome, f"{genome}_defense_finder_systems.tsv")
    df = _read(f, f"DefenseFinder systems table for {genome}")
    if df.empty:
        return 0, 0
    if "type" not in df.columns:
        raise InputError(f"{f} has no 'type' column")
    if "activity" in df.columns:
        df = df[df["activity"].str.strip().str.casefold() != "antidefense"]
    n_total = int(len(df))
    focal = df["type"].str.strip().str.casefold().isin(systems)
    return n_total, int(n_total - focal.sum())


def assembly_size(genomes_dir, genome_id, ext):
    """Total residue count of an assembly (streamed; sequence chars only). '' if absent."""
    if not genomes_dir:
        return ""
    exts = (ext, "fasta", "fna", "fa", "fasta.gz", "fna.gz", "fa.gz")
    cands = [os.path.join(genomes_dir, f"{genome_id}.{e}") for e in exts]
    path = next((p for p in cands if os.path.exists(p)), None)
    if path is None:
        return ""
    n = 0
    # gzip-aware: bacteria.genome_extension may be e.g. 'fna.gz' (geNomad reads it too)
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n if n else ""


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--taxonomy", required=True)
    ap.add_argument("--qc", required=True)
    ap.add_argument("--catalog", required=True, help="bacterial_prometheus_catalog.tsv")
    ap.add_argument("--defense-dir", required=True, help="split per-genome systems dir")
    ap.add_argument("--samples", default="",
                    help="optional analysis-universe TSV (genome_id); default qc ∪ taxonomy")
    ap.add_argument("--genomes-dir", default="", help="optional assemblies dir (for genome_size)")
    ap.add_argument("--genome-ext", default="fasta")
    ap.add_argument("--min-completeness", type=float, default=90.0)
    ap.add_argument("--max-contamination", type=float, default=None,
                    help="pass_primary also needs contamination <= this (default: no bound)")
    ap.add_argument("--system", default="Prometheus",
                    help="focal DefenseFinder type(s), comma-separated, exact match")
    ap.add_argument("--primary-tiers", default="A")
    ap.add_argument("--sensitivity-tiers", default="A,B")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    try:
        run(args)
    except InputError as exc:
        utils.log(f"FATAL: {exc}")
        sys.exit(1)


def run(args) -> None:
    import pandas as pd
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    primary, sensitivity = parse_tiers(args.primary_tiers), parse_tiers(args.sensitivity_tiers)
    systems = {s.strip().casefold() for s in args.system.split(",") if s.strip()}

    qc, has_contam = load_qc(args.qc)
    taxonomy = load_taxonomy(args.taxonomy)
    prometheus, detected = load_catalog(args.catalog)
    if args.max_contamination is not None and not has_contam:
        utils.log(f"WARNING: --max-contamination {args.max_contamination} requested but the QC "
                  f"table has no contamination column; pass_primary uses completeness only")
    genomes = load_samples(args.samples) if args.samples else \
        list(dict.fromkeys(list(qc.keys()) + list(taxonomy.keys())))
    if not genomes:
        utils.write_empty(args.out, "\t".join(HEADER))
        utils.log("No bacterial genomes; header-only table.")
        return
    gset = set(genomes)
    no_qc = [g for g in genomes if g not in qc]
    no_tax = [g for g in genomes if g not in taxonomy]
    if no_qc:
        utils.log(f"WARNING: {len(no_qc)} genome(s) lack QC (pass_primary False), e.g. {no_qc[:3]}")
    if no_tax:
        utils.log(f"WARNING: {len(no_tax)} genome(s) lack taxonomy (no taxon prevalence), "
                  f"e.g. {no_tax[:3]}")
    orphans = sorted(detected - gset)
    if orphans:
        raise InputError(f"{len(orphans)} catalog genome(s) are not in the bacterial genome "
                         f"list (e.g. {orphans[:3]}): ids differ between the inputs")

    rows = []
    for gid in genomes:
        tax = taxonomy.get(gid, {col: "" for col in RANK_COLS})
        q = qc.get(gid, {"completeness": "", "contamination": "", "genome_size": ""})
        comp, cont = _float(q["completeness"]), _float(q["contamination"])
        pass_primary = comp is not None and comp >= args.min_completeness
        if pass_primary and args.max_contamination is not None and has_contam:
            pass_primary = cont is not None and cont <= args.max_contamination
        prom = prometheus.get(gid, {"present": False, "n_loci": 0, "best": ""})
        n_def, n_def_excl = count_defense(args.defense_dir, gid, systems)
        gsize = q.get("genome_size", "") or assembly_size(args.genomes_dir, gid, args.genome_ext)
        rows.append({
            "genome_id": gid, "tree_tip": gid,
            **{col: tax.get(col, "") for col in RANK_COLS},
            "completeness": q["completeness"], "contamination": q["contamination"],
            "genome_size": gsize, "pass_primary": bool(pass_primary),
            "prometheus_present": bool(prom["present"]), "prometheus_tier": prom["best"],
            "prometheus_n_loci": prom["n_loci"],
            "n_defense_systems": n_def, "n_defense_systems_excl_prometheus": n_def_excl,
            "detect_defensefinder": int(gid in detected),
            "prometheus_status_primary": tier_status(prom["best"], primary),
            "prometheus_status_sensitivity": tier_status(prom["best"], sensitivity),
        })

    pd.DataFrame(rows, columns=HEADER).to_csv(args.out, sep="\t", index=False)
    n_prom = sum(1 for r in rows if r["prometheus_present"])
    n_pass = sum(1 for r in rows if r["pass_primary"])
    n_exc = sum(1 for r in rows if r["prometheus_status_primary"] == "excluded")
    utils.log(f"bacterial table: {len(rows)} genomes ({n_pass} pass_primary, {n_prom} with a "
              f"validated Prometheus locus, {n_exc} excluded from the tier-"
              f"{'+'.join(sorted(primary))} analyses).")


if __name__ == "__main__":
    main()
