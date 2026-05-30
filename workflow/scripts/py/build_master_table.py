"""Build the master per-genome table (design §3.5): the single backbone that
left-joins QC, GTDB taxonomy, species-tree tip, gene-based methanogen status,
per-species niche/group presence-absence, Prometheus presence/tier, and total
defense load.

Every genome in genome_metadata.tsv (the Stage-0 registry) appears exactly once;
all joins are left-joins onto that universe so a genome is never dropped for
lacking a taxonomy, tree tip, methanogen call, or ecology row. Absences are
explicit (False / 0 / NA / blank), never silently imputed — completeness,
fragmentation and detection sensitivity must remain visible as covariates.

Ecology (design §3.4): the ecology predictors are the SPECIES-LEVEL binary
group_*/niche_* presence-absence columns produced by Phase-1 preprocessing
(species_metadata.tsv, one row per species representative, value 1 = the species
was observed in that group/niche, unioned across all its genomes). Every column
whose name starts with `group_` or `niche_` is taken as an ecology predictor,
coerced to a 0/1 int (missing -> 0), and left-joined onto the genome universe by
genome_id (genomes absent from the species table -> 0 for all ecology columns).
The legacy free-text/ENVO ecology handling is gone; --envo-map is accepted but
ignored so the rule's shell is unchanged.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

# Backbone columns (everything through n_defense_systems). The dynamic ecology
# columns (sorted group_* then niche_*) are appended at runtime to form the full
# master header; see build_master_header().
MASTER_BACKBONE = [
    "genome_id", "representative", "tree_tip",
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
    "gtdb_family", "gtdb_genus", "gtdb_species",
    "completeness", "contamination", "quality_tier", "pass_primary", "genome_size",
    "mcrA", "n_mcr", "methanogen", "is_anme", "methanogenesis_subtype",
    "prometheus_present", "prometheus_tier", "prometheus_n_loci",
    "detect_defensefinder", "n_defense_systems",
]

GTDB_RANKS = [
    ("d__", "gtdb_domain"), ("p__", "gtdb_phylum"), ("c__", "gtdb_class"),
    ("o__", "gtdb_order"), ("f__", "gtdb_family"), ("g__", "gtdb_genus"),
    ("s__", "gtdb_species"),
]


# --------------------------------------------------------------------------- #
def load_reps(path: str) -> set:
    """representatives.txt: one genome_id per line."""
    reps = set()
    if path and os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith("#"):
                    reps.add(g)
    return reps


def load_tree_tips(path: str) -> set:
    """Parse newick tip labels (species tree). Lazy dendropy; '()'/empty -> none."""
    tips = set()
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return tips
    try:
        import dendropy
    except Exception:
        utils.log("dendropy unavailable; tree_tip left blank for all genomes.")
        return tips
    try:
        tree = dendropy.Tree.get(path=path, schema="newick",
                                 preserve_underscores=True)
        for t in tree.taxon_namespace:
            if t.label:
                tips.add(t.label.strip())
    except Exception as exc:
        utils.log(f"WARNING: could not parse species tree ({exc}); tree_tip blank.")
    return tips


def parse_taxonomy(path: str) -> dict:
    """GTDB-Tk gtdbtk.ar53.summary.tsv -> {genome_id: {gtdb_domain..gtdb_species}}.

    'classification' is 'd__..;p__..;...;s__..'. Prefixes are stripped; blanks
    (or an Unclassified/N/A classification) yield empty rank fields.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return out
    df.columns = [c.strip() for c in df.columns]
    if "user_genome" not in df.columns or "classification" not in df.columns:
        return out
    for _, row in df.iterrows():
        gid = str(row.get("user_genome", "")).strip()
        if not gid:
            continue
        classification = row.get("classification", "")
        ranks = {col: "" for _, col in GTDB_RANKS}
        cstr = "" if pd.isna(classification) else str(classification).strip()
        if cstr and cstr.lower() not in ("unclassified", "n/a", "na", "none"):
            parts = [p.strip() for p in cstr.split(";")]
            for token in parts:
                for prefix, col in GTDB_RANKS:
                    if token.startswith(prefix):
                        ranks[col] = token[len(prefix):].strip()
        out[gid] = ranks
    return out


def load_methanogen(path: str) -> dict:
    """methanogen_calls.tsv -> {genome_id: {mcrA,n_mcr,methanogen,is_anme,subtype}}."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return out
    if df.empty or "genome_id" not in df.columns:
        return out
    for _, row in df.iterrows():
        gid = str(row.get("genome_id", "")).strip()
        if not gid:
            continue
        out[gid] = {
            "mcrA": _as_bool(row.get("mcrA", False)),
            "n_mcr": _as_int(row.get("n_mcr", 0)),
            "methanogen": _as_bool(row.get("methanogen", False)),
            "is_anme": _as_bool(row.get("is_anme", False)),
            "methanogenesis_subtype": _blank(row.get("methanogenesis_subtype", "NA")) or "NA",
        }
    return out


def load_sensitivity(path: str) -> dict:
    """detection_sensitivity.tsv -> {genome_id: {detect_*}} (0/1)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return out
    if df.empty or "genome_id" not in df.columns:
        return out
    for _, row in df.iterrows():
        gid = str(row.get("genome_id", "")).strip()
        if not gid:
            continue
        out[gid] = {
            "detect_defensefinder": _as_int(row.get("detect_defensefinder", 0)),
        }
    return out


def load_prometheus(path: str) -> dict:
    """prometheus_catalog.tsv -> per genome {present, n_loci, tier(best A>B>C)}.

    Counts only non-Rejected loci. tier is the best (A>B>C) among that genome's
    retained loci; '' when the genome has none.
    """
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str)
    except Exception:
        return out
    if df.empty or "genome_id" not in df.columns:
        return out
    rank = {"A": 3, "B": 2, "C": 1}
    for _, row in df.iterrows():
        gid = str(row.get("genome_id", "")).strip()
        tier = str(row.get("confidence_tier", "")).strip()
        if not gid:
            continue
        if tier == "Rejected":
            out.setdefault(gid, {"present": False, "n_loci": 0, "best": ""})
            continue
        rec = out.setdefault(gid, {"present": False, "n_loci": 0, "best": ""})
        rec["present"] = True
        rec["n_loci"] += 1
        if rank.get(tier, 0) > rank.get(rec["best"], 0):
            rec["best"] = tier
    return out


def count_defense_systems(df_dir: str, genome: str) -> int:
    """Data rows (minus header) in {df_dir}/{g}/{g}_defense_finder_systems.tsv."""
    f = os.path.join(df_dir, genome, f"{genome}_defense_finder_systems.tsv")
    if not os.path.exists(f):
        return 0
    try:
        with open(f) as fh:
            lines = [ln for ln in fh if ln.strip()]
        return max(0, len(lines) - 1)  # minus header
    except Exception:
        return 0


def load_ecology(path: str):
    """species_metadata.tsv -> (cols, {genome_id: {col: 0/1}}).

    The ecology predictors are every column whose name starts with `group_` or
    `niche_` (binary species-level presence-absence, 1 = observed in that
    group/niche). Values are coerced to 0/1 ints (missing/blank -> 0). `cols` is
    the harmonized, sorted column list (group_* before niche_*, each block
    alphabetical) used for both the master header and the ecology output.
    Returns ([], {}) when the table is unreadable, empty, lacks genome_id, or has
    no group_/niche_ columns.
    """
    if not path or not os.path.exists(path):
        return [], {}
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except Exception:
        return [], {}
    if df.empty or "genome_id" not in df.columns:
        return [], {}
    df.columns = [c.strip() for c in df.columns]
    eco_cols = ecology_columns(df.columns)
    if not eco_cols:
        return [], {}
    out = {}
    for _, row in df.iterrows():
        gid = str(row.get("genome_id", "")).strip()
        if not gid:
            continue
        out[gid] = {c: _as_int01(row.get(c, 0)) for c in eco_cols}
    return eco_cols, out


def ecology_columns(columns) -> list:
    """Sorted ecology columns: every group_*/niche_* name, group_ before niche_.

    Within each prefix block the columns are sorted alphabetically so the master
    header and ecology output are deterministic regardless of input column order.
    """
    cols = [c for c in columns if c.startswith("group_") or c.startswith("niche_")]
    groups = sorted(c for c in cols if c.startswith("group_"))
    niches = sorted(c for c in cols if c.startswith("niche_"))
    return groups + niches


def build_master_header(eco_cols) -> list:
    """Full master header = backbone columns then the ecology columns."""
    return list(MASTER_BACKBONE) + list(eco_cols)


# ---- small value helpers --------------------------------------------------- #
def _blank(v) -> str:
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none") else s


def _as_bool(v) -> bool:
    s = _blank(v).lower()
    return s in ("true", "1", "yes", "t")


def _as_int(v) -> int:
    s = _blank(v)
    if not s:
        return 0
    try:
        return int(float(s))
    except (ValueError, TypeError):
        return 0


def _as_int01(v) -> int:
    """Coerce a binary presence-absence cell to 0/1 (blank/NA/non-1 -> 0).

    Accepts 1/0, true/false, yes/no, t/f (case-insensitive) and numeric strings;
    any value that parses to a non-zero number is 1, everything else 0.
    """
    s = _blank(v).lower()
    if not s:
        return 0
    if s in ("1", "true", "yes", "t"):
        return 1
    if s in ("0", "false", "no", "f"):
        return 0
    try:
        return 1 if float(s) != 0 else 0
    except (ValueError, TypeError):
        return 0


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", required=True, help="genome_metadata.tsv")
    ap.add_argument("--taxonomy", required=True, help="gtdbtk.ar53.summary.tsv")
    ap.add_argument("--tree", required=True, help="species_tree.treefile")
    ap.add_argument("--methanogen", required=True, help="methanogen_calls.tsv")
    ap.add_argument("--ecology", required=True,
                    help="species_metadata.tsv (binary group_*/niche_* columns)")
    ap.add_argument("--envo-map", required=True,
                    help="accepted but ignored (legacy ENVO map; no longer used)")
    ap.add_argument("--catalog", required=True, help="prometheus_catalog.tsv")
    ap.add_argument("--sensitivity", required=True, help="detection_sensitivity.tsv")
    ap.add_argument("--reps", required=True, help="representatives.txt")
    ap.add_argument("--df-dir", required=True, help="DefenseFinder output root")
    ap.add_argument("--out-master", required=True)
    ap.add_argument("--out-ecology", required=True)
    args = ap.parse_args()

    for p in (args.out_master, args.out_ecology):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    # ---- ecology predictors (species-level group_*/niche_* 0/1) ---------- #
    # Detected first so header-only outputs still carry the ecology columns.
    eco_cols, ecology = load_ecology(args.ecology)
    master_header = build_master_header(eco_cols)
    ecology_header = ["genome_id"] + list(eco_cols)

    # ---- genome universe (Stage-0 registry) ------------------------------ #
    try:
        meta = pd.read_csv(args.meta, sep="\t", dtype=str)
    except Exception:
        utils.write_empty(args.out_master, "\t".join(master_header))
        utils.write_empty(args.out_ecology, "\t".join(ecology_header))
        utils.log("Could not read --meta; wrote header-only outputs.")
        return

    if meta.empty or "genome_id" not in meta.columns:
        utils.write_empty(args.out_master, "\t".join(master_header))
        utils.write_empty(args.out_ecology, "\t".join(ecology_header))
        utils.log("Empty/invalid --meta; wrote header-only outputs.")
        return

    meta.columns = [c.strip() for c in meta.columns]

    # ---- side tables ----------------------------------------------------- #
    reps = load_reps(args.reps)
    tips = load_tree_tips(args.tree)
    taxonomy = parse_taxonomy(args.taxonomy)
    methanogen = load_methanogen(args.methanogen)
    sensitivity = load_sensitivity(args.sensitivity)
    prometheus = load_prometheus(args.catalog)
    # --envo-map is accepted for CLI compatibility but no longer consulted.

    master_rows = []
    ecology_rows = []

    for _, mrow in meta.iterrows():
        gid = str(mrow["genome_id"]).strip()
        if not gid:
            continue

        # --- QC (from metadata) ---
        completeness = _blank(mrow.get("completeness", ""))
        contamination = _blank(mrow.get("contamination", ""))
        quality_tier = _blank(mrow.get("quality_tier", ""))
        pass_primary = _blank(mrow.get("pass_primary", ""))
        genome_size = _blank(mrow.get("genome_size", ""))

        # --- taxonomy ---
        tax = taxonomy.get(gid, {col: "" for _, col in GTDB_RANKS})

        # --- tree tip ---
        tree_tip = gid if gid in tips else ""

        # --- methanogen ---
        m = methanogen.get(gid, {
            "mcrA": False, "n_mcr": 0, "methanogen": False,
            "is_anme": False, "methanogenesis_subtype": "NA",
        })

        # --- prometheus ---
        prom = prometheus.get(gid, {"present": False, "n_loci": 0, "best": ""})

        # --- detection flags ---
        sens = sensitivity.get(gid, {"detect_defensefinder": 0})

        # --- defense load ---
        n_def = count_defense_systems(args.df_dir, gid)

        # --- ecology (species-level group_*/niche_* 0/1; absent genome -> 0) ---
        eco = ecology.get(gid, {})

        row = {
            "genome_id": gid,
            "representative": gid in reps,
            "tree_tip": tree_tip,
            "gtdb_domain": tax.get("gtdb_domain", ""),
            "gtdb_phylum": tax.get("gtdb_phylum", ""),
            "gtdb_class": tax.get("gtdb_class", ""),
            "gtdb_order": tax.get("gtdb_order", ""),
            "gtdb_family": tax.get("gtdb_family", ""),
            "gtdb_genus": tax.get("gtdb_genus", ""),
            "gtdb_species": tax.get("gtdb_species", ""),
            "completeness": completeness,
            "contamination": contamination,
            "quality_tier": quality_tier,
            "pass_primary": pass_primary,
            "genome_size": genome_size,
            "mcrA": m["mcrA"],
            "n_mcr": m["n_mcr"],
            "methanogen": m["methanogen"],
            "is_anme": m["is_anme"],
            "methanogenesis_subtype": m["methanogenesis_subtype"],
            "prometheus_present": prom["present"],
            "prometheus_tier": prom["best"],
            "prometheus_n_loci": prom["n_loci"],
            "detect_defensefinder": sens["detect_defensefinder"],
            "n_defense_systems": n_def,
        }
        # left-join the ecology columns (missing -> 0)
        for c in eco_cols:
            row[c] = eco.get(c, 0)
        master_rows.append(row)

        # ecology subset: genome_id + the group_*/niche_* columns
        erow = {"genome_id": gid}
        for c in eco_cols:
            erow[c] = eco.get(c, 0)
        ecology_rows.append(erow)

    # ---- write master ---------------------------------------------------- #
    if master_rows:
        pd.DataFrame(master_rows, columns=master_header).to_csv(
            args.out_master, sep="\t", index=False)
    else:
        utils.write_empty(args.out_master, "\t".join(master_header))

    # ---- write harmonized ecology --------------------------------------- #
    if ecology_rows:
        pd.DataFrame(ecology_rows, columns=ecology_header).to_csv(
            args.out_ecology, sep="\t", index=False)
    else:
        utils.write_empty(args.out_ecology, "\t".join(ecology_header))

    n_prom = sum(1 for r in master_rows if r["prometheus_present"])
    n_meth = sum(1 for r in master_rows if r["methanogen"])
    utils.log(f"master_table: {len(master_rows)} genomes "
              f"({n_prom} Prometheus-present, {n_meth} methanogen).")


if __name__ == "__main__":
    main()
