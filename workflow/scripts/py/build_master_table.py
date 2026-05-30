"""Build the master per-genome table (design §3.5): the single backbone that
left-joins QC, GTDB taxonomy, species-tree tip, gene-based methanogen status,
harmonized ecology, Prometheus presence/tier, and total defense load.

Every genome in genome_metadata.tsv (the Stage-0 registry) appears exactly once;
all joins are left-joins onto that universe so a genome is never dropped for
lacking a taxonomy, tree tip, methanogen call, or ecology row. Absences are
explicit (False / 0 / NA / blank), never silently imputed — completeness,
fragmentation and detection sensitivity must remain visible as covariates.

Ecology harmonization (design §3.4): biome_envo / material_envo are filled from
the curated free-text -> ENVO map ONLY where blank in the user metadata; existing
values are never overwritten. A coarse salinity_class is derived from a numeric
salinity only when the class field is blank. Per-field *_missing flags are written
to the ecology output so downstream ecology stats run on a known, reported subset.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

MASTER_HEADER = [
    "genome_id", "representative", "tree_tip",
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
    "gtdb_family", "gtdb_genus", "gtdb_species",
    "completeness", "contamination", "quality_tier", "pass_primary", "genome_size",
    "mcrA", "n_mcr", "methanogen", "is_anme", "methanogenesis_subtype",
    "prometheus_present", "prometheus_tier", "prometheus_n_loci",
    "detect_defensefinder", "n_defense_systems",
    "isolation_source", "environment_raw", "biome_envo", "material_envo", "oxygen",
    "temperature_c", "temp_class", "salinity", "salinity_class", "ph",
    "host", "host_associated",
]

# Ecology fields tracked with explicit missingness in the ecology output.
ECOLOGY_FIELDS = [
    "isolation_source", "environment_raw", "biome_envo", "material_envo",
    "oxygen", "temperature_c", "temp_class", "salinity", "salinity_class",
    "ph", "host", "host_associated",
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


def load_envo_map(path: str):
    """free-text environment_raw -> (biome_envo, material_envo). Skip '#' lines."""
    mapping = {}
    if not path or not os.path.exists(path):
        return mapping
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except Exception:
        return mapping
    df.columns = [c.strip() for c in df.columns]
    if "environment_raw" not in df.columns:
        return mapping
    for _, row in df.iterrows():
        key = _blank(row.get("environment_raw", ""))
        if not key:
            continue
        mapping[key.lower()] = (
            _blank(row.get("biome_envo", "")),
            _blank(row.get("material_envo", "")),
        )
    return mapping


def load_ecology(path: str) -> dict:
    """metadata.tsv -> {genome_id: {field: value}} (skip '#' comment lines)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, comment="#")
    except Exception:
        return out
    if df.empty or "genome_id" not in df.columns:
        return out
    for _, row in df.iterrows():
        gid = str(row.get("genome_id", "")).strip()
        if not gid:
            continue
        out[gid] = {fld: _blank(row.get(fld, "")) for fld in ECOLOGY_FIELDS}
    return out


def derive_salinity_class(salinity_raw: str) -> str:
    """Coarse salinity class from a numeric salinity (g/L or PSU-ish).

    Conservative bins; returns '' if the value is non-numeric (we never guess).
      <0.5 -> non-saline ; 0.5-30 -> saline ; >30 -> hypersaline.
    """
    s = _blank(salinity_raw)
    if not s:
        return ""
    try:
        val = float(s)
    except (ValueError, TypeError):
        return ""
    if val < 0.5:
        return "non-saline"
    if val <= 30:
        return "saline"
    return "hypersaline"


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


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", required=True, help="genome_metadata.tsv")
    ap.add_argument("--taxonomy", required=True, help="gtdbtk.ar53.summary.tsv")
    ap.add_argument("--tree", required=True, help="species_tree.treefile")
    ap.add_argument("--methanogen", required=True, help="methanogen_calls.tsv")
    ap.add_argument("--ecology", required=True, help="metadata.tsv")
    ap.add_argument("--envo-map", required=True, help="envo_map.tsv")
    ap.add_argument("--catalog", required=True, help="prometheus_catalog.tsv")
    ap.add_argument("--sensitivity", required=True, help="detection_sensitivity.tsv")
    ap.add_argument("--reps", required=True, help="representatives.txt")
    ap.add_argument("--df-dir", required=True, help="DefenseFinder output root")
    ap.add_argument("--out-master", required=True)
    ap.add_argument("--out-ecology", required=True)
    args = ap.parse_args()

    for p in (args.out_master, args.out_ecology):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    ecology_header = ["genome_id"] + ECOLOGY_FIELDS + [f"{f}_missing" for f in ECOLOGY_FIELDS]

    # ---- genome universe (Stage-0 registry) ------------------------------ #
    try:
        meta = pd.read_csv(args.meta, sep="\t", dtype=str)
    except Exception:
        utils.write_empty(args.out_master, "\t".join(MASTER_HEADER))
        utils.write_empty(args.out_ecology, "\t".join(ecology_header))
        utils.log("Could not read --meta; wrote header-only outputs.")
        return

    if meta.empty or "genome_id" not in meta.columns:
        utils.write_empty(args.out_master, "\t".join(MASTER_HEADER))
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
    envo_map = load_envo_map(args.envo_map)
    ecology = load_ecology(args.ecology)

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

        # --- ecology (harmonized) ---
        eco = ecology.get(gid, {fld: "" for fld in ECOLOGY_FIELDS})
        biome = eco.get("biome_envo", "")
        material = eco.get("material_envo", "")
        env_raw = eco.get("environment_raw", "")
        # Fill ENVO ONLY when blank; never overwrite curated values.
        if (not biome or not material) and env_raw:
            mapped = envo_map.get(env_raw.lower())
            if mapped:
                if not biome and mapped[0]:
                    biome = mapped[0]
                if not material and mapped[1]:
                    material = mapped[1]
        eco["biome_envo"] = biome
        eco["material_envo"] = material
        # Derive salinity_class only when missing.
        if not eco.get("salinity_class", ""):
            eco["salinity_class"] = derive_salinity_class(eco.get("salinity", ""))

        master_rows.append({
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
            "isolation_source": eco.get("isolation_source", ""),
            "environment_raw": env_raw,
            "biome_envo": biome,
            "material_envo": material,
            "oxygen": eco.get("oxygen", ""),
            "temperature_c": eco.get("temperature_c", ""),
            "temp_class": eco.get("temp_class", ""),
            "salinity": eco.get("salinity", ""),
            "salinity_class": eco.get("salinity_class", ""),
            "ph": eco.get("ph", ""),
            "host": eco.get("host", ""),
            "host_associated": eco.get("host_associated", ""),
        })

        # ecology subset with per-field missingness flags
        erow = {"genome_id": gid}
        for fld in ECOLOGY_FIELDS:
            val = eco.get(fld, "")
            erow[fld] = val
            erow[f"{fld}_missing"] = (val == "")
        ecology_rows.append(erow)

    # ---- write master ---------------------------------------------------- #
    if master_rows:
        pd.DataFrame(master_rows, columns=MASTER_HEADER).to_csv(
            args.out_master, sep="\t", index=False)
    else:
        utils.write_empty(args.out_master, "\t".join(MASTER_HEADER))

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
