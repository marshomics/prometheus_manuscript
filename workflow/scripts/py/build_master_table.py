"""Build the master per-genome table (design §3.5): the single backbone that
left-joins QC, GTDB taxonomy, species-tree tip, gene-based methanogen status,
species-level group/niche presence-absence, Prometheus presence/tier and defense load.

Universe and joins
------------------
The universe is genome_metadata.tsv (the Stage-0 registry: the species
representatives, i.e. the genomes of the species table). Every genome appears exactly
once (a duplicated genome_id is an error) and every join is a left join onto it. All
side tables are declared rule inputs: one that is missing or unreadable is an error,
never an empty join, and one that names a genome outside the universe is an error (it
would come from a different genome set). Genomes absent from the methanogen calls, the
Prometheus catalog or the sensitivity table take the explicit no-evidence values
(False / 0 / ''); completeness, fragmentation and detection sensitivity stay visible as
covariates.

Prometheus (Contract v2 item 6, rev 2)
--------------------------------------
prometheus_tier     best NON-REJECTED tier among the genome's catalog loci (A > B > C);
                    '' when it has none.
prometheus_present  True iff the genome has >= 1 non-rejected locus.
prometheus_n_loci   number of non-rejected loci.
A catalog tier outside {A, B, C, Rejected} is an error. For a tier set T a genome is
PRESENT if prometheus_tier is in T, ABSENT if prometheus_tier == '', and EXCLUDED
otherwise (a validated locus whose best tier is outside T is not a confident absence).
The additive columns prometheus_status_primary (--primary-tiers, fed by
stats.confidence_tiers_primary) and prometheus_status_sensitivity (--sensitivity-tiers,
default A,B) carry that status as present | absent | excluded.

Defense load (DefenseFinder 2.0 *_defense_finder_systems.tsv)
-------------------------------------------------------------
Columns sys_id, type, subtype, activity, sys_beg, sys_end, protein_in_syst,
genes_count, name_of_profiles_in_sys (defense_finder_posttreat/df_systems.py, v2.0.0).
With -a DefenseFinder also reports AntiDefenseFinder systems in the same table
(best_solution.py: activity 'Antidefense' for ADF models, 'Defense' otherwise).
n_defense_systems                  systems whose activity is not 'Antidefense'.
n_defense_systems_excl_prometheus  the same minus systems whose type is EXACTLY
                                   --system (case-insensitive equality, no substring;
                                   defense-finder-models 2.0 defines the type as
                                   'Prometheus'), so the focal system is never inside
                                   the defense-richness covariate.

Ecology (design §3.4)
---------------------
The predictors are the species-level binary group_* / niche_* columns of the species
table (Phase 1 unions them over each species' genomes). Each cell is 1 (1, true, yes,
t, any non-zero number; non-binary numbers are counted in the log) or 0 (0, false, no,
f); a blank / NA cell stays EMPTY (NA): an unknown niche is not an absence. Every
universe genome must be in the species table. --out-missingness reports n_1 / n_0 /
n_NA per column. The legacy free-text ENVO harmonisation (config ecology.envo_map,
ecology.predictors) does not apply to these curated binary columns; --envo-map is
accepted for CLI compatibility, ignored, and that is logged.

Column order: the SCHEMA backbone (unchanged), the ecology block (group_* sorted, then
niche_* sorted), then additive columns (ADDITIVE).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

MASTER_BACKBONE = [
    "genome_id", "representative", "tree_tip",
    "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
    "gtdb_family", "gtdb_genus", "gtdb_species",
    "completeness", "contamination", "quality_tier", "pass_primary", "genome_size",
    "mcrA", "n_mcr", "methanogen", "is_anme", "methanogenesis_subtype",
    "prometheus_present", "prometheus_tier", "prometheus_n_loci",
    "detect_defensefinder", "n_defense_systems", "n_defense_systems_excl_prometheus",
]
ADDITIVE = ["prometheus_status_primary", "prometheus_status_sensitivity",
            "n_contigs", "n50", "gc", "tree_status", "mcr_group"]
# Methanogenesis pathway capacity (--pathways, optional; scaffolds.smk methanogenesis_pathways).
# Two-sided calls become 1 / 0 / NA, where NA covers 'unresolved' and 'not_applicable' (not a
# methanogen); one-sided flags are 1 / NA ('not_detected' is not evidence of absence).
MGEN_PATHWAYS = ["co2_reducing", "acetoclastic", "methylotrophic", "methyl_dismutating",
                 "methyl_reducing"]
MGEN_SUBS_2 = ["sub_h2", "sub_formate", "sub_methanol", "sub_methylamine"]
MGEN_SUBS_1 = ["sub_methylthiol", "sub_quatamine", "sub_methoxyl"]
MGEN_ADDITIVE = ([f"mgen_{p}" for p in MGEN_PATHWAYS] +
                 [f"mgen_{x}" for x in MGEN_SUBS_2 + MGEN_SUBS_1] +
                 ["subtype_status", "subtype_taxonomy_concordance"])

GTDB_RANKS = [
    ("d__", "gtdb_domain"), ("p__", "gtdb_phylum"), ("c__", "gtdb_class"),
    ("o__", "gtdb_order"), ("f__", "gtdb_family"), ("g__", "gtdb_genus"),
    ("s__", "gtdb_species"),
]
TIER_RANK = {"A": 3, "B": 2, "C": 1}
VALID_TIERS = {"A", "B", "C", "Rejected"}


# ---- value helpers --------------------------------------------------------- #
def _blank(v) -> str:
    if v is None:
        return ""
    try:
        if pd.isna(v):
            return ""
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return "" if s.lower() in ("nan", "none", "na", "<na>") else s


def _as_bool(v) -> bool:
    return _blank(v).lower() in ("true", "1", "yes", "t")


def _as_int(v) -> int:
    s = _blank(v)
    if not s:
        return 0
    return int(float(s))


def _eco_value(v):
    """-> (value, coerced): 1 / 0 / '' (NA) and whether a non-binary number was coerced."""
    s = _blank(v).lower()
    if not s:
        return "", False
    if s in ("1", "true", "yes", "t", "1.0"):
        return 1, False
    if s in ("0", "false", "no", "f", "0.0"):
        return 0, False
    try:
        x = float(s)
    except ValueError:
        raise ValueError(f"ecology value {v!r} is not binary")
    return (1 if x != 0 else 0), True


def _read_tsv(path: str, what: str) -> pd.DataFrame:
    if not path or not os.path.exists(path):
        raise FileNotFoundError(f"{what} not found: {path!r}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None)
    df.columns = [c.strip() for c in df.columns]
    return df


def tier_status(tier: str, tiers: set) -> str:
    """Contract v2 item 6: present | absent | excluded for one tier set."""
    if tier == "":
        return "absent"
    return "present" if tier in tiers else "excluded"


def parse_tiers(spec: str) -> set:
    out = {t.strip() for t in (spec or "").split(",") if t.strip()}
    bad = out - set(TIER_RANK)
    if not out or bad:
        raise ValueError(f"tier set {spec!r} must be a non-empty subset of A,B,C")
    return out


# ---- side tables ----------------------------------------------------------- #
def load_reps(path: str) -> set:
    """representatives.txt: one genome_id per line (a 'genome_id' header is skipped)."""
    reps = set()
    with open(path) as fh:
        for i, line in enumerate(fh):
            g = line.strip()
            if not g or g.startswith("#") or (i == 0 and g == "genome_id"):
                continue
            reps.add(g)
    return reps


def load_tree_tips(path: str) -> set:
    """Newick tip labels; '();' / empty -> none. dendropy is required."""
    import dendropy
    with open(path) as fh:
        raw = fh.read().strip()
    if raw in ("", ";", "();"):
        return set()
    tree = dendropy.Tree.get(data=raw, schema="newick", preserve_underscores=True)
    tips = [lf.taxon.label.strip() for lf in tree.leaf_node_iter() if lf.taxon is not None]
    if len(set(tips)) != len(tips):
        raise ValueError("species tree has duplicate tip labels")
    return set(tips)


def parse_taxonomy(path: str) -> dict:
    """taxonomy.tsv (user_genome, classification) -> {genome_id: {gtdb_*: name}}."""
    df = _read_tsv(path, "taxonomy")
    if df.empty:
        return {}
    if "user_genome" not in df.columns or "classification" not in df.columns:
        raise ValueError(f"taxonomy {path} lacks user_genome/classification")
    out = {}
    for gid, cls in zip(df["user_genome"], df["classification"]):
        gid = gid.strip()
        if not gid:
            continue
        ranks = {col: "" for _, col in GTDB_RANKS}
        cstr = cls.strip()
        if cstr and not cstr.lower().startswith("unclassified") and cstr.lower() not in ("n/a", "na"):
            for token in (p.strip() for p in cstr.split(";")):
                for prefix, col in GTDB_RANKS:
                    if token.startswith(prefix):
                        ranks[col] = token[len(prefix):].strip()
        out[gid] = ranks
    return out


def load_methanogen(path: str) -> dict:
    df = _read_tsv(path, "methanogen calls")
    if df.empty:
        return {}
    if "genome_id" not in df.columns:
        raise ValueError(f"{path} has no genome_id column")
    out = {}
    for _, row in df.iterrows():
        gid = row["genome_id"].strip()
        if not gid:
            continue
        out[gid] = {
            "mcrA": _as_bool(row.get("mcrA", "")),
            "n_mcr": _as_int(row.get("n_mcr", "")),
            "methanogen": _as_bool(row.get("methanogen", "")),
            "is_anme": _as_bool(row.get("is_anme", "")),
            "methanogenesis_subtype": _blank(row.get("methanogenesis_subtype", "")) or "NA",
            "mcr_group": _blank(row.get("mcr_group", "")) or "NA",
        }
    return out


def load_sensitivity(path: str) -> dict:
    df = _read_tsv(path, "detection sensitivity")
    if df.empty:
        return {}
    if "genome_id" not in df.columns:
        raise ValueError(f"{path} has no genome_id column")
    return {r["genome_id"].strip(): {"detect_defensefinder": _as_int(r.get("detect_defensefinder", ""))}
            for _, r in df.iterrows() if r["genome_id"].strip()}


def load_prometheus(path: str) -> dict:
    """prometheus_catalog.tsv -> {genome: {present, n_loci, best}} over NON-rejected loci."""
    df = _read_tsv(path, "prometheus catalog")
    out = {}
    if df.empty:
        return out
    for col in ("genome_id", "confidence_tier"):
        if col not in df.columns:
            raise ValueError(f"{path} has no {col} column")
    for gid, tier in zip(df["genome_id"], df["confidence_tier"]):
        gid, tier = gid.strip(), tier.strip()
        if not gid:
            continue
        if tier not in VALID_TIERS:
            raise ValueError(f"catalog tier {tier!r} for {gid} is not one of {sorted(VALID_TIERS)}")
        rec = out.setdefault(gid, {"present": False, "n_loci": 0, "best": ""})
        if tier == "Rejected":
            continue
        rec["present"] = True
        rec["n_loci"] += 1
        if TIER_RANK[tier] > TIER_RANK.get(rec["best"], 0):
            rec["best"] = tier
    return out


def load_membership(path: str) -> dict:
    if not path:
        return {}
    df = _read_tsv(path, "species-tree membership")
    if df.empty:
        return {}
    return {g.strip(): s.strip() for g, s in zip(df["genome_id"], df["status"])}


def count_defense_systems(df_dir: str, genome: str, system: str = "Prometheus"):
    """(n_defense, n_defense_excl_focal) for one genome's DefenseFinder systems table."""
    f = os.path.join(df_dir, genome, f"{genome}_defense_finder_systems.tsv")
    df = _read_tsv(f, f"DefenseFinder systems for {genome}")
    if df.empty:
        return 0, 0
    if "type" not in df.columns:
        raise ValueError(f"{f} has no 'type' column")
    if "activity" in df.columns:
        df = df[df["activity"].str.strip().str.casefold() != "antidefense"]
    n_total = int(len(df))
    focal = df["type"].str.strip().str.casefold() == system.strip().casefold()
    return n_total, int(n_total - focal.sum())


def ecology_columns(columns) -> list:
    groups = sorted(c for c in columns if c.startswith("group_"))
    niches = sorted(c for c in columns if c.startswith("niche_"))
    return groups + niches


def load_ecology(path: str):
    """species table -> (eco_cols, {genome: {col: 1|0|''}}, n_coerced_cells)."""
    df = _read_tsv(path, "species/ecology table")
    if "genome_id" not in df.columns:
        raise ValueError(f"{path} has no genome_id column")
    cols = ecology_columns(df.columns)
    out, coerced = {}, 0
    for _, row in df.iterrows():
        gid = row["genome_id"].strip()
        if not gid:
            continue
        if gid in out:
            raise ValueError(f"duplicate genome_id {gid} in {path}")
        vals = {}
        for c in cols:
            v, k = _eco_value(row[c])
            vals[c] = v
            coerced += k
        out[gid] = vals
    return cols, out, coerced


def build_master_header(eco_cols, with_pathways: bool = False) -> list:
    return (list(MASTER_BACKBONE) + list(eco_cols) + list(ADDITIVE) +
            (list(MGEN_ADDITIVE) if with_pathways else []))


def load_pathways(path: str) -> dict:
    """methanogenesis_pathways.tsv -> {genome_id: {mgen_* column: '1' | '0' | 'NA'}}.

    present -> 1; absent -> 0; unresolved / not_applicable -> NA. A one-sided flag's
    'not_detected' is NA, never 0: those systems have enzymes outside the KO set."""
    df = _read_tsv(path, "methanogenesis pathways")
    if df.empty:
        return {}
    need = ["genome_id"] + MGEN_PATHWAYS + MGEN_SUBS_2 + MGEN_SUBS_1 + \
           ["subtype_status", "methanogenesis_subtype", "taxonomy_concordance"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"{path} lacks column(s) {miss}")
    two = {"present": "1", "absent": "0"}
    one = {"present": "1"}
    out = {}
    for _, row in df.iterrows():
        gid = str(row["genome_id"]).strip()
        if not gid:
            continue
        rec = {f"mgen_{p}": two.get(str(row[p]).strip(), "NA")
               for p in MGEN_PATHWAYS + MGEN_SUBS_2}
        rec.update({f"mgen_{x}": one.get(str(row[x]).strip(), "NA") for x in MGEN_SUBS_1})
        rec["subtype_status"] = _blank(row["subtype_status"]) or "NA"
        rec["subtype_taxonomy_concordance"] = _blank(row["taxonomy_concordance"]) or "NA"
        rec["methanogenesis_subtype"] = _blank(row["methanogenesis_subtype"]) or "NA"
        out[gid] = rec
    return out


def _check_subset(names, universe, what):
    extra = sorted(set(names) - universe)
    if extra:
        raise ValueError(f"{what} names {len(extra)} genome(s) outside the genome universe "
                         f"(e.g. {extra[:5]}): inputs come from different genome sets")


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True, help="genome_metadata.tsv")
    ap.add_argument("--taxonomy", required=True, help="taxonomy.tsv (user_genome, classification)")
    ap.add_argument("--tree", required=True, help="species_tree.treefile")
    ap.add_argument("--methanogen", required=True, help="methanogen_calls.tsv")
    ap.add_argument("--pathways", default="",
                    help="methanogenesis_pathways.tsv (adds the mgen_* columns and replaces "
                         "methanogenesis_subtype with the composite pathway label)")
    ap.add_argument("--ecology", required=True, help="species_metadata.tsv (group_*/niche_*)")
    ap.add_argument("--envo-map", default="", help="accepted for compatibility; ignored")
    ap.add_argument("--catalog", required=True, help="prometheus_catalog.tsv")
    ap.add_argument("--sensitivity", required=True, help="detection_sensitivity.tsv")
    ap.add_argument("--reps", required=True, help="representatives.txt")
    ap.add_argument("--df-dir", required=True, help="DefenseFinder output root")
    ap.add_argument("--system", default="Prometheus",
                    help="DefenseFinder type of the focal system (exact match)")
    ap.add_argument("--primary-tiers", default="A", help="stats.confidence_tiers_primary")
    ap.add_argument("--sensitivity-tiers", default="A,B", help="sensitivity tier set")
    ap.add_argument("--tree-membership", default="", help="species_tree_membership.tsv")
    ap.add_argument("--out-master", required=True)
    ap.add_argument("--out-ecology", required=True)
    ap.add_argument("--out-missingness", default="")
    args = ap.parse_args(argv)

    for p in (args.out_master, args.out_ecology, args.out_missingness):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    primary, sens_tiers = parse_tiers(args.primary_tiers), parse_tiers(args.sensitivity_tiers)
    if args.envo_map:
        utils.log(f"--envo-map {args.envo_map} ignored: ecology predictors are the curated "
                  "binary group_*/niche_* columns of the species table (no free text to map).")

    eco_cols, ecology, n_coerced = load_ecology(args.ecology)
    master_header = build_master_header(eco_cols, with_pathways=bool(args.pathways))
    ecology_header = ["genome_id"] + list(eco_cols)

    meta = _read_tsv(args.meta, "genome metadata")
    if "genome_id" not in meta.columns:
        raise ValueError(f"{args.meta} has no genome_id column")
    meta["genome_id"] = meta["genome_id"].str.strip()
    meta = meta[meta["genome_id"] != ""]
    dup = meta["genome_id"][meta["genome_id"].duplicated()].unique().tolist()
    if dup:
        raise ValueError(f"duplicated genome_id in {args.meta}: {dup[:5]}")
    universe = set(meta["genome_id"])

    reps = load_reps(args.reps)
    tips = load_tree_tips(args.tree)
    taxonomy = parse_taxonomy(args.taxonomy)
    methanogen = load_methanogen(args.methanogen)
    pathways = load_pathways(args.pathways) if args.pathways else {}
    sensitivity = load_sensitivity(args.sensitivity)
    prometheus = load_prometheus(args.catalog)
    membership = load_membership(args.tree_membership)
    for names, what in ((tips, "the species tree"), (taxonomy, "the taxonomy table"),
                        (methanogen, "methanogen_calls"), (sensitivity, "detection_sensitivity"),
                        (pathways, "methanogenesis_pathways"),
                        (prometheus, "the Prometheus catalog"), (reps, "representatives.txt"),
                        (membership, "species-tree membership")):
        _check_subset(names, universe, what)
    missing_eco = sorted(universe - set(ecology))
    if missing_eco:
        raise ValueError(f"{len(missing_eco)} genome(s) are missing from the species table "
                         f"{args.ecology} (e.g. {missing_eco[:5]})")
    missing_tax = sorted(universe - set(taxonomy))
    if missing_tax:
        raise ValueError(f"{len(missing_tax)} genome(s) are missing from the taxonomy table "
                         f"(e.g. {missing_tax[:5]})")

    master_rows, ecology_rows = [], []
    no_m = {"mcrA": False, "n_mcr": 0, "methanogen": False, "is_anme": False,
            "methanogenesis_subtype": "NA", "mcr_group": "NA"}
    no_p = {c: "NA" for c in MGEN_ADDITIVE}
    if pathways:
        missing_p = sorted(g for g in universe if g not in pathways)
        if missing_p:
            raise ValueError(f"{len(missing_p)} genome(s) are missing from --pathways "
                             f"(e.g. {missing_p[:5]}): it must cover the genome universe")
    for _, mrow in meta.iterrows():
        gid = mrow["genome_id"]
        tax = taxonomy[gid]
        m = methanogen.get(gid, no_m)
        prom = prometheus.get(gid, {"present": False, "n_loci": 0, "best": ""})
        sens = sensitivity.get(gid, {"detect_defensefinder": 0})
        n_def, n_def_excl = count_defense_systems(args.df_dir, gid, args.system)
        eco = ecology[gid]
        row = {
            "genome_id": gid,
            "representative": gid in reps,
            "tree_tip": gid if gid in tips else "",
            **{col: tax.get(col, "") for _, col in GTDB_RANKS},
            "completeness": _blank(mrow.get("completeness", "")),
            "contamination": _blank(mrow.get("contamination", "")),
            "quality_tier": _blank(mrow.get("quality_tier", "")),
            "pass_primary": _blank(mrow.get("pass_primary", "")),
            "genome_size": _blank(mrow.get("genome_size", "")),
            "mcrA": m["mcrA"], "n_mcr": m["n_mcr"], "methanogen": m["methanogen"],
            "is_anme": m["is_anme"],
            "methanogenesis_subtype": (pathways[gid]["methanogenesis_subtype"] if pathways
                                       else m["methanogenesis_subtype"]),
            "prometheus_present": prom["present"], "prometheus_tier": prom["best"],
            "prometheus_n_loci": prom["n_loci"],
            "detect_defensefinder": sens["detect_defensefinder"],
            "n_defense_systems": n_def, "n_defense_systems_excl_prometheus": n_def_excl,
            "prometheus_status_primary": tier_status(prom["best"], primary),
            "prometheus_status_sensitivity": tier_status(prom["best"], sens_tiers),
            "n_contigs": _blank(mrow.get("n_contigs", "")),
            "n50": _blank(mrow.get("n50", "")),
            "gc": _blank(mrow.get("gc", "")),
            "tree_status": ("in_tree" if gid in tips else membership.get(gid, "not_in_tree"))
            if (tips or membership) else "no_species_tree",
            "mcr_group": m["mcr_group"],
        }
        for c in eco_cols:
            row[c] = eco[c]
        if pathways:
            pr = pathways.get(gid, no_p)
            for c in MGEN_ADDITIVE:
                row[c] = pr.get(c, "NA")
        master_rows.append(row)
        ecology_rows.append({"genome_id": gid, **{c: eco[c] for c in eco_cols}})

    if master_rows:
        pd.DataFrame(master_rows, columns=master_header).to_csv(
            args.out_master, sep="\t", index=False)
        pd.DataFrame(ecology_rows, columns=ecology_header).to_csv(
            args.out_ecology, sep="\t", index=False)
    else:
        utils.write_empty(args.out_master, "\t".join(master_header))  # nonfatal: empty universe
        utils.write_empty(args.out_ecology, "\t".join(ecology_header))
    if args.out_missingness:
        lines = ["column\tn_genomes\tn_1\tn_0\tn_NA"]
        for c in eco_cols:
            vals = [r[c] for r in ecology_rows]
            lines.append(f"{c}\t{len(vals)}\t{sum(v == 1 for v in vals)}\t"
                         f"{sum(v == 0 for v in vals)}\t{sum(v == '' for v in vals)}")
        with open(args.out_missingness, "w") as fh:
            fh.write("\n".join(lines) + "\n")

    n = len(master_rows)
    st = {}
    for r in master_rows:
        st[r["prometheus_status_primary"]] = st.get(r["prometheus_status_primary"], 0) + 1
    n_tip = sum(1 for r in master_rows if r["tree_tip"])
    utils.log(f"master_table: {n} genomes; {sum(r['prometheus_present'] for r in master_rows)} "
              f"with a non-rejected Prometheus locus; primary tiers {sorted(primary)}: {st}; "
              f"{sum(r['methanogen'] for r in master_rows)} methanogens; {n_tip}/{n} on the "
              f"species tree; {len(eco_cols)} ecology columns"
              + (f"; {n_coerced} non-binary ecology cell(s) coerced to 1" if n_coerced else ""))
    if n and n_tip < n and tips:
        utils.log(f"WARNING: {n - n_tip} genome(s) are not species-tree tips (see tree_status).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
