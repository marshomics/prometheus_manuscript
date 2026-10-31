"""Quality-filter the annotated genome table, then collapse it to one representative per
GTDB species, unioning niche/group presence-absence across each species' genomes.

Steps
-----
1. Domain filter: keep genomes whose GTDB domain equals --keep-domain (case-insensitive;
   'Archaea' and 'd__Archaea' are equivalent). '' disables the filter. Genomes GTDB-Tk
   could not classify have no domain and are excluded here (reported, never silent).
2. quality_score = completeness - 5 * contamination (recomputed, rounded to 1e-6, so a
   CheckM2 value with two decimals is compared exactly against the thresholds).
3. passing = completeness >= MIN and contamination <= MAX and n_contigs < MAXC and
   n50 > MINN50 and quality_score > MINQS (config preprocess.quality). A missing value
   fails the corresponding test.
4. Species key = gtdb_species. Genomes without a species (novel GTDB species, s__ empty)
   are DROPPED (--unclassified drop, default) or kept as one-genome species
   (--unclassified singleton).
5. Representative = the passing genome with the highest quality_score; ties are broken by
   higher completeness, lower contamination, higher n50, fewer contigs, then genome_id
   (ascending), so the choice is fully deterministic. A species with no passing genome is
   dropped.
6. Niche union: group_*/niche_* columns (0/1; blank/NA = unknown; anything else is an
   error) are OR-ed over every genome of the species (--niche-aggregation all, default:
   quality filtering cannot erase where a species was observed) or over its passing
   genomes only (passing). Unknown is not absence (SCHEMA Contract v2; design §3.4): the
   species value is 1 if any aggregated genome is 1, 0 if none is 1 and at least one is
   0, and blank (NA) if every aggregated genome is unknown, so build_master_table's
   ecology_missingness.tsv reports real missingness. n_genomes_in_species / n_passing_in_species always count all / passing
   genomes of the species, independent of the aggregation mode.
7. Outputs: species_metadata.tsv (the single analysis table: samples + metadata),
   genome_metadata.tsv, taxonomy.tsv (user_genome, classification), representatives.txt
   (one genome_id per line, no header), the .staged marker, the staged representative
   FASTAs ({stage_dir}/{genome_id}.{stage_ext} symlinks to the input FASTAs), and with
   --report a per-genome audit table (status + reasons for every input genome).

MAG/isolate (genome_type) is deliberately NOT produced or used anywhere; quality is
represented by completeness / contamination / n_contigs / n50 instead.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

RANK_COLS = ["gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
             "gtdb_family", "gtdb_genus", "gtdb_species"]
GENOME_META_COLS = ["genome_id", "source", "completeness", "contamination", "genome_size",
                    "n_contigs", "n50", "gc", "quality_tier", "pass_primary"]
GENOME_META_HEADER = "\t".join(GENOME_META_COLS)
QUAL_COLS = ["completeness", "contamination", "gc", "genome_size", "n_contigs", "n50", "quality_score"]
INT_COLS = ["genome_size", "n_contigs", "n50"]
REPORT_COLS = ["genome_id", "gtdbtk_status", "gtdb_domain", "gtdb_species", "species_key",
               "completeness", "contamination", "quality_score", "n_contigs", "n50",
               "passing", "status", "reasons", "species_representative"]
LEDGER = ".staged_links"        # links this script created (never touches other files)
TRUE_TOK = {"1", "1.0", "true", "t", "yes", "y"}
FALSE_TOK = {"0", "0.0", "false", "f", "no", "n"}
NA_TOK = {"", "na", "nan", "n/a", "none", "null"}


class FilterError(Exception):
    pass


def to01(series, col):
    """0/1 presence column -> nullable Int64 (1 / 0 / <NA>); blank/NA tokens stay <NA>
    (unknown, never imputed 0); any other token is an error."""
    import pandas as pd
    s = series.fillna("").astype(str).str.strip().str.lower()
    ok = s.isin(TRUE_TOK | FALSE_TOK | NA_TOK)
    if not ok.all():
        bad = series[~ok].astype(str).unique().tolist()[:5]
        raise FilterError(f"niche/group column '{col}' must be 0/1 (blank = unknown); found {bad}")
    out = pd.Series(pd.NA, index=series.index, dtype="Int64")
    out[s.isin(TRUE_TOK)] = 1
    out[s.isin(FALSE_TOK)] = 0
    return out, int(s.isin(NA_TOK).sum())


def quality_tier(comp, cont, hi_comp, hi_cont, floor, ceil):
    """MIMAG-style completeness/contamination tier (no rRNA/tRNA criterion):
    A = high (>= hi_comp, <= hi_cont), B = medium (>= floor, <= ceil), C = otherwise."""
    if comp >= hi_comp and cont <= hi_cont:
        return "A"
    if comp >= floor and cont <= ceil:
        return "B"
    return "C"


def stage_link(dst, src):
    if os.path.islink(dst):
        if os.readlink(dst) == src:
            return
        os.remove(dst)
    elif os.path.lexists(dst):
        raise FilterError(f"refusing to overwrite non-symlink file {dst} in the staging directory")
    os.symlink(src, dst)


def _fmt_int(df, cols):
    import pandas as pd
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").round().astype("Int64")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--genomes-dir", default="", help="raw input FASTA dir (used when no source_path column)")
    ap.add_argument("--input-ext", default="fasta", help="extension of the raw input FASTAs")
    ap.add_argument("--file-col", default="filename")
    ap.add_argument("--src-col", default="source")
    ap.add_argument("--min-completeness", type=float, required=True)
    ap.add_argument("--max-contamination", type=float, required=True)
    ap.add_argument("--max-contigs", type=float, required=True)
    ap.add_argument("--min-n50", type=float, required=True)
    ap.add_argument("--min-quality-score", type=float, required=True)
    ap.add_argument("--niche-aggregation", default="all", choices=["all", "passing"])
    ap.add_argument("--unclassified", default="drop", choices=["singleton", "drop"])
    ap.add_argument("--keep-domain", default="Archaea",
                    help="keep only genomes whose GTDB domain == this (e.g. Archaea); '' disables")
    ap.add_argument("--stage-dir", required=True, help="where Phase 2 reads genomes (config genomes_dir)")
    ap.add_argument("--stage-ext", default="fna", help="extension Phase 2 expects (config genome_extension)")
    ap.add_argument("--tier-high-completeness", type=float, default=90.0)
    ap.add_argument("--tier-high-contamination", type=float, default=5.0)
    ap.add_argument("--completeness-floor", type=float, default=50.0)
    ap.add_argument("--contamination-ceiling", type=float, default=10.0)
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-reps-list", required=True)
    ap.add_argument("--out-genome-meta", required=True)
    ap.add_argument("--out-taxonomy", required=True)
    ap.add_argument("--out-staged", required=True)
    ap.add_argument("--report", default="", help="optional per-genome filter report TSV")
    args = ap.parse_args()
    try:
        _run(args)
    except FilterError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)


def _run(args) -> None:
    import numpy as np
    import pandas as pd

    outs = [args.out_species, args.out_reps_list, args.out_genome_meta, args.out_taxonomy,
            args.out_staged] + ([args.report] if args.report else [])
    for p in outs:
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    os.makedirs(args.stage_dir, exist_ok=True)

    df = pd.read_csv(args.annotated, sep="\t", dtype=str, keep_default_na=False)
    if "genome_id" not in df.columns:
        raise FilterError(f"{args.annotated} has no genome_id column")
    df["genome_id"] = df["genome_id"].astype(str).str.strip()
    if df["genome_id"].duplicated().any():
        raise FilterError(f"duplicated genome_id in {args.annotated}")

    niche_cols = [c for c in df.columns if c.startswith("group_") or c.startswith("niche_")]
    n_na_total = 0
    for c in niche_cols:
        df[c], n_na = to01(df[c], c)
        n_na_total += n_na
    if n_na_total:
        utils.log(f"{n_na_total} blank/NA niche/group value(s) kept as unknown (NA, not 0).")
    for c in ("completeness", "contamination", "n_contigs", "n50", "gc", "genome_size", "quality_score"):
        df[c] = pd.to_numeric(df[c], errors="coerce") if c in df.columns else np.nan
    for c in ["classification", "gtdbtk_status"] + RANK_COLS:
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].fillna("").astype(str).str.strip()

    qs = (df["completeness"] - 5.0 * df["contamination"]).round(6)
    given = df["quality_score"]
    diff = given.notna() & qs.notna() & ((given - qs).abs() > 1e-6)
    if diff.any():
        utils.log(f"WARNING: {int(diff.sum())} quality_score value(s) in the annotated table differ "
                  "from completeness - 5*contamination; using the recomputed score.")
    df["quality_score"] = qs

    # ---- 1. domain ------------------------------------------------------------
    keep = args.keep_domain.strip()
    if keep.lower().startswith("d__"):
        keep = keep[3:]
    dom = df["gtdb_domain"].str.lower()
    in_domain = (dom == keep.lower()) if keep else pd.Series(True, index=df.index)

    # ---- 3. quality -----------------------------------------------------------
    tests = [
        ("completeness", df["completeness"] >= args.min_completeness, f"completeness<{args.min_completeness:g}"),
        ("contamination", df["contamination"] <= args.max_contamination, f"contamination>{args.max_contamination:g}"),
        ("n_contigs", df["n_contigs"] < args.max_contigs, f"n_contigs>={args.max_contigs:g}"),
        ("n50", df["n50"] > args.min_n50, f"n50<={args.min_n50:g}"),
        ("quality_score", df["quality_score"] > args.min_quality_score, f"quality_score<={args.min_quality_score:g}"),
    ]
    passing = pd.Series(True, index=df.index)
    for _c, ok, _lab in tests:
        passing &= ok.fillna(False).astype(bool)
    df["passing"] = passing

    # ---- 4. species key -------------------------------------------------------
    sp = df["gtdb_species"]
    has_sp = (sp != "") & (sp.str.lower() != "s__")
    if args.unclassified == "singleton":
        key = np.where(has_sp, sp, "unclassified__" + df["genome_id"])
    else:
        key = np.where(has_sp, sp, "")
    df["species_key"] = np.where(in_domain, key, "")
    keyed = df[df["species_key"] != ""].copy()

    # ---- 5. representative ----------------------------------------------------
    cand = keyed[keyed["passing"]].copy()
    cand = cand.sort_values(
        ["species_key", "quality_score", "completeness", "contamination", "n50", "n_contigs", "genome_id"],
        ascending=[True, False, False, True, False, True, True], kind="mergesort", na_position="last")
    reps = cand.drop_duplicates(subset="species_key", keep="first").reset_index(drop=True)
    rep_of = dict(zip(reps["species_key"], reps["genome_id"]))

    # ---- 6. niche union + counts ----------------------------------------------
    agg_set = keyed if args.niche_aggregation == "all" else keyed[keyed["passing"]]
    # max() skips <NA>: any 1 -> 1, else any 0 -> 0, all unknown -> <NA>
    union = agg_set.groupby("species_key")[niche_cols].max() if niche_cols else None
    n_all = keyed.groupby("species_key").size()
    n_pass = keyed[keyed["passing"]].groupby("species_key").size()

    # ---- per-genome report ----------------------------------------------------
    report_rows = []
    for i, r in df.iterrows():
        reasons = []
        if not in_domain.loc[i]:
            st = r["gtdbtk_status"]
            if st == "missing":
                reasons.append("gtdbtk_missing")
            elif st == "unclassified" or not r["gtdb_domain"]:
                reasons.append("gtdbtk_unclassified")
            else:
                reasons.append(f"domain={r['gtdb_domain']}")
        elif r["species_key"] == "":
            reasons.append("no_gtdb_species")
        for c, ok, lab in tests:
            if pd.isna(r[c]):
                reasons.append(f"{c}_missing")
            elif not bool(ok.loc[i]):
                reasons.append(lab)
        k = r["species_key"]
        rep = rep_of.get(k, "") if k else ""
        if not in_domain.loc[i]:
            status = "excluded_domain"
        elif k == "":
            status = "no_species"
        elif not r["passing"]:
            status = "failed_quality"
        elif rep == r["genome_id"]:
            status = "representative"
        else:
            status = "species_member"
        report_rows.append({
            "genome_id": r["genome_id"], "gtdbtk_status": r["gtdbtk_status"],
            "gtdb_domain": r["gtdb_domain"], "gtdb_species": r["gtdb_species"],
            "species_key": k,
            "completeness": r["completeness"], "contamination": r["contamination"],
            "quality_score": r["quality_score"], "n_contigs": r["n_contigs"], "n50": r["n50"],
            "passing": bool(r["passing"]), "status": status, "reasons": ";".join(reasons),
            "species_representative": rep})
    if args.report:
        rep_df = _fmt_int(pd.DataFrame(report_rows, columns=REPORT_COLS), ["n_contigs", "n50"])
        rep_df.to_csv(args.report, sep="\t", index=False)

    # ---- 7a. the single species-level table -----------------------------------
    sp_rows = []
    n_unknown_cells = 0          # species x niche cells with no observation (written blank)
    for _, r in reps.iterrows():
        k = r["species_key"]
        gid = r["genome_id"]
        row = {"genome_id": gid,
               "filename": r[args.file_col] if args.file_col in reps.columns else "",
               "source": r[args.src_col] if args.src_col in reps.columns else "",
               "classification": r["classification"]}
        for extra in (args.file_col, args.src_col):
            if extra in reps.columns and extra not in row:
                row[extra] = r[extra]
        for c in RANK_COLS:
            row[c] = r[c]
        row["assembly_path"] = os.path.join(args.stage_dir, f"{gid}.{args.stage_ext}")
        for c in QUAL_COLS:
            row[c] = r[c]
        row["n_genomes_in_species"] = int(n_all.get(k, 0))
        row["n_passing_in_species"] = int(n_pass.get(k, 0))
        for c in niche_cols:
            v = union.loc[k, c] if (union is not None and k in union.index) else r[c]
            row[c] = "" if pd.isna(v) else int(v)
            n_unknown_cells += int(pd.isna(v))
        sp_rows.append(row)
    sp_cols = (["genome_id", "filename", "source"]
               + [c for c in (args.file_col, args.src_col) if c not in ("filename", "source")
                  and c in df.columns]
               + ["classification"] + RANK_COLS + ["assembly_path"] + QUAL_COLS
               + ["n_genomes_in_species", "n_passing_in_species"] + niche_cols)
    species = _fmt_int(pd.DataFrame(sp_rows, columns=sp_cols), INT_COLS)
    species.to_csv(args.out_species, sep="\t", index=False)

    # ---- 7b. staging ----------------------------------------------------------
    ledger = os.path.join(args.stage_dir, LEDGER)
    prev = set()
    if os.path.exists(ledger):
        with open(ledger) as fh:
            prev = {ln.strip() for ln in fh if ln.strip()}
    staged_ids, problems = [], []
    for _, r in reps.iterrows():
        gid = r["genome_id"]
        src = str(r["source_path"]).strip() if "source_path" in reps.columns else ""
        if not src:
            fname = str(r[args.file_col]).strip() if args.file_col in reps.columns else ""
            fname = fname or f"{gid}.{args.input_ext}"
            src = os.path.join(args.genomes_dir, fname) if args.genomes_dir else fname
        src = os.path.abspath(src)
        if not os.path.isfile(src):
            problems.append(f"{gid}: {src}")
            continue
        stage_link(os.path.join(args.stage_dir, f"{gid}.{args.stage_ext}"), src)
        staged_ids.append(gid)
    if problems:
        raise FilterError(f"{len(problems)} representative FASTA(s) not found: {problems[:5]}")
    removed = 0
    for gid in sorted(prev - set(staged_ids)):
        p = os.path.join(args.stage_dir, f"{gid}.{args.stage_ext}")
        if os.path.islink(p):
            os.remove(p)
            removed += 1
    with open(ledger, "w") as fh:
        fh.write("".join(f"{g}\n" for g in staged_ids))
    if removed:
        utils.log(f"removed {removed} staged link(s) of genomes that are no longer representatives")

    # ---- 7c. derived helper files ---------------------------------------------
    with open(args.out_reps_list, "w") as fh:
        fh.write("".join(f"{g}\n" for g in staged_ids))
    gm = []
    for _, r in reps.iterrows():
        comp, cont = float(r["completeness"]), float(r["contamination"])
        gm.append({"genome_id": r["genome_id"],
                   "source": r[args.src_col] if args.src_col in reps.columns else "",
                   "completeness": comp, "contamination": cont,
                   "genome_size": r["genome_size"], "n_contigs": r["n_contigs"], "n50": r["n50"],
                   "gc": r["gc"],
                   "quality_tier": quality_tier(comp, cont, args.tier_high_completeness,
                                                args.tier_high_contamination,
                                                args.completeness_floor, args.contamination_ceiling),
                   "pass_primary": bool(comp >= args.completeness_floor
                                        and cont <= args.contamination_ceiling)})
    _fmt_int(pd.DataFrame(gm, columns=GENOME_META_COLS), INT_COLS).to_csv(
        args.out_genome_meta, sep="\t", index=False)
    pd.DataFrame({"user_genome": reps["genome_id"], "classification": reps["classification"]},
                 columns=["user_genome", "classification"]).to_csv(args.out_taxonomy, sep="\t", index=False)
    with open(args.out_staged, "w") as fh:
        fh.write("genome_id\n" + "".join(f"{g}\n" for g in staged_ids))

    st = pd.Series([r["status"] for r in report_rows]).value_counts().to_dict() if report_rows else {}
    utils.log(f"{len(df)} input -> {int(in_domain.sum())} in domain '{keep or '(any)'}' -> "
              f"{int(keyed['passing'].sum())} passing with a species key -> {len(reps)} species "
              f"representatives; status counts {st}; {len(niche_cols)} niche/group column(s) "
              f"unioned over '{args.niche_aggregation}' genomes ({n_unknown_cells} species x column "
              f"cell(s) unknown, written blank); unclassified species: {args.unclassified}.")


if __name__ == "__main__":
    main()
