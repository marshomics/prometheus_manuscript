"""Quality-filter the annotated genome table, then collapse to one representative
per GTDB species, unioning niche/group presence-absence across each species' genomes.

Steps (= user spec)
-------------------
1. quality_score = completeness - 5*contamination (already in the annotated table).
2. Keep genomes with completeness >= MIN, contamination <= MAX, n_contigs < MAXC,
   contig-N50 (BBMap ctg_L50) > MINN50, and quality_score > MINQS. (NaN fails.)
3. Group by GTDB species (s__...). Genomes with no species assignment (novel /
   unassigned per GTDB-Tk) are DROPPED by default (--unclassified drop); set
   --unclassified singleton to keep each as its own species instead.
4. Representative = the passing genome of highest quality_score (ties: completeness,
   then genome_id). A species with no passing genome is dropped.
5. NICHE TRACKING: group_*/niche_* 0/1 columns are unioned (max) across the
   species' genomes so quality-filtering and collapsing don't erase where the
   species was observed (--niche-aggregation all | passing).
6. Stage the representative FASTAs as {stage_dir}/{genome_id}.fna and emit ONE
   filtered table (species_metadata.tsv = both samples and metadata for the
   analysis), plus derived helper files (genome_metadata, taxonomy, reps list).

MAG/isolate (genome_type) is deliberately NOT produced or used anywhere; quality is
represented by Completeness / Contamination / n_contigs / N50 instead.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

RANK_COLS = ["gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
             "gtdb_family", "gtdb_genus", "gtdb_species"]
GENOME_META_HEADER = ("genome_id\tsource\tcompleteness\tcontamination\tgenome_size\t"
                      "n_contigs\tn50\tgc\tquality_tier\tpass_primary")


def _empty_outputs(args):
    utils.write_empty(args.out_species, "genome_id")
    utils.write_empty(args.out_reps_list, "genome_id")
    utils.write_empty(args.out_genome_meta, GENOME_META_HEADER)
    utils.write_empty(args.out_taxonomy, "user_genome\tclassification")
    utils.write_empty(args.out_staged, "genome_id")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--genomes-dir", default="")
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
    ap.add_argument("--stage-dir", required=True)
    ap.add_argument("--out-species", required=True)
    ap.add_argument("--out-reps-list", required=True)
    ap.add_argument("--out-genome-meta", required=True)
    ap.add_argument("--out-taxonomy", required=True)
    ap.add_argument("--out-staged", required=True)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    for p in (args.out_species, args.out_reps_list, args.out_genome_meta,
              args.out_taxonomy, args.out_staged):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    os.makedirs(args.stage_dir, exist_ok=True)

    df = pd.read_csv(args.annotated, sep="\t", dtype=str)
    if "genome_id" not in df.columns or df.empty:
        _empty_outputs(args)
        utils.log("annotated table empty; wrote header-only preprocessing outputs.")
        return

    niche_cols = [c for c in df.columns if c.startswith("group_") or c.startswith("niche_")]
    for c in niche_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).clip(0, 1).astype(int)
    for c in ("completeness", "contamination", "n_contigs", "n50", "gc", "genome_size", "quality_score"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        else:
            df[c] = np.nan

    # domain filter: after GTDB-Tk, keep only the requested domain (e.g. Archaea).
    if args.keep_domain:
        if "gtdb_domain" in df.columns:
            before = len(df)
            dom = df["gtdb_domain"].fillna("").astype(str).str.strip().str.lower()
            df = df[dom == args.keep_domain.strip().lower()].copy()
            utils.log(f"domain filter: kept {len(df)}/{before} genomes with domain "
                      f"'{args.keep_domain}' (d__{args.keep_domain}).")
            if df.empty:
                _empty_outputs(args)
                utils.log("no genomes left after the domain filter.")
                return
        else:
            utils.log(f"WARNING: --keep-domain '{args.keep_domain}' set but no gtdb_domain "
                      "column; domain filter skipped.")

    # species key (novel/unassigned -> dropped by default)
    sp = df["gtdb_species"].fillna("").astype(str).str.strip() if "gtdb_species" in df.columns \
        else pd.Series([""] * len(df))
    def key_for(i):
        s = sp.iloc[i]
        if s and s.lower() not in ("nan", "s__"):
            return s
        return None if args.unclassified == "drop" else f"unclassified__{df['genome_id'].iloc[i]}"
    df["species_key"] = [key_for(i) for i in range(len(df))]

    df["passing"] = (
        (df["completeness"] >= args.min_completeness)
        & (df["contamination"] <= args.max_contamination)
        & (df["n_contigs"] < args.max_contigs)
        & (df["n50"] > args.min_n50)
        & (df["quality_score"] > args.min_quality_score)
    ).fillna(False)

    keyed = df[df["species_key"].notna()].copy()
    agg_set = keyed if args.niche_aggregation == "all" else keyed[keyed["passing"]]

    niche_union = (agg_set.groupby("species_key")[niche_cols].max() if niche_cols
                   else pd.DataFrame(index=agg_set["species_key"].unique()))
    n_in_species = agg_set.groupby("species_key").size()
    n_passing = keyed[keyed["passing"]].groupby("species_key").size()

    passing = keyed[keyed["passing"]].copy()
    if passing.empty:
        _empty_outputs(args)
        utils.log("no genome passed the quality filter; nothing to collapse.")
        return
    passing = passing.sort_values(["species_key", "quality_score", "completeness", "genome_id"],
                                  ascending=[True, False, False, True])
    # Take the ACTUAL first row per species (highest quality_score after the sort).
    # NOTE: GroupBy.first() returns the first NON-NULL value per column independently,
    # which would Frankenstein a representative's metadata together from several
    # genomes when the top row has a NaN (e.g. a missing BBMap stat). drop_duplicates
    # on the (already sorted) species_key keeps the whole top row intact instead.
    reps = passing.drop_duplicates(subset="species_key", keep="first").reset_index(drop=True)

    # one filtered table (= samples + metadata): rep identity + taxonomy + quality
    # + niche union + counts.
    id_cols = [c for c in (["genome_id", args.file_col, args.src_col, "classification"] + RANK_COLS)
               if c in reps.columns]
    qual_cols = ["completeness", "contamination", "gc", "genome_size", "n_contigs", "n50", "quality_score"]
    sp_rows = []
    for _, r in reps.iterrows():
        k = r["species_key"]
        row = {c: r.get(c, "") for c in id_cols}
        # path to the staged assembly used by the analysis (Prodigal input)
        row["assembly_path"] = os.path.join(args.stage_dir, f"{str(r['genome_id'])}.fna")
        for c in qual_cols:
            row[c] = r.get(c, "")
        row["n_genomes_in_species"] = int(n_in_species.get(k, 1))
        row["n_passing_in_species"] = int(n_passing.get(k, 1))
        for c in niche_cols:
            row[c] = int(niche_union.loc[k, c]) if (k in niche_union.index and c in niche_union.columns) else int(r.get(c, 0))
        sp_rows.append(row)
    pd.DataFrame(sp_rows).to_csv(args.out_species, sep="\t", index=False)

    # derived helper files + stage genomes
    src_col = args.src_col if args.src_col in reps.columns else None
    gmeta_rows, staged_ids = [], []
    for _, r in reps.iterrows():
        gid = str(r["genome_id"])
        src = str(r[src_col]) if src_col else ""
        fname = str(r[args.file_col]) if args.file_col in reps.columns and not pd.isna(r.get(args.file_col)) \
            else f"{gid}.{args.input_ext}"
        gmeta_rows.append({
            "genome_id": gid, "source": src,
            "completeness": r.get("completeness", ""), "contamination": r.get("contamination", ""),
            "genome_size": r.get("genome_size", ""), "n_contigs": r.get("n_contigs", ""),
            "n50": r.get("n50", ""), "gc": r.get("gc", ""), "quality_tier": "A", "pass_primary": True})
        src_path = os.path.abspath(os.path.join(args.genomes_dir, fname)) if args.genomes_dir \
            else os.path.abspath(fname)
        dst = os.path.join(args.stage_dir, f"{gid}.fna")
        try:
            if os.path.islink(dst) or os.path.exists(dst):
                os.remove(dst)
            os.symlink(src_path, dst)
        except OSError as exc:
            utils.log(f"WARNING: could not stage {dst} -> {src_path} ({exc})")
        staged_ids.append(gid)

    with open(args.out_reps_list, "w") as fh:
        fh.write("\n".join(staged_ids) + "\n")
    pd.DataFrame(gmeta_rows, columns=GENOME_META_HEADER.split("\t")).to_csv(
        args.out_genome_meta, sep="\t", index=False)
    reps_out = reps.copy()
    reps_out["user_genome"] = reps_out["genome_id"]
    reps_out[["user_genome", "classification"]].to_csv(args.out_taxonomy, sep="\t", index=False)
    with open(args.out_staged, "w") as fh:
        fh.write("genome_id\n" + "\n".join(staged_ids) + "\n")

    utils.log(f"{len(df)} input -> {int(df['passing'].sum())} passing -> {len(reps)} species "
              f"representatives ({len(niche_cols)} niche/group columns unioned over "
              f"'{args.niche_aggregation}' genomes; novel/unassigned species: {args.unclassified}).")


if __name__ == "__main__":
    main()
