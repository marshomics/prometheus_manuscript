"""Genome x defense-system presence/absence matrix from DefenseFinder systems tables.

One row per master genome (every genome, in master order); one 0/1 column per
DefenseFinder system `type` (the model family, e.g. RM, CBASS, Gabija, Prometheus --
DefenseFinder v2 derives `type` from the model's family directory and `subtype` from
the model name). Family level is the unit of the co-occurrence analysis: subtypes of
one family are near-exclusive variants and would multiply tests. A `Prometheus`
column is always present. Rows whose DefenseFinder `activity` is "Antidefense"
(AntiDefenseFinder models, `defense-finder run -a`) are not defense systems and are
excluded. The Prometheus column records raw DefenseFinder calls; the curated
tier-based presence lives in the master table and is what cooccurrence.R uses.

Per-genome input: {df_dir}/{g}/{g}_defense_finder_systems.tsv. A header-only table
means "no systems" (all zeros). With --require-all a missing table is an error (the
archaeal rule: every genome has a DefenseFinder output); without it a missing table
counts as zero systems (the bacterial branch writes per-genome tables only for genomes
that have systems) and the number of such genomes is logged.

Empty master (zero rows) -> header-only matrix. Unreadable/empty (0-byte) tables or a
table without a `type` column -> exit 1.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

FOCAL = "Prometheus"


def system_types(path: str):
    """Return (set of defense types, n_antidefense_rows) for one systems table."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    if "type" not in df.columns:
        sys.exit(f"[build_defense_matrix] {path}: no 'type' column (not a DefenseFinder systems table)")
    n_anti = 0
    if "activity" in df.columns:
        anti = df["activity"].astype(str).str.strip().str.lower() == "antidefense"
        n_anti = int(anti.sum())
        df = df[~anti]
    types = {t.strip() for t in df["type"].astype(str) if t.strip()}
    if "genome_id" in types:
        sys.exit(f"[build_defense_matrix] {path}: a system type is literally 'genome_id'")
    return types, n_anti


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--master", required=True)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--require-all", action="store_true",
                    help="a missing per-genome systems table is an error (archaeal rule)")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    master = pd.read_csv(args.master, sep="\t", dtype=str, keep_default_na=False)
    if "genome_id" not in master.columns:
        sys.exit("[build_defense_matrix] master lacks a genome_id column")
    genomes = master["genome_id"].astype(str).str.strip().tolist()
    if not genomes:
        utils.write_empty(args.out, f"genome_id\t{FOCAL}")  # nonfatal: empty input (zero genomes)
        utils.log("master has zero genomes; wrote header-only matrix")
        return
    if any(not g for g in genomes):
        sys.exit("[build_defense_matrix] master has blank genome_id values")
    dup = pd.Series(genomes)[pd.Series(genomes).duplicated()]
    if len(dup):
        sys.exit(f"[build_defense_matrix] duplicated genome_id in master: {dup.iloc[0]}")

    per_genome, all_types = {}, {FOCAL}
    n_missing = n_anti_total = 0
    for gid in genomes:
        path = os.path.join(args.df_dir, gid, f"{gid}_defense_finder_systems.tsv")
        if not os.path.exists(path):
            if args.require_all:
                sys.exit(f"[build_defense_matrix] missing DefenseFinder systems table: {path}")
            n_missing += 1
            per_genome[gid] = set()
            continue
        types, n_anti = system_types(path)
        n_anti_total += n_anti
        per_genome[gid] = types
        all_types |= types

    cols = sorted(all_types)
    out = pd.DataFrame([[gid] + [int(t in per_genome[gid]) for t in cols] for gid in genomes],
                       columns=["genome_id"] + cols)
    out.to_csv(args.out, sep="\t", index=False)
    if n_missing:
        utils.log(f"{n_missing} genome(s) without a systems table counted as zero systems")
    if n_anti_total:
        utils.log(f"excluded {n_anti_total} Antidefense system row(s)")
    utils.log(f"defense matrix: {len(genomes)} genomes x {len(cols)} system types "
              f"({int(out[FOCAL].sum())} with a DefenseFinder {FOCAL} call) -> {args.out}")


if __name__ == "__main__":
    main()
