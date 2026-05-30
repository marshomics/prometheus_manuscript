"""Aggregate proviral load, CheckV quality, and R-loop data per genome."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

HEADER = (
    "genome_id\tn_proviruses\ttotal_provirus_kb\t"
    "mean_checkv_completeness\tn_rloop_proviruses\tfrac_rloop_proviruses"
)


def read_virus_summary(path: str):
    """Return (n_proviruses, total_provirus_kb) from geNomad virus summary TSV."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, 0.0
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return 0, 0.0
        n = len(df)
        # Look for a length column
        len_col = None
        for c in df.columns:
            if c.lower() in ("length", "seq_length", "virus_length", "len"):
                len_col = c
                break
        total_kb = 0.0
        if len_col is not None:
            try:
                total_kb = pd.to_numeric(df[len_col], errors="coerce").sum() / 1000.0
            except Exception:
                total_kb = 0.0
        return n, total_kb
    except Exception:
        return 0, 0.0


def read_checkv_quality(path: str):
    """Return mean completeness (numeric) from CheckV quality_summary.tsv."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return None
        comp_col = None
        for c in df.columns:
            if c.lower() in ("completeness", "completeness_high", "checkv_completeness"):
                comp_col = c
                break
        if comp_col is None:
            return None
        vals = pd.to_numeric(df[comp_col], errors="coerce").dropna()
        if vals.empty:
            return None
        return float(vals.mean())
    except Exception:
        return None


def read_rloop(path: str, n_proviruses: int):
    """Return (n_rloop_proviruses, frac_rloop_proviruses) from rloop.tsv."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0, 0.0
    try:
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return 0, 0.0
        # Find rloop_prone boolean column
        prone_col = None
        for c in df.columns:
            if "prone" in c.lower():
                prone_col = c
                break
        if prone_col is None:
            return 0, 0.0
        # Coerce to bool-like
        prone_vals = df[prone_col]
        # Handle True/False strings
        if prone_vals.dtype == object:
            prone_vals = prone_vals.map(
                lambda x: str(x).strip().lower() in ("true", "1", "yes")
            )
        else:
            prone_vals = prone_vals.astype(bool)
        n_prone = int(prone_vals.sum())
        frac = (n_prone / n_proviruses) if n_proviruses > 0 else 0.0
        return n_prone, frac
    except Exception:
        return 0, 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--viral-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    try:
        samples = pd.read_csv(args.samples, sep="\t")
    except Exception:
        utils.write_empty(args.out, HEADER)
        return

    if samples.empty or "genome_id" not in samples.columns:
        utils.write_empty(args.out, HEADER)
        return

    rows = []
    for _, srow in samples.iterrows():
        gid = str(srow["genome_id"]).strip()

        virus_summary_path = os.path.join(
            args.viral_dir, gid, f"{gid}_summary", f"{gid}_virus_summary.tsv"
        )
        checkv_path = os.path.join(
            args.viral_dir, gid, "checkv", "quality_summary.tsv"
        )
        rloop_path = os.path.join(args.viral_dir, gid, "rloop.tsv")

        n_proviruses, total_provirus_kb = read_virus_summary(virus_summary_path)
        mean_checkv = read_checkv_quality(checkv_path)
        n_rloop, frac_rloop = read_rloop(rloop_path, n_proviruses)

        rows.append({
            "genome_id": gid,
            "n_proviruses": n_proviruses,
            "total_provirus_kb": round(total_provirus_kb, 3),
            "mean_checkv_completeness": round(mean_checkv, 3) if mean_checkv is not None else "",
            "n_rloop_proviruses": n_rloop,
            "frac_rloop_proviruses": round(frac_rloop, 4),
        })

    if not rows:
        utils.write_empty(args.out, HEADER)
        return

    out_df = pd.DataFrame(rows, columns=HEADER.split("\t"))
    out_df.to_csv(args.out, sep="\t", index=False)
    utils.log(f"Wrote viral summary for {len(rows)} genomes -> {args.out}")


if __name__ == "__main__":
    main()
