"""Join GTDB-Tk taxonomy + CheckM2 completeness/contamination + BBMap assembly
stats onto the input metadata (one row per genome), and add the quality score
quality_score = completeness - 5*contamination.

Every input genome is retained (left join); genomes that GTDB/CheckM2/BBMap did not
return are left blank and will simply fail the Stage-1 quality filter. Taxonomy is
split into per-rank columns (gtdb_domain..gtdb_species) plus the full classification.
"""
from __future__ import annotations
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

RANKS = [("d__", "gtdb_domain"), ("p__", "gtdb_phylum"), ("c__", "gtdb_class"),
         ("o__", "gtdb_order"), ("f__", "gtdb_family"), ("g__", "gtdb_genus"),
         ("s__", "gtdb_species")]


def split_classification(s: str) -> dict:
    out = {col: "" for _, col in RANKS}
    out["classification"] = s if isinstance(s, str) and s and s.lower() != "nan" else ""
    if not out["classification"]:
        return out
    for field in out["classification"].split(";"):
        field = field.strip()
        for pre, col in RANKS:
            if field.startswith(pre):
                out[col] = field[len(pre):].strip()
    return out


def _read_concat(pattern, **kw):
    import pandas as pd
    frames = []
    for f in sorted(glob.glob(pattern)):
        try:
            d = pd.read_csv(f, sep="\t", dtype=str, **kw)
            if not d.empty:
                frames.append(d)
        except Exception as exc:
            utils.log(f"WARNING: could not read {f} ({exc})")
    return pd.concat(frames, ignore_index=True) if frames else None


def _pick(df, name):
    """Case-insensitive column lookup."""
    if df is None:
        return None
    low = {c.lower(): c for c in df.columns}
    return low.get(name.lower())


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--batch-dir", required=True)
    ap.add_argument("--id-col", default="genome_id")
    ap.add_argument("--gc-col", default="gc_avg")
    ap.add_argument("--length-col", default="contig_bp")  # matches config preprocess.bbmap.length_col
    ap.add_argument("--contigs-col", default="n_contigs")
    ap.add_argument("--n50-col", default="ctg_L50")  # BBMap ctg_L50 = N50 length (bp); ctg_N50 = L50 count
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import numpy as np
    import pandas as pd

    meta = pd.read_csv(args.metadata, sep="\t", dtype=str)
    idc = args.id_col
    meta[idc] = meta[idc].astype(str)

    # ---- GTDB-Tk taxonomy ----
    gt = _read_concat(os.path.join(args.batch_dir, "batch_*", "gtdbtk_summary.tsv"))
    tax = {}
    if gt is not None:
        ug = _pick(gt, "user_genome") or gt.columns[0]
        cl = _pick(gt, "classification")
        for _, r in gt.iterrows():
            tax[str(r[ug])] = split_classification(str(r[cl]) if cl else "")
    tax_df = pd.DataFrame.from_dict(tax, orient="index")
    if not tax_df.empty:
        tax_df.index.name = idc
        tax_df = tax_df.reset_index()
    else:
        tax_df = pd.DataFrame(columns=[idc, "classification"] + [c for _, c in RANKS])

    # ---- CheckM2 ----
    cm = _read_concat(os.path.join(args.batch_dir, "batch_*", "checkm2_quality.tsv"))
    if cm is not None:
        nm = _pick(cm, "Name") or cm.columns[0]
        comp = _pick(cm, "Completeness"); cont = _pick(cm, "Contamination")
        cm = cm.rename(columns={nm: idc, comp: "completeness", cont: "contamination"})
        cm = cm[[idc, "completeness", "contamination"]].drop_duplicates(subset=[idc])
    else:
        cm = pd.DataFrame(columns=[idc, "completeness", "contamination"])

    # ---- BBMap ----
    bb = _read_concat(os.path.join(args.batch_dir, "batch_*", "bbmap_stats.tsv"))
    if bb is not None:
        gid = _pick(bb, "genome_id") or bb.columns[0]
        ren = {gid: idc}
        for src, dst in [(args.gc_col, "gc"), (args.length_col, "genome_size"),
                         (args.contigs_col, "n_contigs"), (args.n50_col, "n50")]:
            col = _pick(bb, src)
            if col:
                ren[col] = dst
            else:
                utils.log(f"WARNING: BBMap column '{src}' not found; '{dst}' will be blank")
        bb = bb.rename(columns=ren)
        keep = [idc] + [c for c in ("gc", "genome_size", "n_contigs", "n50") if c in bb.columns]
        bb = bb[keep].drop_duplicates(subset=[idc])
    else:
        bb = pd.DataFrame(columns=[idc, "gc", "genome_size", "n_contigs", "n50"])

    # ---- merge (left: keep every input genome) ----
    out = meta.merge(tax_df, on=idc, how="left").merge(cm, on=idc, how="left").merge(bb, on=idc, how="left")

    # numeric + quality score
    for c in ("completeness", "contamination", "gc", "genome_size", "n_contigs", "n50"):
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
        else:
            out[c] = np.nan
    out["quality_score"] = out["completeness"] - 5.0 * out["contamination"]

    # ---- BBMap N50 unit sanity check (defensible; non-destructive) ----
    # format=6 should emit the chosen N50 column (ctg_L50) in BASE PAIRS. If a
    # BBMap build/locale instead reported it in kilobases, every value is ~1000x
    # too small and the downstream `n50 > min_n50` (default 5000 bp) filter would
    # over-reject nearly every genome. We do NOT silently rescale (that would hide
    # a real problem and could corrupt correctly-bp data); we emit a PROMINENT
    # warning so the operator can inspect the BBMap output and the --n50-col choice.
    n50 = out["n50"].dropna() if "n50" in out.columns else np.array([])
    gsize = out["genome_size"].dropna() if "genome_size" in out.columns else np.array([])
    if len(n50) and len(gsize):
        max_n50 = float(np.nanmax(n50))
        med_size = float(np.nanmedian(gsize))
        if max_n50 < 1000 and med_size > 100000:
            utils.log(
                "WARNING: ***POSSIBLE BBMap UNIT MISMATCH*** the chosen N50 column "
                f"('{args.n50_col}' -> n50) has max={max_n50:g} while median genome_size="
                f"{med_size:g} bp. statswrapper format=6 should report N50 in BASE PAIRS; "
                "these values look like KILOBASES. The >5000 bp N50 quality filter will "
                "OVER-REJECT. Inspect the raw BBMap output and the --n50-col mapping; "
                "values are NOT auto-rescaled.")

    out.to_csv(args.out, sep="\t", index=False)
    n_tax = out["classification"].fillna("").astype(str).str.len().gt(0).sum() if "classification" in out else 0
    n_cm = out["completeness"].notna().sum()
    utils.log(f"annotated {len(out)} genomes; {n_tax} classified, {n_cm} with CheckM2.")


if __name__ == "__main__":
    main()
