"""Join GTDB-Tk taxonomy + CheckM2 completeness/contamination + BBMap assembly stats onto
the Phase-1 input metadata (one row per genome) and add the quality score

    quality_score = completeness - 5 * contamination      (rounded to 1e-6)

Inputs are the EXACT per-batch files of the current DAG (--gtdbtk/--checkm2/--bbmap
lists; a glob would also pick up stale batch directories of an earlier, larger run).
Every input genome is kept (left join). Nothing is silently lost: each genome carries
  gtdbtk_status  classified | unclassified (GTDB-Tk 'Unclassified*': failed gene calling,
                 too few markers) | missing (absent from every GTDB-Tk summary)
  checkm2_status ok | missing (CheckM2 1.0.2 skips bins whose protein file is empty)
  bbmap_status   ok | missing
and genomes with missing values simply fail the Phase-1 quality filter, which reports the
reason per genome (filter_collapse.py --report).

Hard errors (exit 1): duplicated ids, a tool row for a genome that is not in the input
metadata (= the file-stem/genome_id contract is broken), a manifest that does not match
the metadata, a configured BBMap column that is absent, or a non-numeric BBMap/CheckM2
value (e.g. a 'KB'-formatted N50 from a non-machine statswrapper format).

Column mapping (config preprocess.bbmap): gc = gc_avg (FRACTION 0-1), genome_size =
contig_bp (assembled bp excluding scaffold gaps), n_contigs = n_contigs, n50 = ctg_L50
(BBMap's ctg_L50 is the contig N50 LENGTH in bp; its ctg_N50 is the L50 COUNT).
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
RANK_COLS = [c for _, c in RANKS]
# GTDB-Tk 2.4.0 summary columns (gtdbtk/files/classify_summary.py) carried through.
GTDBTK_EXTRA = [("marker_set", "gtdbtk_marker_set"),
                ("classification_method", "gtdbtk_classification_method"),
                ("closest_genome_reference", "gtdbtk_closest_genome_reference"),
                ("closest_genome_ani", "gtdbtk_closest_genome_ani"),
                ("msa_percent", "gtdbtk_msa_percent"),
                ("red_value", "gtdbtk_red_value"),
                ("translation_table", "gtdbtk_translation_table"),
                ("note", "gtdbtk_note"),
                ("warnings", "gtdbtk_warnings")]
# CheckM2 1.0.2 quality_report.tsv columns (checkm2/predictQuality.py) carried through.
CHECKM2_EXTRA = [("Completeness_Model_Used", "checkm2_model"),
                 ("Translation_Table_Used", "checkm2_translation_table"),
                 ("Contig_N50", "checkm2_contig_n50"),
                 ("Genome_Size", "checkm2_genome_size"),
                 ("Additional_Notes", "checkm2_notes")]
NA_TOKENS = {"", "n/a", "na", "nan", "none", "null"}


class MergeError(Exception):
    pass


def split_classification(s) -> dict:
    """GTDB-Tk classification string -> rank columns + classification + gtdbtk_status."""
    out = {col: "" for col in RANK_COLS}
    c = "" if s is None else str(s).strip()
    if c.lower() in NA_TOKENS:
        c = ""
    out["classification"] = c
    if c.lower().startswith("unclassified") or not c:
        out["gtdbtk_status"] = "unclassified"
        return out
    out["gtdbtk_status"] = "classified"
    for field in c.split(";"):
        field = field.strip()
        for pre, col in RANKS:
            if field.startswith(pre):
                out[col] = field[len(pre):].strip()
    if not out["gtdb_domain"]:
        out["gtdbtk_status"] = "unclassified"
    return out


def _pick(df, name):
    """Case-insensitive column lookup (None when absent)."""
    low = {str(c).strip().lower(): c for c in df.columns}
    return low.get(name.lower())


def _read_tables(paths, label):
    """Read and concatenate TSVs (all values as text, 'N/A' kept verbatim)."""
    import pandas as pd
    frames = []
    for p in paths:
        try:
            d = pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            raise MergeError(f"{label} file {p} is empty (no header)")
        d.columns = [str(c).strip() for c in d.columns]
        d["__file"] = p
        frames.append(d)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True)


def _check_ids(ids, known, label):
    from collections import Counter
    ids = [str(i).strip() for i in ids]
    dup = sorted(i for i, n in Counter(ids).items() if n > 1)
    if dup:
        raise MergeError(f"{label}: {len(dup)} genome id(s) reported more than once: {dup[:5]}")
    unknown = sorted(set(ids) - known)
    if unknown:
        raise MergeError(f"{label}: {len(unknown)} id(s) are not genomes of the input metadata "
                         f"(file stem != genome_id?): {unknown[:5]}")


def _numeric(series, label, col):
    """Text -> float; a non-empty, non-NA token that does not parse is an error."""
    import numpy as np
    import pandas as pd
    s = series.fillna("").astype(str).str.strip()
    blank = s.str.lower().isin(NA_TOKENS)
    num = pd.to_numeric(s.where(~blank, np.nan), errors="coerce")
    bad = s[(~blank) & num.isna()]
    if len(bad):
        raise MergeError(f"{label} column '{col}' has {len(bad)} non-numeric value(s), e.g. "
                         f"{bad.unique().tolist()[:3]} (BBMap must run with format=6, which "
                         "prints plain integers)")
    return num


def read_metadata(path, idc):
    """Input metadata parsed like common.smk (comment='#', blank lines skipped)."""
    import pandas as pd
    with open(path) as fh:            # a mid-line '#' would silently truncate the row
        for n, line in enumerate(fh, 1):
            s = line.strip()
            if s and not s.startswith("#") and "#" in line:
                raise MergeError(f"{path}:{n}: '#' inside a data row (read with comment='#')")
    meta = pd.read_csv(path, sep="\t", dtype=str, comment="#", skip_blank_lines=True)
    if idc not in meta.columns:
        raise MergeError(f"metadata {path} has no '{idc}' column")
    meta[idc] = meta[idc].fillna("").astype(str).str.strip()
    meta = meta[meta[idc] != ""].reset_index(drop=True)
    dup = sorted(meta.loc[meta[idc].duplicated(), idc].unique().tolist())
    if dup:
        raise MergeError(f"metadata {path}: duplicated genome id(s) {dup[:5]}")
    return meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--gtdbtk", nargs="*", default=None, help="per-batch gtdbtk_summary.tsv files")
    ap.add_argument("--checkm2", nargs="*", default=None, help="per-batch checkm2_quality.tsv files")
    ap.add_argument("--bbmap", nargs="*", default=None, help="per-batch bbmap_stats.tsv files")
    ap.add_argument("--manifest", default="", help="batch_genomes.py manifest.tsv (adds batch, source_path)")
    ap.add_argument("--batch-dir", default="",
                    help="legacy: glob batch_*/ when the explicit file lists are not given")
    ap.add_argument("--id-col", default="genome_id")
    ap.add_argument("--gc-col", default="gc_avg")
    ap.add_argument("--length-col", default="contig_bp")
    ap.add_argument("--contigs-col", default="n_contigs")
    ap.add_argument("--n50-col", default="ctg_L50")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    try:
        _run(args)
    except MergeError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)


def _files(explicit, batch_dir, name):
    if explicit is not None:
        return list(explicit)
    if batch_dir:
        return sorted(glob.glob(os.path.join(batch_dir, "batch_*", name)))
    return []


def _run(args) -> None:
    import numpy as np
    import pandas as pd

    idc = args.id_col
    meta = read_metadata(args.metadata, idc)
    known = set(meta[idc])

    # ---- batch manifest ---------------------------------------------------- #
    if args.manifest:
        man = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
        if set(man["genome_id"]) != known or len(man) != len(known):
            raise MergeError(f"manifest {args.manifest} ({len(man)} rows) does not match the "
                             f"{len(known)} metadata genomes; re-run batch_genomes")
        meta = meta.merge(man[["genome_id", "batch", "source_path"]].rename(columns={"genome_id": idc}),
                          on=idc, how="left")

    # ---- GTDB-Tk ----------------------------------------------------------- #
    gt = _read_tables(_files(args.gtdbtk, args.batch_dir, "gtdbtk_summary.tsv"), "GTDB-Tk")
    tax_rows = {}
    if gt is not None and len(gt):
        ug, cl = _pick(gt, "user_genome"), _pick(gt, "classification")
        if ug is None or cl is None:
            raise MergeError(f"GTDB-Tk summary lacks user_genome/classification: {list(gt.columns)}")
        _check_ids(gt[ug].tolist(), known, "GTDB-Tk")
        extra = [(src, dst, _pick(gt, src)) for src, dst in GTDBTK_EXTRA]
        for _, r in gt.iterrows():
            rec = split_classification(r[cl])
            for _src, dst, col in extra:
                v = str(r[col]).strip() if col else ""
                rec[dst] = "" if v.lower() in NA_TOKENS else v
            tax_rows[str(r[ug]).strip()] = rec
    tax_cols = ["gtdbtk_status", "classification"] + RANK_COLS + [d for _, d in GTDBTK_EXTRA]
    tax = pd.DataFrame([{idc: g, **rec} for g, rec in tax_rows.items()], columns=[idc] + tax_cols)

    # ---- CheckM2 ----------------------------------------------------------- #
    cm = _read_tables(_files(args.checkm2, args.batch_dir, "checkm2_quality.tsv"), "CheckM2")
    cm_cols = ["completeness", "contamination"] + [d for _, d in CHECKM2_EXTRA]
    if cm is not None and len(cm):
        nm, comp, cont = _pick(cm, "Name"), _pick(cm, "Completeness"), _pick(cm, "Contamination")
        if None in (nm, comp, cont):
            raise MergeError(f"CheckM2 report lacks Name/Completeness/Contamination: {list(cm.columns)}")
        _check_ids(cm[nm].tolist(), known, "CheckM2")
        cmo = pd.DataFrame({idc: cm[nm].astype(str).str.strip(),
                            "completeness": _numeric(cm[comp], "CheckM2", comp),
                            "contamination": _numeric(cm[cont], "CheckM2", cont)})
        for src, dst in CHECKM2_EXTRA:
            col = _pick(cm, src)
            cmo[dst] = cm[col].astype(str).str.strip() if col else ""
    else:
        cmo = pd.DataFrame(columns=[idc] + cm_cols)

    # ---- BBMap ------------------------------------------------------------- #
    bb = _read_tables(_files(args.bbmap, args.batch_dir, "bbmap_stats.tsv"), "BBMap")
    bb_map = [(args.gc_col, "gc"), (args.length_col, "genome_size"),
              (args.contigs_col, "n_contigs"), (args.n50_col, "n50")]
    if bb is not None and len(bb):
        gid = _pick(bb, "genome_id")
        if gid is None:
            raise MergeError(f"BBMap table lacks genome_id: {list(bb.columns)}")
        _check_ids(bb[gid].tolist(), known, "BBMap")
        bbo = pd.DataFrame({idc: bb[gid].astype(str).str.strip()})
        for src, dst in bb_map:
            col = _pick(bb, src)
            if col is None:
                raise MergeError(f"BBMap column '{src}' (config preprocess.bbmap) not in "
                                 f"statswrapper output {[c for c in bb.columns if c != '__file']}")
            bbo[dst] = _numeric(bb[col], "BBMap", src)
    else:
        bbo = pd.DataFrame(columns=[idc, "gc", "genome_size", "n_contigs", "n50"])

    # ---- merge (left: keep every input genome) ------------------------------ #
    out = meta.merge(tax, on=idc, how="left").merge(cmo, on=idc, how="left").merge(bbo, on=idc, how="left")
    out["gtdbtk_status"] = out["gtdbtk_status"].fillna("missing")
    for c in tax_cols[1:]:
        out[c] = out[c].fillna("")
    out["checkm2_status"] = np.where(out[idc].isin(set(cmo[idc])), "ok", "missing")
    out["bbmap_status"] = np.where(out[idc].isin(set(bbo[idc])), "ok", "missing")
    for c in ("completeness", "contamination", "gc", "genome_size", "n_contigs", "n50"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    out["quality_score"] = (out["completeness"] - 5.0 * out["contamination"]).round(6)

    # ---- consistency checks on the BBMap numbers ---------------------------- #
    both = out["n50"].notna() & out["genome_size"].notna()
    over = out[both & (out["n50"] > out["genome_size"])]
    if len(over):
        raise MergeError(f"{len(over)} genome(s) have n50 > genome_size (e.g. "
                         f"{over[idc].tolist()[:3]}); check config preprocess.bbmap column mapping")
    if "checkm2_contig_n50" in out.columns:
        ck = pd.to_numeric(out["checkm2_contig_n50"], errors="coerce")
        m = ck.notna() & out["n50"].notna() & (ck > 0)
        if m.any():
            ratio = (out.loc[m, "n50"] / ck[m])
            far = int(((ratio < 0.5) | (ratio > 2.0)).sum())
            utils.log(f"BBMap {args.n50_col} vs CheckM2 Contig_N50 (scaffold-level): median ratio "
                      f"{ratio.median():.3f}; {far}/{int(m.sum())} genome(s) differ > 2-fold.")
            if far > 0.5 * int(m.sum()):
                utils.log("WARNING: most genomes disagree > 2-fold; check that n50_col is the "
                          "N50 LENGTH column (BBMap ctg_L50), not the count (ctg_N50).")

    lead = list(meta.columns)
    rest = [c for c in (tax_cols + ["completeness", "contamination", "quality_score", "checkm2_status"]
                        + [d for _, d in CHECKM2_EXTRA] + ["gc", "genome_size", "n_contigs", "n50",
                                                           "bbmap_status"]) if c not in lead]
    out = out[lead + rest]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    vc = out["gtdbtk_status"].value_counts().to_dict()
    n_cm_missing = int((out["checkm2_status"] == "missing").sum())
    n_bb_missing = int((out["bbmap_status"] == "missing").sum())
    utils.log(f"annotated {len(out)} genomes: GTDB-Tk {vc}; CheckM2 missing {n_cm_missing}; "
              f"BBMap missing {n_bb_missing}.")
    for label, n in (("GTDB-Tk", int(vc.get("missing", 0))), ("CheckM2", n_cm_missing),
                     ("BBMap", n_bb_missing)):
        if n:
            utils.log(f"WARNING: {n} genome(s) have no {label} result; they fail the quality "
                      "filter and are listed with the reason in filter_report.tsv.")


if __name__ == "__main__":
    main()
