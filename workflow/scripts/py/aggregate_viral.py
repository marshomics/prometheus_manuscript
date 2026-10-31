"""Per-genome viral-pressure summary: geNomad 1.8 proviruses, CheckV 1.0.3 quality and
QmRLFS R-loop prediction (rloop_predict.py).

Inputs per genome, under --viral-dir/{genome}/ (the layout the genomad/checkv/rloop rules
write; the bacterial branch reuses it):
  {genome}_summary/{genome}_virus_summary.tsv  geNomad `end-to-end` virus summary. geNomad
      1.8.0 always writes it (header-only when nothing is found): seq_name, length,
      topology, coordinates, n_genes, genetic_code, virus_score, fdr, n_hallmarks,
      marker_enrichment, taxonomy. topology is 'Provirus' for a host-flanked region excised
      from a contig (seq_name '<contig>|provirus_<start>_<end>') and 'DTR' / 'ITR' /
      'No terminal repeats' for a WHOLE contig classified as viral.
  checkv/quality_summary.tsv  CheckV `end_to_end` (contig_id, ..., checkv_quality,
      miuvig_quality, completeness ['NA' when undetermined], ...); header-only when geNomad
      found nothing.
  rloop.tsv  rloop_predict.py (provirus_id, length, n_rlfs, rloop_prone, ..., gc).

What is counted (SCHEMA ecology/viral/viral_summary.tsv; first six columns are the contract)
  n_proviruses            integrated proviruses only (topology == 'Provirus')
  total_provirus_kb       summed provirus length / 1000
  mean_checkv_completeness  mean CheckV completeness of the proviruses with an estimate
                          ('' when none has one)
  n_rloop_proviruses      proviruses with >= 1 predicted RLFS (rloop_prone)
  frac_rloop_proviruses   n_rloop_proviruses / n_proviruses; '' (undefined, NOT 0) for a
                          genome without proviruses, so "no provirus" never masquerades as
                          "proviruses without R-loops" in the association test
  + n_viral_contigs, total_viral_contig_kb (whole contigs called viral: free/extrachromosomal
    viruses, or proviruses split off by assembly fragmentation -- reported separately,
    never mixed into the proviral load), n_viral_sequences (= proviruses + viral contigs),
    n_rloop_viral_sequences, frac_rloop_viral_sequences, n_proviruses_checkv_mq (CheckV
    Complete/High/Medium-quality proviruses), mean_provirus_kb, mean_provirus_gc.
Zeros therefore mean "geNomad ran and found nothing": the rules fail when geNomad or
CheckV fail, and this script fails (never silently zero-fills) when a per-genome file is
missing, malformed, or inconsistent with the others (ids or lengths that do not match).
"""
import argparse
import math
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

COLUMNS = [
    "genome_id", "n_proviruses", "total_provirus_kb", "mean_checkv_completeness",
    "n_rloop_proviruses", "frac_rloop_proviruses",
    "n_viral_contigs", "total_viral_contig_kb", "n_viral_sequences",
    "n_rloop_viral_sequences", "frac_rloop_viral_sequences", "n_proviruses_checkv_mq",
    "mean_provirus_kb", "mean_provirus_gc",
]
HEADER = "\t".join(COLUMNS)

GENOMAD_REQUIRED = ("seq_name", "length", "topology")
GENOMAD_TOPOLOGIES = {"Provirus", "DTR", "ITR", "No terminal repeats"}
CHECKV_REQUIRED = ("contig_id", "checkv_quality", "completeness")
CHECKV_MQ_PLUS = {"Complete", "High-quality", "Medium-quality"}
RLOOP_REQUIRED = ("provirus_id", "rloop_prone")


class ViralInputError(RuntimeError):
    pass


def _read(path: str, required, what: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise ViralInputError(f"{what} missing: {path}")
    if os.path.getsize(path) == 0:
        raise ViralInputError(f"{what} is an empty file (no header): {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None)
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise ViralInputError(f"{what} {path} lacks column(s) {miss}; got {list(df.columns)[:12]}")
    return df


def _ints(series: pd.Series, what: str) -> pd.Series:
    out = pd.to_numeric(series, errors="coerce")
    if out.isna().any() or (out < 0).any():
        bad = series[out.isna() | (out < 0)].tolist()[:3]
        raise ViralInputError(f"{what}: non-integer length(s) {bad}")
    return out.astype(int)


def _bools(series: pd.Series, what: str) -> pd.Series:
    s = series.str.strip().str.lower()
    ok = s.isin(["true", "false"])
    if not ok.all():
        raise ViralInputError(f"{what}: rloop_prone must be True/False, got {series[~ok].tolist()[:3]}")
    return s == "true"


def _same_ids(a, b, what_a, what_b, gid):
    a, b = set(a), set(b)
    if a != b:
        only_a, only_b = sorted(a - b)[:3], sorted(b - a)[:3]
        raise ViralInputError(
            f"{gid}: {what_a} and {what_b} disagree on the viral sequence ids "
            f"(only in {what_a}: {only_a}; only in {what_b}: {only_b}) -- stale or partial outputs")


def _mean(values):
    vals = [v for v in values if v is not None and not math.isnan(v)]
    return (sum(vals) / len(vals)) if vals else None


def summarize_genome(viral_dir: str, gid: str) -> dict:
    gdir = os.path.join(viral_dir, gid)
    vs = _read(os.path.join(gdir, f"{gid}_summary", f"{gid}_virus_summary.tsv"),
               GENOMAD_REQUIRED, "geNomad virus summary")
    cv = _read(os.path.join(gdir, "checkv", "quality_summary.tsv"),
               CHECKV_REQUIRED, "CheckV quality summary")
    rl = _read(os.path.join(gdir, "rloop.tsv"), RLOOP_REQUIRED, "R-loop table")

    if vs["seq_name"].duplicated().any():
        raise ViralInputError(f"{gid}: duplicated seq_name in the geNomad virus summary")
    bad_top = sorted(set(vs["topology"]) - GENOMAD_TOPOLOGIES)
    if bad_top:
        raise ViralInputError(f"{gid}: unexpected geNomad topology value(s) {bad_top} "
                              f"(geNomad 1.8 writes {sorted(GENOMAD_TOPOLOGIES)})")
    vs["length_i"] = _ints(vs["length"], f"{gid} geNomad length")
    vs["is_provirus"] = vs["topology"] == "Provirus"
    named = vs["seq_name"].str.contains("|provirus_", regex=False)
    if (named != vs["is_provirus"]).any():
        utils.log(f"WARNING {gid}: {(named != vs['is_provirus']).sum()} sequence(s) whose "
                  "'|provirus_' name and topology disagree; topology is used")

    for df, col, what in ((cv, "contig_id", "CheckV"), (rl, "provirus_id", "rloop.tsv")):
        if df[col].duplicated().any():
            raise ViralInputError(f"{gid}: duplicated {col} in {what}: "
                                  f"{df.loc[df[col].duplicated(), col].tolist()[:3]}")
    _same_ids(vs["seq_name"], cv["contig_id"], "geNomad", "CheckV", gid)
    _same_ids(vs["seq_name"], rl["provirus_id"], "geNomad", "rloop.tsv", gid)

    cv = cv.set_index("contig_id")
    rl = rl.set_index("provirus_id")
    rl["prone"] = _bools(rl["rloop_prone"], f"{gid} rloop.tsv")
    if "length" in rl.columns and len(rl):
        rl_len = _ints(rl["length"], f"{gid} rloop.tsv length")
        vs_len = vs.set_index("seq_name")["length_i"]
        diff = (rl_len.reindex(vs_len.index) != vs_len)
        if diff.any():
            raise ViralInputError(f"{gid}: sequence lengths differ between geNomad and rloop.tsv "
                                  f"for {list(vs_len.index[diff])[:3]} -- stale outputs")

    prov = vs[vs["is_provirus"]]
    cont = vs[~vs["is_provirus"]]
    n_prov, n_cont = len(prov), len(cont)
    prov_ids = prov["seq_name"].tolist()
    all_ids = vs["seq_name"].tolist()

    comp = pd.to_numeric(cv.loc[prov_ids, "completeness"], errors="coerce") if n_prov else pd.Series(dtype=float)
    mean_comp = _mean(comp.tolist())
    n_mq = int(cv.loc[prov_ids, "checkv_quality"].isin(CHECKV_MQ_PLUS).sum()) if n_prov else 0
    n_rl_prov = int(rl.loc[prov_ids, "prone"].sum()) if n_prov else 0
    n_rl_all = int(rl.loc[all_ids, "prone"].sum()) if all_ids else 0
    gc = None
    if n_prov and "gc" in rl.columns:
        gc = _mean(pd.to_numeric(rl.loc[prov_ids, "gc"], errors="coerce").tolist())
    tot_prov_kb = prov["length_i"].sum() / 1000.0
    return {
        "genome_id": gid,
        "n_proviruses": n_prov,
        "total_provirus_kb": round(float(tot_prov_kb), 3),
        "mean_checkv_completeness": "" if mean_comp is None else round(mean_comp, 3),
        "n_rloop_proviruses": n_rl_prov,
        "frac_rloop_proviruses": round(n_rl_prov / n_prov, 4) if n_prov else "",
        "n_viral_contigs": n_cont,
        "total_viral_contig_kb": round(float(cont["length_i"].sum() / 1000.0), 3),
        "n_viral_sequences": n_prov + n_cont,
        "n_rloop_viral_sequences": n_rl_all,
        "frac_rloop_viral_sequences": round(n_rl_all / (n_prov + n_cont), 4) if (n_prov + n_cont) else "",
        "n_proviruses_checkv_mq": n_mq,
        "mean_provirus_kb": round(float(tot_prov_kb) / n_prov, 3) if n_prov else "",
        "mean_provirus_gc": "" if gc is None else round(gc, 5),
    }


def read_genome_ids(path: str):
    """Genome ids from the samples table: genome_id (archaeal species table), else
    user_genome, else the first column (the bacterial qc fallback of common.smk)."""
    if not os.path.exists(path):
        raise ViralInputError(f"samples table missing: {path}")
    if os.path.getsize(path) == 0:
        return []
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")
    if df.shape[1] == 0:
        return []
    col = "genome_id" if "genome_id" in df.columns else (
        "user_genome" if "user_genome" in df.columns else df.columns[0])
    ids = [g.strip() for g in df[col].tolist() if g.strip()]
    dup = sorted({g for g in ids if ids.count(g) > 1})
    if dup:
        raise ViralInputError(f"duplicated genome id(s) in {path}: {dup[:5]}")
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--viral-dir", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    try:
        ids = read_genome_ids(args.samples)
        rows = [summarize_genome(args.viral_dir, gid) for gid in ids]
    except ViralInputError as exc:
        sys.exit(f"[aggregate_viral] ERROR: {exc}")
    if not rows:
        utils.write_empty(args.out, HEADER)
        utils.log(f"no genomes in {args.samples}; header-only {args.out}")
        return
    pd.DataFrame(rows, columns=COLUMNS).to_csv(args.out, sep="\t", index=False)
    n_with = sum(1 for r in rows if r["n_proviruses"] > 0)
    utils.log(f"viral summary for {len(rows)} genome(s) ({n_with} with >= 1 provirus, "
              f"{sum(r['n_proviruses'] for r in rows)} proviruses, "
              f"{sum(r['n_viral_contigs'] for r in rows)} whole viral contigs) -> {args.out}")


if __name__ == "__main__":
    main()
