"""Ingest pre-computed bacterial DefenseFinder output and extract Prometheus candidates.

The bacterial side is provided as COMBINED DefenseFinder tables (one systems TSV and
one genes TSV covering all genomes) plus per-genome {genome_id}.faa proteomes. This
script:

  1. splits the combined SYSTEMS table into the per-genome layout the rest of the
     pipeline already understands -- {out_defense_dir}/{g}/{g}_defense_finder_systems.tsv
     -- so `validate_prometheus.py` (defense-context) and `build_defense_matrix.py`
     are reused UNCHANGED on bacteria;
  2. selects Prometheus systems (type OR subtype matches the configured system name,
     case-insensitive substring -- handles 'Prometheus', 'proA', model-name variants);
  3. pulls the member proteins of those systems from the per-genome faa and writes the
     candidate FASTA + sources, with the SAME combined id scheme as archaea
     (genome__origProteinId) so every downstream script and the cross-domain tree line
     up by id.

Detection is taken at face value here (DefenseFinder already applied its mandatory-gene
architecture rule); the architecture/placement/tiering curation is done by the reused
`validate_prometheus.py`, exactly as for archaea.
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SOURCES_HEADER = "protein_id\tgenome_id\tdetect_defensefinder\tdf_score\tdf_subtype"
RAW_HEADER = "protein_id\tgenome_id\tmethod\tscore\tsubtype"
SYS_FALLBACK_COLS = ["genome_id", "genome", "user_genome", "sample", "assembly", "replicon"]


def _pick(df, *candidates):
    low = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand and cand.lower() in low:
            return low[cand.lower()]
    return None


def split_systems_by_genome(systems_df, gcol, out_defense_dir):
    """Write {out_defense_dir}/{g}/{g}_defense_finder_systems.tsv per genome.

    Returns the set of genome_ids written. The per-genome files carry the original
    DefenseFinder systems columns so the reused validate/defense-matrix code is happy.
    """
    import pandas as pd  # noqa: F401
    genomes = set()
    for gid, sub in systems_df.groupby(gcol):
        gid = str(gid).strip()
        if not gid or gid.lower() == "nan":
            continue
        gdir = os.path.join(out_defense_dir, gid)
        os.makedirs(gdir, exist_ok=True)
        sub.drop(columns=[gcol]).to_csv(
            os.path.join(gdir, f"{gid}_defense_finder_systems.tsv"), sep="\t", index=False)
        genomes.add(gid)
    return genomes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--df-systems", required=True)
    ap.add_argument("--df-genes", required=True)
    ap.add_argument("--faa-dir", required=True)
    ap.add_argument("--cds-dir", default="", help="optional per-genome {g}.fna CDS dir; "
                    "when set, the matching Prometheus CDS are written to --out-cds")
    ap.add_argument("--systems-genome-col", default="genome_id")
    ap.add_argument("--genes-genome-col", default="genome_id")
    ap.add_argument("--system", default="Prometheus")
    ap.add_argument("--out-defense-dir", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-sources", required=True)
    ap.add_argument("--out-raw", required=True)
    ap.add_argument("--out-cds", default="", help="combined Prometheus CDS (genome__hit_id)")
    args = ap.parse_args()

    import pandas as pd
    for p in (args.out_candidates, args.out_sources, args.out_raw):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    os.makedirs(args.out_defense_dir, exist_ok=True)

    # ---- read combined tables -------------------------------------------- #
    try:
        sysdf = pd.read_csv(args.df_systems, sep="\t", dtype=str)
        sysdf.columns = [c.strip() for c in sysdf.columns]
    except Exception as exc:
        utils.log(f"FATAL: could not read --df-systems ({exc})")
        sys.exit(1)
    try:
        gdf = pd.read_csv(args.df_genes, sep="\t", dtype=str)
        gdf.columns = [c.strip() for c in gdf.columns]
    except Exception as exc:
        utils.log(f"FATAL: could not read --df-genes ({exc})")
        sys.exit(1)

    sgcol = _pick(sysdf, args.systems_genome_col, *SYS_FALLBACK_COLS)
    ggcol = _pick(gdf, args.genes_genome_col, *SYS_FALLBACK_COLS)
    if sgcol is None:
        utils.log(f"FATAL: no genome column in --df-systems (tried "
                  f"'{args.systems_genome_col}', {SYS_FALLBACK_COLS})")
        sys.exit(1)
    if ggcol is None:
        utils.log(f"FATAL: no genome column in --df-genes (tried "
                  f"'{args.genes_genome_col}', {SYS_FALLBACK_COLS})")
        sys.exit(1)

    # ---- split systems per genome (reused layout) ------------------------ #
    split_systems_by_genome(sysdf, sgcol, args.out_defense_dir)

    # ---- select Prometheus systems --------------------------------------- #
    sl = args.system.lower()
    tcol = _pick(sysdf, "type"); stcol = _pick(sysdf, "subtype")
    sidcol = _pick(sysdf, "sys_id")
    if sidcol is None:
        utils.log("FATAL: --df-systems has no 'sys_id' column")
        sys.exit(1)
    prom_keys = set()       # (genome, sys_id)
    subtype_by_key = {}
    for _, r in sysdf.iterrows():
        t = str(r.get(tcol, "")).lower() if tcol else ""
        st = str(r.get(stcol, "")).lower() if stcol else ""
        if sl in t or sl in st:
            key = (str(r[sgcol]).strip(), str(r[sidcol]).strip())
            prom_keys.add(key)
            subtype_by_key[key] = str(r.get(stcol, "")) if stcol else ""

    if not prom_keys:
        open(args.out_candidates, "w").close()
        utils.write_empty(args.out_sources, SOURCES_HEADER)
        utils.write_empty(args.out_raw, RAW_HEADER)
        utils.log(f"No bacterial {args.system} systems found in {args.df_systems}.")
        return

    # ---- map Prometheus genes -> (genome, hit_id, score, subtype) -------- #
    g_sidcol = _pick(gdf, "sys_id")
    hitcol = _pick(gdf, "hit_id")
    scorecol = _pick(gdf, "hit_score", "score")
    if g_sidcol is None or hitcol is None:
        utils.log("FATAL: --df-genes needs 'sys_id' and 'hit_id' columns")
        sys.exit(1)

    per_genome_hits = {}    # genome -> list of (hit_id, score, subtype)
    raw_rows = []
    for _, r in gdf.iterrows():
        genome = str(r[ggcol]).strip()
        sid = str(r[g_sidcol]).strip()
        key = (genome, sid)
        if key not in prom_keys:
            continue
        hit_id = str(r[hitcol]).strip()
        if not hit_id or hit_id.lower() == "nan":
            continue
        try:
            score = float(r[scorecol]) if scorecol and str(r.get(scorecol, "")).strip() else None
        except (ValueError, TypeError):
            score = None
        subtype = subtype_by_key.get(key, "")
        per_genome_hits.setdefault(genome, []).append((hit_id, score, subtype))
        pid = utils.make_protein_id(genome, hit_id)
        raw_rows.append((pid, genome, "defensefinder", score, subtype))

    # ---- extract candidate proteins (and optional CDS) from per-genome files - #
    cand_records = []
    cds_records = []
    source_rows = []
    want_cds = bool(args.cds_dir and args.out_cds)
    n_missing_faa = n_missing_prot = 0
    for genome, hits in per_genome_hits.items():
        faa = os.path.join(args.faa_dir, f"{genome}.faa")
        if not os.path.exists(faa):
            # tolerate a couple of common alternative extensions
            alt = [os.path.join(args.faa_dir, f"{genome}.{e}") for e in ("fasta", "faa.gz", "fa")]
            faa = next((a for a in alt if os.path.exists(a)), faa)
        seqs = dict(utils.read_fasta(faa)) if os.path.exists(faa) else {}
        if not seqs:
            n_missing_faa += 1
            utils.log(f"WARNING: no proteome for bacterial genome '{genome}' "
                      f"({args.faa_dir}/{genome}.faa); its {args.system} proteins are skipped")
            continue
        best = {}   # hit_id -> (score, subtype)
        for hit_id, score, subtype in hits:
            if hit_id not in best:
                best[hit_id] = (score, subtype)
        cds_seqs = {}
        if want_cds:
            cfna = os.path.join(args.cds_dir, f"{genome}.fna")
            if not os.path.exists(cfna):
                cfna = next((os.path.join(args.cds_dir, f"{genome}.{e}")
                             for e in ("ffn", "fasta", "fa")
                             if os.path.exists(os.path.join(args.cds_dir, f"{genome}.{e}"))), cfna)
            cds_seqs = dict(utils.read_fasta(cfna)) if os.path.exists(cfna) else {}
        for hit_id, (score, subtype) in best.items():
            if hit_id not in seqs:
                n_missing_prot += 1
                continue
            pid = utils.make_protein_id(genome, hit_id)
            cand_records.append((pid, seqs[hit_id]))
            if want_cds and hit_id in cds_seqs:
                cds_records.append((pid, cds_seqs[hit_id]))
            source_rows.append({
                "protein_id": pid, "genome_id": genome, "detect_defensefinder": 1,
                "df_score": score if score is not None else "", "df_subtype": subtype})

    utils.write_fasta(sorted(set(cand_records)), args.out_candidates)
    if want_cds:
        utils.write_fasta(sorted(set(cds_records)), args.out_cds)
    if source_rows:
        pd.DataFrame(source_rows, columns=SOURCES_HEADER.split("\t")).to_csv(
            args.out_sources, sep="\t", index=False)
    else:
        utils.write_empty(args.out_sources, SOURCES_HEADER)
    if raw_rows:
        pd.DataFrame([{"protein_id": p, "genome_id": g, "method": m,
                       "score": s if s is not None else "", "subtype": st}
                      for p, g, m, s, st in raw_rows],
                     columns=RAW_HEADER.split("\t")).to_csv(args.out_raw, sep="\t", index=False)
    else:
        utils.write_empty(args.out_raw, RAW_HEADER)

    utils.log(f"Bacterial {args.system}: {len(cand_records)} candidate proteins across "
              f"{len(per_genome_hits)} genome(s) (missing faa={n_missing_faa}, "
              f"missing protein={n_missing_prot}).")


if __name__ == "__main__":
    main()
