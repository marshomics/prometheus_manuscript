"""Extract Prometheus candidates from DefenseFinder's Prometheus system calls.

Detection is DefenseFinder-only (no PADLOC, no custom HMM sweep). DefenseFinder applies
its MacSyFinder genetic-architecture rule before emitting a system, so a called
Prometheus system already satisfies the model's mandatory-gene requirement
(defense-finder-models Prometheus.xml: one mandatory gene `Prometheus__ProA`,
min_mandatory_genes_required=1).

Output format (defense-finder v2.0.x/v3.0.0, defense_finder_posttreat):
  *_defense_finder_systems.tsv : sys_id, type, subtype, activity, sys_beg, sys_end,
                                 protein_in_syst, genes_count, name_of_profiles_in_sys
  *_defense_finder_genes.tsv   : MacSyFinder 2.1.4 TsvSystemSerializer columns
                                 (replicon, hit_id, gene_name, hit_pos, model_fqn, sys_id,
                                 ..., hit_i_eval, hit_score, hit_profile_cov, hit_seq_cov,
                                 ...) + type, subtype, activity
  where type/subtype are the last two components of model_fqn
  (defense-finder-models/DefenseFinder/Prometheus/Prometheus -> Prometheus/Prometheus).

Selection: a system matches when its type OR subtype EQUALS --system (case-insensitive;
no substring matching). Candidates are the proteins of its genes whose gene_name equals
--gene (default "<system>__ProA"; pass --gene "" to take every gene of the system).
Genomes are taken from --samples (the analysis genome set), never from whatever
directories happen to exist under --df-dir. Every candidate protein must exist in the
combined proteome; a missing one is an id mismatch (e.g. ingested DefenseFinder output
run on a different proteome) and fails the step.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

SOURCES_COLUMNS = ["protein_id", "genome_id", "detect_defensefinder", "df_score",
                   "df_subtype",
                   # additive (design 1.1: system instance, profile, i-evalue, coverage)
                   "df_sys_id", "df_gene_name", "df_i_eval", "df_profile_cov", "df_seq_cov"]
SOURCES_HEADER = "\t".join(SOURCES_COLUMNS[:5])     # SCHEMA base columns (compatibility)
RAW_HEADER = "protein_id\tgenome_id\tmethod\tscore\tsubtype"
SYSTEMS_REQUIRED = {"sys_id", "type", "subtype"}
GENES_REQUIRED = {"hit_id", "sys_id"}


class InputError(RuntimeError):
    pass


def load_pmap(pmap_path: str) -> dict:
    """protein_map.tsv -> {(genome_id, orig_id): protein_id} (required, never optional)."""
    df = pd.read_csv(pmap_path, sep="\t", dtype=str, keep_default_na=False)
    needed = {"genome_id", "orig_id", "protein_id"}
    if not needed.issubset(df.columns):
        raise InputError(f"{pmap_path} lacks columns {sorted(needed - set(df.columns))}")
    return {(g, o): p for g, o, p in zip(df["genome_id"], df["orig_id"], df["protein_id"])}


def genome_list(samples: str, df_dir: str) -> list:
    if samples:
        df = pd.read_csv(samples, sep="\t", dtype=str, keep_default_na=False, comment="#")
        if "genome_id" not in df.columns:
            raise InputError(f"{samples} has no genome_id column")
        return [g.strip() for g in df["genome_id"] if g.strip()]
    # legacy behaviour (no --samples): every genome directory under df_dir
    if not df_dir or not os.path.isdir(df_dir):
        return []
    return sorted(d for d in os.listdir(df_dir) if os.path.isdir(os.path.join(df_dir, d)))


def read_table(path: str, required: set, kind: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise InputError(f"missing DefenseFinder {kind} table {path}")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    if not required.issubset(df.columns):
        raise InputError(f"{path} lacks required {kind} columns "
                         f"{sorted(required - set(df.columns))} (header: {list(df.columns)})")
    return df


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def parse_defensefinder(df_dir: str, system: str, gene: str, pmap: dict, genomes: list):
    """Return (hits, total_systems_all_types, n_genomes).

    hits: dicts with protein_id, genome_id, score, subtype, sys_id, gene_name, i_eval,
    profile_cov, seq_cov -- one per matching DefenseFinder gene row.
    """
    sys_l, gene_l = system.strip().lower(), gene.strip().lower()
    hits, total, missing = [], 0, []
    for g in genomes:
        gdir = os.path.join(df_dir, g)
        sdf = read_table(os.path.join(gdir, f"{g}_defense_finder_systems.tsv"),
                         SYSTEMS_REQUIRED, "systems")
        total += len(sdf)
        match = sdf[(sdf["type"].str.strip().str.lower() == sys_l)
                    | (sdf["subtype"].str.strip().str.lower() == sys_l)]
        if match.empty:
            continue
        sub_of = dict(zip(match["sys_id"].str.strip(), match["subtype"].str.strip()))
        gdf = read_table(os.path.join(gdir, f"{g}_defense_finder_genes.tsv"),
                         GENES_REQUIRED, "genes")
        rows = gdf[gdf["sys_id"].str.strip().isin(sub_of)]
        if gene_l and "gene_name" in rows.columns:
            rows = rows[rows["gene_name"].str.strip().str.lower() == gene_l]
        if rows.empty:
            raise InputError(f"genome {g}: {len(sub_of)} {system} system(s) in the systems "
                             f"table but no matching gene rows (gene_name == '{gene}') in "
                             "the genes table")
        for r in rows.to_dict("records"):
            hid = str(r.get("hit_id", "")).strip()
            if not hid:
                raise InputError(f"genome {g}: {system} gene row without hit_id")
            pid = pmap.get((g, hid))
            if pid is None:
                missing.append(f"{g}:{hid}")
                continue
            hits.append({
                "protein_id": pid, "genome_id": g,
                "score": as_float(r.get("hit_score")),
                "subtype": sub_of.get(str(r["sys_id"]).strip(), ""),
                "sys_id": str(r["sys_id"]).strip(),
                "gene_name": str(r.get("gene_name", "")).strip(),
                "i_eval": as_float(r.get("hit_i_eval")),
                "profile_cov": as_float(r.get("hit_profile_cov")),
                "seq_cov": as_float(r.get("hit_seq_cov")),
            })
    if missing:
        raise InputError(
            f"{len(missing)} DefenseFinder {system} hit_id(s) are absent from the pipeline "
            f"proteome (protein_map.tsv), e.g. {missing[:5]}. The DefenseFinder output was "
            "produced on different protein ids (ingested output from another gene call, or "
            "a renamed proteome). Re-run DefenseFinder on genes/{g}/{g}.faa.")
    return hits, total, len(genomes)


def fmt(v):
    return "" if v is None else v


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", required=True)
    ap.add_argument("--system", required=True)
    ap.add_argument("--strict", default="True",
                    help="FAIL (exit 1) if the proteome is non-empty but DefenseFinder "
                         "produced zero systems of ANY type (broken detector / model set)")
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-sources", required=True)
    ap.add_argument("--out-raw", required=True)
    ap.add_argument("--samples", default="",
                    help="analysis genome table (genome_id); default: every dir in --df-dir")
    ap.add_argument("--gene", default=None,
                    help="gene_name of the candidate protein within the system "
                         "(default '<system>__ProA'; '' = every gene of the system)")
    args = ap.parse_args()
    gene = f"{args.system}__ProA" if args.gene is None else args.gene

    for path in (args.out_candidates, args.out_sources, args.out_raw):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    try:
        pmap = load_pmap(args.pmap)
        genomes = genome_list(args.samples, args.df_dir)
        df_hits, total_systems, n_genomes = parse_defensefinder(
            args.df_dir, args.system, gene, pmap, genomes)
    except InputError as exc:
        raise SystemExit(f"[parse_defensefinder] FATAL: {exc}")

    # --- strict detection sanity (complements the defensefinder_check preflight) ---
    strict = str(args.strict).strip().lower() in ("true", "1", "yes")
    n_prot = len(pmap)
    if strict and n_prot > 0 and n_genomes > 0 and total_systems == 0:
        raise SystemExit(
            "[parse_defensefinder] FATAL: DefenseFinder produced ZERO systems of any type "
            f"across {n_genomes} genome(s) / {n_prot} proteins. This indicates a broken "
            "detector or wrong model set, not a Prometheus-free dataset. Check the "
            "defensefinder_check sentinel and the per-genome logs. Set "
            "prometheus.strict_detection=false to override (not recommended).")

    # one row per protein: best (highest-scoring) DefenseFinder hit
    best = {}
    for h in df_hits:
        cur = best.get(h["protein_id"])
        if cur is None or (h["score"] or float("-inf")) > (cur["score"] or float("-inf")):
            best[h["protein_id"]] = h

    if not best:
        open(args.out_candidates, "w").close()
        utils.write_empty(args.out_sources, "\t".join(SOURCES_COLUMNS))
        utils.write_empty(args.out_raw, RAW_HEADER)
        utils.log(f"No DefenseFinder {args.system} candidates in {n_genomes} genome(s) "
                  f"({total_systems} systems of any type).")
        return

    proteome = dict(utils.read_fasta(args.proteome))
    absent = [p for p in best if p not in proteome]
    if absent:
        raise SystemExit(f"[parse_defensefinder] FATAL: {len(absent)} candidate(s) are in "
                         f"protein_map.tsv but not in {args.proteome} (e.g. {absent[:3]}); "
                         "the combined proteome is stale or truncated.")
    utils.write_fasta([(p, proteome[p]) for p in sorted(best)], args.out_candidates)

    pd.DataFrame([{
        "protein_id": p, "genome_id": best[p]["genome_id"], "detect_defensefinder": 1,
        "df_score": fmt(best[p]["score"]), "df_subtype": best[p]["subtype"],
        "df_sys_id": best[p]["sys_id"], "df_gene_name": best[p]["gene_name"],
        "df_i_eval": fmt(best[p]["i_eval"]), "df_profile_cov": fmt(best[p]["profile_cov"]),
        "df_seq_cov": fmt(best[p]["seq_cov"]),
    } for p in sorted(best)], columns=SOURCES_COLUMNS).to_csv(
        args.out_sources, sep="\t", index=False)

    pd.DataFrame([{
        "protein_id": h["protein_id"], "genome_id": h["genome_id"], "method": "defensefinder",
        "score": fmt(h["score"]), "subtype": h["subtype"],
    } for h in df_hits], columns=RAW_HEADER.split("\t")).to_csv(
        args.out_raw, sep="\t", index=False)

    utils.log(f"DefenseFinder {args.system} candidates: {len(best)} proteins across "
              f"{len(set(h['genome_id'] for h in best.values()))} of {n_genomes} genomes "
              f"({total_systems} systems of any type).")


if __name__ == "__main__":
    main()
