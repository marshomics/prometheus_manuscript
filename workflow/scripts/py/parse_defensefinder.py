"""Extract Prometheus candidates from DefenseFinder's Prometheus system calls.

Detection is DefenseFinder-only (no PADLOC, no custom HMM sweep). DefenseFinder
applies its MacSyFinder genetic-architecture rule before emitting a system, so a
called Prometheus system already satisfies the model's mandatory-gene requirement;
this script simply collects the proteins of those systems for downstream
validation (architecture, phylogenetic placement, context, tiering).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

SOURCES_HEADER = "protein_id\tgenome_id\tdetect_defensefinder\tdf_score\tdf_subtype"
RAW_HEADER = "protein_id\tgenome_id\tmethod\tscore\tsubtype"


def load_pmap(pmap_path: str):
    """Load protein_map.tsv -> dict (genome_id, orig_id) -> protein_id."""
    mapping = {}
    try:
        df = pd.read_csv(pmap_path, sep="\t")
    except Exception as exc:
        utils.log(f"WARNING: could not read pmap ({exc}); falling back to make_protein_id")
        return mapping
    needed = {"genome_id", "orig_id", "protein_id"}
    if not needed.issubset(set(df.columns)):
        utils.log(f"WARNING: pmap missing columns {needed - set(df.columns)}; "
                  "falling back to make_protein_id")
        return mapping
    for _, row in df.iterrows():
        mapping[(str(row["genome_id"]), str(row["orig_id"]))] = str(row["protein_id"])
    return mapping


def lookup_pid(pmap, genome_id, orig_id):
    key = (genome_id, orig_id)
    return pmap.get(key, utils.make_protein_id(genome_id, orig_id))


def parse_defensefinder(df_dir: str, system: str, pmap: dict):
    """Return list of (protein_id, genome_id, score, subtype) for matching systems."""
    hits = []
    if not df_dir or not os.path.isdir(df_dir):
        return hits

    for genome_id in os.listdir(df_dir):
        gdir = os.path.join(df_dir, genome_id)
        if not os.path.isdir(gdir):
            continue
        sys_file = os.path.join(gdir, f"{genome_id}_defense_finder_systems.tsv")
        gene_file = os.path.join(gdir, f"{genome_id}_defense_finder_genes.tsv")

        matching_sys_ids = set()
        try:
            sdf = pd.read_csv(sys_file, sep="\t")
            sdf.columns = [c.strip() for c in sdf.columns]
        except Exception as exc:
            # A genome with a real Prometheus must never be silently dropped:
            # warn loudly so a malformed systems file is visible, not swallowed.
            utils.log(f"WARNING: could not read systems file for '{genome_id}' "
                      f"({exc}); this genome is skipped — inspect {sys_file}")
            continue
        # type OR subtype, case-insensitive substring (handles 'Prometheus',
        # 'PrometheusTypeI', etc. across DefenseFinder model-name variants).
        sys_lower = system.lower()
        for _, row in sdf.iterrows():
            t = str(row.get("type", "")).lower()
            st = str(row.get("subtype", "")).lower()
            if sys_lower in t or sys_lower in st:
                matching_sys_ids.add(str(row.get("sys_id", "")).strip())

        if not matching_sys_ids:
            continue

        try:
            gdf = pd.read_csv(gene_file, sep="\t")
            gdf.columns = [c.strip() for c in gdf.columns]
        except Exception as exc:
            # The systems file DID call a Prometheus system here, so losing the
            # genes file silently would drop a real locus. Warn explicitly.
            utils.log(f"WARNING: '{genome_id}' has a matching {system} system but its "
                      f"genes file could not be read ({exc}); inspect {gene_file}")
            continue

        if "hit_id" not in gdf.columns:
            utils.log(f"WARNING: '{genome_id}' genes file lacks a 'hit_id' column; "
                      f"cannot map proteins for its {system} system(s) — inspect {gene_file}")
            continue

        for _, row in gdf.iterrows():
            sid = str(row.get("sys_id", "")).strip()
            if sid not in matching_sys_ids:
                continue
            orig_id = str(row.get("hit_id", "")).strip()
            if not orig_id or orig_id.lower() == "nan":
                continue  # no protein id -> nothing to extract; do not emit genome__
            score = row.get("hit_score", "")
            try:
                score = float(score)
            except (ValueError, TypeError):
                score = None
            subtype = ""
            try:
                smatch = sdf[sdf["sys_id"].astype(str).str.strip() == sid]
                if not smatch.empty:
                    subtype = str(smatch.iloc[0].get("subtype", ""))
            except Exception:
                pass
            pid = lookup_pid(pmap, genome_id, orig_id)
            hits.append((pid, genome_id, score, subtype))

    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", required=True)
    ap.add_argument("--system", required=True)
    ap.add_argument("--out-candidates", required=True)
    ap.add_argument("--out-sources", required=True)
    ap.add_argument("--out-raw", required=True)
    args = ap.parse_args()

    for path in (args.out_candidates, args.out_sources, args.out_raw):
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    pmap = load_pmap(args.pmap)
    df_hits = parse_defensefinder(args.df_dir, args.system, pmap)

    # best (score, subtype) per protein id
    df_by_pid = {}
    for pid, gid, score, subtype in df_hits:
        if pid not in df_by_pid:
            df_by_pid[pid] = (score, subtype)

    if not df_by_pid:
        open(args.out_candidates, "w").close()
        utils.write_empty(args.out_sources, SOURCES_HEADER)
        utils.write_empty(args.out_raw, RAW_HEADER)
        utils.log("No DefenseFinder Prometheus candidates found.")
        return

    # extract candidate sequences from the combined proteome
    proteome = dict(utils.read_fasta(args.proteome))
    candidate_records = [(pid, proteome[pid]) for pid in sorted(df_by_pid) if pid in proteome]
    utils.write_fasta(candidate_records, args.out_candidates)

    # candidate_sources.tsv (DefenseFinder is the sole detector)
    source_rows = []
    for pid in sorted(df_by_pid):
        score, subtype = df_by_pid[pid]
        source_rows.append({
            "protein_id": pid,
            "genome_id": utils.genome_of(pid),
            "detect_defensefinder": 1,
            "df_score": score if score is not None else "",
            "df_subtype": subtype,
        })
    pd.DataFrame(source_rows, columns=SOURCES_HEADER.split("\t")).to_csv(
        args.out_sources, sep="\t", index=False)

    # raw_hits.tsv (one row per DefenseFinder gene hit)
    raw_rows = [{
        "protein_id": pid, "genome_id": gid, "method": "defensefinder",
        "score": score if score is not None else "", "subtype": subtype,
    } for pid, gid, score, subtype in df_hits]
    if raw_rows:
        pd.DataFrame(raw_rows, columns=RAW_HEADER.split("\t")).to_csv(
            args.out_raw, sep="\t", index=False)
    else:
        utils.write_empty(args.out_raw, RAW_HEADER)

    utils.log(f"DefenseFinder Prometheus candidates: {len(df_by_pid)} proteins "
              f"across {len(set(utils.genome_of(p) for p in df_by_pid))} genomes.")


if __name__ == "__main__":
    main()
