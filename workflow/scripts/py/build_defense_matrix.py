"""Build genome x defense-system presence/absence matrix from DefenseFinder outputs."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", required=True)
    ap.add_argument("--df-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # Load master table to get genome list
    try:
        master = pd.read_csv(args.master, sep="\t")
    except Exception:
        utils.write_empty(args.out, "genome_id\tPrometheus")
        return

    if master.empty or "genome_id" not in master.columns:
        utils.write_empty(args.out, "genome_id\tPrometheus")
        return

    genomes = master["genome_id"].astype(str).tolist()

    # Collect all system types across all genomes
    genome_systems = {}  # genome_id -> set of system type strings
    all_system_types = set()

    for gid in genomes:
        gdir = os.path.join(args.df_dir, gid)
        sys_file = os.path.join(gdir, f"{gid}_defense_finder_systems.tsv")
        types = set()
        if os.path.exists(sys_file):
            try:
                df = pd.read_csv(sys_file, sep="\t")
                if not df.empty and "type" in df.columns:
                    for t in df["type"].dropna().astype(str):
                        t = t.strip()
                        if t:
                            types.add(t)
            except Exception:
                pass
        genome_systems[gid] = types
        all_system_types.update(types)

    # Ensure Prometheus column always present
    all_system_types.add("Prometheus")
    sorted_types = sorted(all_system_types)

    # Build matrix
    rows = []
    for gid in genomes:
        row = {"genome_id": gid}
        present = genome_systems.get(gid, set())
        for stype in sorted_types:
            row[stype] = 1 if stype in present else 0
        rows.append(row)

    if not rows:
        utils.write_empty(args.out, "genome_id\tPrometheus")
        return

    cols = ["genome_id"] + sorted_types
    out_df = pd.DataFrame(rows, columns=cols)
    out_df.to_csv(args.out, sep="\t", index=False)

    utils.log(
        f"Defense matrix: {len(rows)} genomes x {len(sorted_types)} systems -> {args.out}"
    )


if __name__ == "__main__":
    main()
