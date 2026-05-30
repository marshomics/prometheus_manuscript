"""Concatenate per-genome proteomes into a single FASTA with genome-prefixed IDs."""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

MAP_HEADER = "protein_id\tgenome_id\torig_id"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", required=True)
    ap.add_argument("--genes-dir", required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--suffix", default="faa")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out_faa)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_map)), exist_ok=True)

    try:
        samples = pd.read_csv(args.samples, sep="\t")
    except Exception:
        utils.write_empty(args.out_faa, "")
        utils.write_empty(args.out_map, MAP_HEADER)
        return

    if samples.empty or "genome_id" not in samples.columns:
        utils.write_empty(args.out_faa, "")
        utils.write_empty(args.out_map, MAP_HEADER)
        return

    all_records = []
    map_rows = []

    for _, row in samples.iterrows():
        gid = str(row["genome_id"]).strip()
        faa_path = os.path.join(args.genes_dir, gid, f"{gid}.{args.suffix}")
        for orig_id, seq in utils.read_fasta(faa_path):
            pid = utils.make_protein_id(gid, orig_id)
            all_records.append((pid, seq))
            map_rows.append({"protein_id": pid, "genome_id": gid, "orig_id": orig_id})

    # Write combined FASTA
    if all_records:
        utils.write_fasta(all_records, args.out_faa)
    else:
        # Write empty file (valid empty FASTA)
        os.makedirs(os.path.dirname(os.path.abspath(args.out_faa)), exist_ok=True)
        open(args.out_faa, "w").close()

    # Write map TSV
    if map_rows:
        df = pd.DataFrame(map_rows, columns=["protein_id", "genome_id", "orig_id"])
        df.to_csv(args.out_map, sep="\t", index=False)
    else:
        utils.write_empty(args.out_map, MAP_HEADER)

    utils.log(f"Wrote {len(all_records)} proteins from {len(samples)} genomes -> {args.out_faa}")


if __name__ == "__main__":
    main()
