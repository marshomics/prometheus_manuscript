"""Split the raw input genome set into deterministic batches of symlinks.

Each batch directory gets `genomes/{genome_id}.fna` symlinks pointing at the real
FASTAs ({genomes_dir}/{filename}). Normalising names to {genome_id}.fna means
GTDB-Tk / CheckM2 / BBMap all report stable genome_ids, with no filename->id
remapping downstream. Batching keeps GTDB-Tk runs to ~20k genomes each.

The batch assignment is deterministic (genomes sorted by id, chunked by size), so
the number of batches matches common.smk's N_BATCHES.
"""
from __future__ import annotations
import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--genomes-dir", default="")
    ap.add_argument("--id-col", default="genome_id")
    ap.add_argument("--file-col", default="filename")
    ap.add_argument("--input-ext", default="fasta", help="extension of raw input FASTAs")
    ap.add_argument("--batch-size", type=int, default=20000)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    import pandas as pd
    df = pd.read_csv(args.metadata, sep="\t", dtype=str)
    if args.id_col not in df.columns:
        utils.log(f"FATAL: id column '{args.id_col}' not in metadata"); sys.exit(1)
    fcol = args.file_col if args.file_col in df.columns else None
    df = df[[c for c in (args.id_col, fcol) if c]].dropna(subset=[args.id_col])
    df = df.drop_duplicates(subset=[args.id_col]).sort_values(args.id_col).reset_index(drop=True)

    n = len(df)
    nb = max(1, math.ceil(n / args.batch_size)) if n else 0
    utils.log(f"{n} genomes -> {nb} batch(es) of <= {args.batch_size}")
    if args.genomes_dir == "":
        utils.log("WARNING: input.genomes_dir is empty; symlinks will be created but "
                  "will dangle until you set it. Real runs need it pointed at the FASTAs.")

    made = dangling = 0
    for i in range(nb):
        b = f"{i:04d}"
        gdir = os.path.join(args.outdir, f"batch_{b}", "genomes")
        os.makedirs(gdir, exist_ok=True)
        sub = df.iloc[i * args.batch_size:(i + 1) * args.batch_size]
        for _, row in sub.iterrows():
            gid = utils.sanitize(str(row[args.id_col]))
            fname = str(row[fcol]) if fcol and not pd.isna(row[fcol]) else f"{gid}.{args.input_ext}"
            src = os.path.abspath(os.path.join(args.genomes_dir, fname)) if args.genomes_dir \
                else os.path.abspath(fname)
            dst = os.path.join(gdir, f"{gid}.fna")
            try:
                if os.path.islink(dst) or os.path.exists(dst):
                    os.remove(dst)
                os.symlink(src, dst)
                made += 1
                if not os.path.exists(src):
                    dangling += 1
            except OSError as exc:
                utils.log(f"WARNING: could not link {dst} -> {src} ({exc})")
        open(os.path.join(args.outdir, f"batch_{b}", ".ready"), "w").close()

    utils.log(f"staged {made} symlinks ({dangling} currently dangling); wrote {nb} batch .ready files")


if __name__ == "__main__":
    main()
