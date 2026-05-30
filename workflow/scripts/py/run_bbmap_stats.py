"""Run BBMap statswrapper.sh (format=6) over a batch of genomes, in chunks.

statswrapper.sh processes many FASTAs in one JVM and prints a tab-delimited table
with a '#'-prefixed header (format=6) plus a `filename` column. We chunk the file
list (default 2000/JVM) to stay under the shell argument limit while keeping JVM
startups few, parse the output BY HEADER NAME (robust to column-order/version
differences), recover genome_id from the filename basename, and write one combined
TSV (genome_id + every stat column). `n_=t` is passed so contig/scaffold counts are
named n_contigs / n_scaffolds.
"""
from __future__ import annotations
import argparse
import glob
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


def parse_statswrapper(text: str):
    """Return (colnames_excluding_filename, list_of_(genome_id, {col:val}))."""
    header = None
    rows = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if header is None and ("n_contigs" in line or "n_scaffolds" in line
                               or line.lstrip().startswith("#")):
            header = line.lstrip("#").strip().split("\t")
            continue
        if header is None or line.lstrip().startswith("#"):
            continue
        vals = line.split("\t")
        if len(vals) < len(header):
            continue
        rec = dict(zip(header, vals))
        # filename column -> genome_id
        fcol = next((c for c in header if c.lower() == "filename"), None)
        if fcol is None:
            continue
        gid = os.path.basename(rec.pop(fcol)).rsplit(".", 1)[0]
        rows.append((gid, rec))
    cols = [c for c in (header or []) if c.lower() != "filename"]
    return cols, rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--ext", default="fna")
    ap.add_argument("--chunk", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.genomes_dir, f"*.{args.ext}")))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if not files:
        utils.write_empty(args.out, "genome_id")
        utils.log("no genomes in batch; wrote header-only bbmap stats")
        return

    all_cols: list[str] = []
    records: dict[str, dict] = {}
    for i in range(0, len(files), args.chunk):
        chunk = files[i:i + args.chunk]
        cmd = ["statswrapper.sh", "in=" + ",".join(chunk), "format=6", "n_=t"]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            utils.log("statswrapper.sh not on PATH; writing header-only output")
            utils.write_empty(args.out, "genome_id")
            return
        if res.returncode != 0:
            utils.log(f"WARNING: statswrapper chunk {i//args.chunk} rc={res.returncode}; "
                      f"stderr tail: {res.stderr[-300:]}")
        cols, rows = parse_statswrapper(res.stdout)
        for c in cols:
            if c not in all_cols:
                all_cols.append(c)
        for gid, rec in rows:
            records[gid] = rec

    with open(args.out, "w") as out:
        out.write("genome_id\t" + "\t".join(all_cols) + "\n")
        for gid in sorted(records):
            rec = records[gid]
            out.write(gid + "\t" + "\t".join(rec.get(c, "") for c in all_cols) + "\n")
    utils.log(f"BBMap stats for {len(records)} genomes ({len(all_cols)} columns)")


if __name__ == "__main__":
    main()
