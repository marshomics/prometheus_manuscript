"""Run BBMap statswrapper.sh (format=6) over one batch of genomes, in chunks.

Verified against BBMap 39.06 (jgi/AssemblyStats2.java, jgi/AssemblyStatsWrapper.java and
a real run of statswrapper.sh):
  * format=6 prints ONE header line per JVM, prefixed with '#':
      #n_scaffolds n_contigs scaf_bp contig_bp gap_pct scaf_N50 scaf_L50 ctg_N50 ctg_L50
       scaf_N90 scaf_L90 ctg_N90 ctg_L90 scaf_max ctg_max scaf_n_gt50K scaf_pct_gt50K
       gc_avg gc_std filename
    (`n_` prefix: N_UNDERSCORE defaults to true; we still pass n_=t explicitly; the
    wrapper adds addname=t, hence the trailing `filename` column).
  * BBMap's N/L naming is the reverse of the usual convention: ctg_L50 / scaf_L50 are the
    N50 LENGTHS and ctg_N50 / scaf_N50 the L50 COUNTS. formatKB() prints plain integer bp
    whenever FORMAT >= 2, so format=6 lengths are bp, never "KB"/"MB".
  * contigs = scaffolds split at runs of >= 10 N (n=10 default); contig_bp excludes gaps.
  * gc_avg is a FRACTION (0-1), not a percentage.
  * a missing/unreadable input aborts the JVM with a RuntimeException and exit code 1
    after printing the rows of the files before it.

The script therefore FAILS (exit 1) when statswrapper.sh is not on PATH, when a chunk
exits non-zero, when the output lacks the header/filename column, or when any genome of
the batch has no (or more than one) stats row. Output: one TSV, genome_id + every
statswrapper column, one row per genome, sorted by genome_id. A batch directory with no
genome at all (genuinely empty input) yields a header-only table.
"""
from __future__ import annotations
import argparse
import glob
import os
import shutil
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


class StatsError(Exception):
    pass


def _gid_from_filename(fname: str, ext) -> str:
    base = os.path.basename(fname.strip())
    if ext:
        suf = "." + ext.lstrip(".")
        if not base.endswith(suf):
            raise StatsError(f"statswrapper filename '{fname}' does not end with '{suf}'")
        return base[: -len(suf)]
    return base.rsplit(".", 1)[0]


def parse_statswrapper(text: str, ext=None):
    """Parse statswrapper format=6 (or 3) stdout.

    Returns (columns_excluding_filename, [(genome_id, {column: value})]). The header is
    the first non-blank line; it must contain a `filename` column. Any later line that
    starts with '#' (a second header) is skipped. A data line whose field count differs
    from the header is an error, not a silent skip.
    """
    header = None
    rows = []
    for ln_no, line in enumerate(text.splitlines(), 1):
        if not line.strip():
            continue
        if header is None:
            header = [h.strip() for h in line.lstrip().lstrip("#").split("\t")]
            if not any(h.lower() == "filename" for h in header):
                raise StatsError(f"statswrapper header has no 'filename' column: {header}")
            continue
        if line.lstrip().startswith("#"):
            continue
        vals = line.rstrip("\n").split("\t")
        if len(vals) != len(header):
            raise StatsError(f"statswrapper line {ln_no} has {len(vals)} fields, header has "
                             f"{len(header)}: {line[:200]!r}")
        rec = dict(zip(header, vals))
        fcol = next(c for c in header if c.lower() == "filename")
        gid = _gid_from_filename(rec.pop(fcol), ext)
        rows.append((gid, rec))
    cols = [c for c in (header or []) if c.lower() != "filename"]
    return cols, rows


def run_chunk(names, genomes_dir: str, ext: str, exe: str):
    cmd = [exe, "in=" + ",".join(names), "format=6", "n_=t"]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd=genomes_dir)
    if res.returncode != 0:
        raise StatsError(f"statswrapper.sh exited {res.returncode} on a chunk of {len(names)} "
                         f"genome(s) (first: {names[0]}); stderr tail:\n{res.stderr[-1500:]}")
    return parse_statswrapper(res.stdout, ext)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--ext", default="fna")
    ap.add_argument("--chunk", type=int, default=2000, help="genomes per JVM")
    ap.add_argument("--statswrapper", default="statswrapper.sh")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    try:
        _run(args)
    except StatsError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)


def _run(args) -> None:
    if args.chunk < 1:
        raise StatsError(f"--chunk must be >= 1 (got {args.chunk})")
    files = sorted(glob.glob(os.path.join(args.genomes_dir, f"*.{args.ext}")))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if not files:
        # Genuinely empty input (no genome in this batch): header-only table.
        utils.write_empty(args.out, "genome_id")
        utils.log(f"no *.{args.ext} genome in {args.genomes_dir}; wrote header-only BBMap stats")
        return
    exe = shutil.which(args.statswrapper)
    if exe is None:
        raise StatsError(f"'{args.statswrapper}' not found on PATH (is the bbmap env active?)")

    names = [os.path.basename(f) for f in files]
    expected = [_gid_from_filename(nm, args.ext) for nm in names]
    all_cols: list = []
    records: dict = {}
    for i in range(0, len(names), args.chunk):
        chunk = names[i:i + args.chunk]
        cols, rows = run_chunk(chunk, args.genomes_dir, args.ext, exe)
        for c in cols:
            if c not in all_cols:
                all_cols.append(c)
        want = {_gid_from_filename(nm, args.ext) for nm in chunk}
        got = [g for g, _ in rows]
        dups = sorted(g for g, k in Counter(got).items() if k > 1)
        if dups:
            raise StatsError(f"statswrapper returned duplicated rows for {dups[:5]}")
        missing = sorted(want - set(got))
        unexpected = sorted(set(got) - want)
        if missing or unexpected:
            raise StatsError(f"chunk {i // args.chunk}: {len(missing)} genome(s) without stats "
                             f"(e.g. {missing[:5]}), {len(unexpected)} unexpected (e.g. {unexpected[:5]})")
        records.update(rows)

    assert len(records) == len(expected)
    with open(args.out, "w") as out:
        out.write("genome_id\t" + "\t".join(all_cols) + "\n")
        for gid in sorted(records):
            rec = records[gid]
            out.write(gid + "\t" + "\t".join(rec.get(c, "") for c in all_cols) + "\n")
    utils.log(f"BBMap stats for {len(records)} genome(s) ({len(all_cols)} columns) -> {args.out}")


if __name__ == "__main__":
    main()
