"""Split the Phase-1 input genome set into deterministic batches of symlinks.

Genome set = the data rows of the input metadata with a non-empty id, read exactly as
common.smk counts them (`#` comment lines and blank lines ignored, ids stripped,
duplicated ids are an error). For every genome the script

  * validates the id: `[A-Za-z0-9][A-Za-z0-9_.-]*` and no `__`. The id becomes a file
    stem, a Snakemake wildcard, a newick tip and the prefix of every gene id
    (`{genome}_{contig#}_{gene#}`); `__` is the combined-id separator (utils.ID_SEP);
  * resolves its FASTA (`{genomes_dir}/{filename}`, or `{genomes_dir}/{id}.{input_ext}`
    when the filename column is absent or blank) and FAILS if the file is missing,
    empty, compressed or not FASTA. A genome is never silently dropped;
  * assigns the genomes, sorted by id, to batches of --batch-size and stages
    `batch_XXXX/genomes/{genome_id}.fna` symlinks, so CheckM2 (`Name` = file stem) and
    BBMap (`filename` column) report the genome_id;
  * writes `batch_XXXX/gtdbtk_batchfile.tsv` (`<staged fasta>\\tuser_{genome_id}`).
    GTDB-Tk 2.4.0 aborts the whole run when a user genome id equals a GTDB reference
    accession (e.g. `GCA_000008085.1`, see gtdbtk/main.py `_genomes_to_process`); the
    `user_` prefix makes every id legal and the gtdbtk_batch rule strips it again;
  * writes `manifest.tsv` (genome_id, batch, source_path, staged_path, gtdbtk_id) and one
    `.ready` marker per batch.

With --expected-n (the count common.smk computed) the script fails when its own count
differs, so the `.ready` files it writes are exactly the DAG's BATCHES.
"""
from __future__ import annotations
import argparse
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

GTDBTK_ID_PREFIX = "user_"
STAGED_EXT = "fna"
ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-]*$")
MANIFEST_COLS = ["genome_id", "batch", "source_path", "staged_path", "gtdbtk_id"]


class InputError(Exception):
    pass


def check_no_inline_hash(path: str, err=None):
    """pandas comment='#' (as common.smk reads this table) drops everything after a '#'
    ANYWHERE in a line, silently truncating ids, filenames or the niche columns of that
    row. Only whole-line comments are allowed: a '#' inside a data row is an error."""
    err = err or InputError
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            s = line.strip()
            if s and not s.startswith("#") and "#" in line:
                raise err(f"{path}:{n}: '#' inside a data row would silently truncate it "
                          f"(the table is read with comment='#'): {line.rstrip()[:120]!r}")


def read_genome_ids(path: str, id_col: str, file_col: str):
    """Return [(genome_id, filename_or_'')] in file order, parsed exactly like
    common.smk._count_input_genomes (comment='#', blank lines skipped, ids stripped,
    empty ids excluded, duplicates an error)."""
    import pandas as pd
    check_no_inline_hash(path)
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#", skip_blank_lines=True)
    if id_col not in df.columns:
        raise InputError(f"input metadata {path} has no '{id_col}' column (config input.id_column)")
    ids = df[id_col].fillna("").astype(str).str.strip()
    n_blank = int((ids == "").sum())
    if n_blank:
        utils.log(f"WARNING: {n_blank} metadata row(s) have an empty '{id_col}' and are not genomes "
                  "(excluded, as in common.smk).")
    fnames = (df[file_col].fillna("").astype(str).str.strip() if file_col in df.columns
              else pd.Series([""] * len(df), index=df.index))
    keep = ids != ""
    ids, fnames = ids[keep], fnames[keep]
    dup = sorted(ids[ids.duplicated()].unique().tolist())
    if dup:
        raise InputError(f"{len(dup)} duplicated genome id(s) in {path}: {dup[:5]}")
    return list(zip(ids.tolist(), fnames.tolist()))


def bad_id_reason(gid: str):
    if not ID_RE.match(gid):
        return "only [A-Za-z0-9_.-] allowed and it must start with a letter or digit"
    if utils.ID_SEP in gid:
        return f"contains '{utils.ID_SEP}' (the combined protein-id separator)"
    return None


def fasta_problem(path: str):
    """None if `path` looks like an uncompressed, non-empty FASTA file, else a reason."""
    if not os.path.exists(path):
        return "missing" + (" (dangling symlink)" if os.path.islink(path) else "")
    if not os.path.isfile(path):
        return "not a regular file"
    if os.path.getsize(path) == 0:
        return "empty file"
    with open(path, "rb") as fh:
        head = fh.read(4096)
    if head[:2] == b"\x1f\x8b":
        return "gzip-compressed (decompress first: Prodigal cannot read .gz)"
    if head[:3] == b"BZh" or head[:6] == b"\xfd7zXZ\x00" or head[:4] == b"\x28\xb5\x2f\xfd":
        return "compressed (bzip2/xz/zstd); decompress first"
    first = head.lstrip()[:1]
    if first != b">":
        return "not FASTA (first non-blank character is not '>')"
    return None


def stage_link(dst: str, src: str) -> None:
    """Create/refresh the symlink dst -> src. Never overwrites a non-symlink file."""
    if os.path.islink(dst):
        if os.readlink(dst) == src:
            return
        os.remove(dst)
    elif os.path.lexists(dst):
        raise InputError(f"refusing to overwrite non-symlink file {dst}")
    os.symlink(src, dst)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--genomes-dir", default="")
    ap.add_argument("--id-col", default="genome_id")
    ap.add_argument("--file-col", default="filename")
    ap.add_argument("--input-ext", default="fasta", help="extension of raw input FASTAs")
    ap.add_argument("--batch-size", type=int, default=20000)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--expected-n", type=int, default=None,
                    help="genome count computed by common.smk; a mismatch is an error")
    ap.add_argument("--manifest", default=None, help="default: {outdir}/manifest.tsv")
    ap.add_argument("--gtdbtk-id-prefix", default=GTDBTK_ID_PREFIX,
                    help="prefix of the GTDB-Tk genome ids in the batchfile (stripped by gtdbtk_batch)")
    args = ap.parse_args()

    try:
        _run(args)
    except InputError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)


def _run(args) -> None:
    if args.batch_size < 1:
        raise InputError(f"--batch-size must be >= 1 (got {args.batch_size})")
    pre = args.gtdbtk_id_prefix
    if not ID_RE.match(pre) or pre.upper().startswith(("RS_", "GB_", "GCA", "GCF")):
        raise InputError(f"--gtdbtk-id-prefix {pre!r} must be a non-empty [A-Za-z0-9_.-] string "
                         "that cannot collide with GTDB accessions (RS_/GB_/GCA_/GCF_)")
    outdir = os.path.abspath(args.outdir)
    if " " in outdir:
        raise InputError(f"batch directory '{outdir}' contains a space; GTDB-Tk 2.4.0 rejects "
                         "genome paths with spaces (main.py _verify_file_path)")
    manifest = args.manifest or os.path.join(args.outdir, "manifest.tsv")

    rows = read_genome_ids(args.metadata, args.id_col, args.file_col)
    n = len(rows)
    if args.expected_n is not None and n != args.expected_n:
        raise InputError(f"{n} genomes in {args.metadata} but common.smk expected {args.expected_n}; "
                         "the metadata changed after the DAG was built -- re-run snakemake")

    bad_ids = [(g, bad_id_reason(g)) for g, _ in rows if bad_id_reason(g)]
    if bad_ids:
        ex = "; ".join(f"'{g}': {r}" for g, r in bad_ids[:10])
        raise InputError(f"{len(bad_ids)} invalid genome id(s) (rename them in the input metadata): {ex}")

    if not args.genomes_dir:
        utils.log("NOTE: input.genomes_dir is empty; FASTA paths are resolved relative to the "
                  "working directory.")
    resolved, problems = [], []
    for gid, fname in rows:
        fname = fname or f"{gid}.{args.input_ext}"
        src = os.path.abspath(os.path.join(args.genomes_dir, fname) if args.genomes_dir else fname)
        why = fasta_problem(src)
        if why:
            problems.append(f"{gid}: {src} ({why})")
        resolved.append((gid, src))
    if problems:
        raise InputError(f"{len(problems)} input genome FASTA(s) unusable:\n  " + "\n  ".join(problems[:25])
                         + ("\n  ..." if len(problems) > 25 else ""))

    resolved.sort(key=lambda t: t[0])  # deterministic: genomes sorted by id
    nb = math.ceil(n / args.batch_size) if n else 0
    utils.log(f"{n} genomes -> {nb} batch(es) of <= {args.batch_size}")
    os.makedirs(args.outdir, exist_ok=True)

    man_rows = []
    for i in range(nb):
        b = f"{i:04d}"
        bdir = os.path.join(args.outdir, f"batch_{b}")
        gdir = os.path.join(bdir, "genomes")
        os.makedirs(gdir, exist_ok=True)
        sub = resolved[i * args.batch_size:(i + 1) * args.batch_size]
        want = {f"{gid}.{STAGED_EXT}" for gid, _ in sub}
        stale = 0
        for name in sorted(os.listdir(gdir)):
            p = os.path.join(gdir, name)
            if name in want:
                continue
            if os.path.islink(p):
                os.remove(p)        # a genome no longer in this batch (pipeline-owned symlink)
                stale += 1
            else:
                raise InputError(f"unexpected non-symlink entry {p} in a pipeline staging directory")
        if stale:
            utils.log(f"batch_{b}: removed {stale} stale symlink(s) from a previous batch layout")
        with open(os.path.join(bdir, "gtdbtk_batchfile.tsv"), "w") as bf:
            for gid, src in sub:
                dst = os.path.join(gdir, f"{gid}.{STAGED_EXT}")
                stage_link(dst, src)
                gtdb_id = args.gtdbtk_id_prefix + gid
                bf.write(f"{os.path.abspath(dst)}\t{gtdb_id}\n")
                man_rows.append((gid, b, src, os.path.abspath(dst), gtdb_id))
        with open(os.path.join(bdir, ".ready"), "w") as fh:
            fh.write(f"{len(sub)}\t{sub[0][0]}\t{sub[-1][0]}\n")

    # Report (never delete) batch directories left over from a larger earlier run; the
    # merge step reads only the batches of the current DAG, so they cannot leak in.
    extra = sorted(d for d in os.listdir(args.outdir)
                   if re.fullmatch(r"batch_\d{4}", d) and int(d[6:]) >= nb)
    if extra:
        utils.log(f"WARNING: {len(extra)} stale batch dir(s) from an earlier run are ignored: {extra[:5]}")

    os.makedirs(os.path.dirname(os.path.abspath(manifest)), exist_ok=True)
    with open(manifest, "w") as fh:
        fh.write("\t".join(MANIFEST_COLS) + "\n")
        for r in man_rows:
            fh.write("\t".join(r) + "\n")
    utils.log(f"staged {len(man_rows)} symlinks in {nb} batch(es); manifest -> {manifest}")


if __name__ == "__main__":
    main()
