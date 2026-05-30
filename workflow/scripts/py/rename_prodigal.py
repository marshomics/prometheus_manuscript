"""Rewrite Prodigal .faa and .fna headers to >{genome}_{N} (N = 1..), consistently.

Prodigal emits the protein (-a) and CDS nucleotide (-d) records for the same genes,
and the gene id (the first whitespace token of the header, e.g. `k141_3_12`) is
identical in both files. We number the proteins in file order as `{genome}_{i}`
(i from 1) and apply the SAME new id to the matching CDS record (matched by the
original gene id), so the i-th protein and its CDS share `{genome}_{i}`. The new
header consists of ONLY `{genome}_{i}` (no Prodigal coordinate description).
"""
from __future__ import annotations
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--fna", required=True)
    ap.add_argument("--genome", required=True)
    args = ap.parse_args()

    g = utils.sanitize(args.genome)
    faa = list(utils.read_fasta(args.faa))   # [(orig_gene_id, seq)] in gene order
    fna = list(utils.read_fasta(args.fna))

    # number proteins in order; remember orig_gene_id -> new id
    new_by_orig = {}
    new_faa = []
    for i, (orig, seq) in enumerate(faa, start=1):
        nid = f"{g}_{i}"
        new_by_orig[orig] = nid
        new_faa.append((nid, seq))
    utils.write_fasta(new_faa, args.faa)

    # apply the SAME mapping to the CDS file, matched by the original gene id, so
    # corresponding protein/CDS records carry identical headers.
    new_fna, missing = [], 0
    for orig, seq in fna:
        nid = new_by_orig.get(orig)
        if nid is None:
            missing += 1
            continue
        new_fna.append((nid, seq))
    utils.write_fasta(new_fna, args.fna)

    if missing:
        utils.log(f"WARNING: {missing} CDS record(s) in {args.fna} had no matching "
                  "protein id and were dropped (faa/fna kept consistent).")
    utils.log(f"renamed {len(new_faa)} proteins / {len(new_fna)} CDS to '{g}_N' "
              f"for genome {args.genome}")


if __name__ == "__main__":
    main()
