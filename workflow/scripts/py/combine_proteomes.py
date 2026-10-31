"""Concatenate FASTA files into one FASTA with newick-safe, collision-checked ids.

Proteome mode (rules combined_proteome; combined_cds with --suffix fna):
    --samples TSV(genome_id) --genes-dir DIR --out-faa F --out-map M [--suffix faa]
  reads {genes-dir}/{g}/{g}.{suffix} for every genome_id and writes ids
  `genome__orig` (utils.make_protein_id) plus protein_map.tsv (protein_id, genome_id,
  orig_id). A missing per-genome file, an unreadable samples table or a duplicated id
  is an ERROR (Contract v2 item 4); an empty per-genome file (no genes) contributes 0.

Reference mode (rule prometheus_validation_tree):
    --candidates FAA --ref Prometheus_ref=FAA --ref Paralog_ref=FAA
    --out-faa ALL --out-labels LABELS [--out-map MAP]
  writes candidates (ids kept; they must already be newick-safe) + references with ids
  passed through utils.sanitize ([^A-Za-z0-9_.-] -> '_'), and LABELS as the headerless
  `<id>\\t<group>` table for the references. IQ-TREE 2.3.6 itself rewrites ':' '=' ','
  ';' '(' ')' and quotes in names (it keeps '|' and '/'), so sanitising up front keeps
  the alignment, the tree tips and the label table identical. Collisions after
  sanitising and an empty reference set (with candidates present) are errors.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import pandas as pd

MAP_HEADER = "protein_id\tgenome_id\torig_id"
REF_MAP_HEADER = "tree_id\tgroup\torig_id"


def proteome_mode(args) -> None:
    if not args.samples or not args.genes_dir:
        raise SystemExit("[combine_proteomes] FATAL: proteome mode needs --samples and --genes-dir")
    try:
        samples = pd.read_csv(args.samples, sep="\t", dtype=str, keep_default_na=False,
                              comment="#")
    except Exception as exc:  # noqa: BLE001 - an unreadable sample sheet is fatal
        raise SystemExit(f"[combine_proteomes] FATAL: cannot read --samples {args.samples}: {exc}")
    if "genome_id" not in samples.columns:
        raise SystemExit(f"[combine_proteomes] FATAL: {args.samples} has no genome_id column")
    genomes = [g.strip() for g in samples["genome_id"] if g.strip()]
    if len(set(genomes)) != len(genomes):
        raise SystemExit(f"[combine_proteomes] FATAL: duplicated genome_id in {args.samples}")

    seen = {}
    n = 0
    os.makedirs(os.path.dirname(os.path.abspath(args.out_faa)), exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_map)), exist_ok=True)
    with open(args.out_faa, "w") as faa, open(args.out_map, "w") as mp:
        mp.write(MAP_HEADER + "\n")
        for gid in genomes:
            path = os.path.join(args.genes_dir, gid, f"{gid}.{args.suffix}")
            if not os.path.exists(path):
                raise SystemExit(f"[combine_proteomes] FATAL: missing {path}")
            for orig_id, seq in utils.read_fasta(path):
                pid = utils.make_protein_id(gid, orig_id)
                if pid in seen:
                    raise SystemExit(f"[combine_proteomes] FATAL: duplicate id {pid} "
                                     f"({seen[pid]} and {gid}:{orig_id})")
                seen[pid] = f"{gid}:{orig_id}"
                faa.write(f">{pid}\n")
                for i in range(0, len(seq), 60):
                    faa.write(seq[i:i + 60] + "\n")
                mp.write(f"{pid}\t{gid}\t{orig_id}\n")
                n += 1
    utils.log(f"Wrote {n} sequences from {len(genomes)} genomes -> {args.out_faa}")


def reference_mode(args) -> None:
    groups = []
    for spec in args.ref:
        if "=" not in spec:
            raise SystemExit(f"[combine_proteomes] FATAL: --ref expects GROUP=PATH, got {spec}")
        grp, path = spec.split("=", 1)
        if not os.path.exists(path):
            raise SystemExit(f"[combine_proteomes] FATAL: reference FASTA {path} missing")
        groups.append((grp.strip(), path))
    cands = list(utils.read_fasta(args.candidates)) if args.candidates else []
    records, labels, mapping, owner = [], [], [], {}
    for cid, seq in cands:
        if utils.sanitize(cid) != cid:
            raise SystemExit(f"[combine_proteomes] FATAL: candidate id '{cid}' is not "
                             "newick-safe (upstream id scheme broken)")
        if cid in owner:
            raise SystemExit(f"[combine_proteomes] FATAL: duplicate candidate id {cid}")
        owner[cid] = "candidate"
        records.append((cid, seq))
    counts = {}
    for grp, path in groups:
        counts[grp] = 0
        for orig, seq in utils.read_fasta(path):
            tid = utils.sanitize(orig)
            if tid in owner:
                raise SystemExit(f"[combine_proteomes] FATAL: reference id '{orig}' -> '{tid}' "
                                 f"collides with an existing {owner[tid]} id; rename it")
            owner[tid] = grp
            records.append((tid, seq))
            labels.append(f"{tid}\t{grp}")
            mapping.append(f"{tid}\t{grp}\t{orig}")
            counts[grp] += 1
    if cands and any(c == 0 for c in counts.values()):
        empty = [g for g, c in counts.items() if c == 0]
        raise SystemExit(f"[combine_proteomes] FATAL: reference set(s) {empty} are empty but "
                         f"{len(cands)} candidate(s) need placement (refs.prometheus_faa / "
                         "refs.paralog_faa are required inputs)")
    utils.write_fasta(records, args.out_faa)
    with open(args.out_labels, "w") as fh:
        fh.write("".join(l + "\n" for l in labels))
    if args.out_map:
        with open(args.out_map, "w") as fh:
            fh.write(REF_MAP_HEADER + "\n" + "".join(m + "\n" for m in mapping))
    utils.log(f"validation-tree input: {len(cands)} candidates + references {counts} "
              f"-> {args.out_faa}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="")
    ap.add_argument("--genes-dir", default="")
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-map", default="")
    ap.add_argument("--suffix", default="faa")
    ap.add_argument("--ref", action="append", default=[],
                    help="reference mode: GROUP=FASTA (repeatable)")
    ap.add_argument("--candidates", default="", help="reference mode: candidate FASTA")
    ap.add_argument("--out-labels", default="", help="reference mode: headerless id<TAB>group")
    args = ap.parse_args()
    if args.ref or args.out_labels:
        if not args.out_labels:
            raise SystemExit("[combine_proteomes] FATAL: reference mode needs --out-labels")
        reference_mode(args)
    else:
        if not args.out_map:
            raise SystemExit("[combine_proteomes] FATAL: proteome mode needs --out-map")
        proteome_mode(args)


if __name__ == "__main__":
    main()
