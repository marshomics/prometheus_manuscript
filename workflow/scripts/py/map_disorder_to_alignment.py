"""Map per-residue consensus disorder onto the Prometheus family alignment.

Whole-protein %disorder across divergent sequences is confounded by length and domain
content; anchoring disorder to the (UNtrimmed) family alignment compares HOMOLOGOUS
positions. For every alignment column that carries a residue of at least one scored
locus we report:

  alignment_col          1-based column of the untrimmed MSA (evolution/align/family.aln)
  n_loci                 scored loci with a residue in that column
  mean_disorder          mean consensus call (0/1) of those loci
  frac_loci_disordered   fraction of those loci called disordered -- identical to
                         mean_disorder by construction (the consensus is binary); kept
                         as a separate SCHEMA column
  domain                 domain category of the column (NA = linker/tail), see below
  occupancy              n_loci / number of scored loci in the alignment   (additive)
  domain_support         fraction of the column's residue-bearing annotated loci whose
                         residue lies in that category                     (additive)
  mean_meta_score        mean continuous metapredict score                 (additive)

Domain track: each aligned locus' own hmmscan hits are turned into per-residue
categories with utils.load_domain_spec / utils.domain_intervals (exact Pfam matching,
Contract v2 item 3); each column takes the plurality category among the aligned loci
that have >=1 categorised hit and a residue there, if that category covers at least
--min-domain-support of them (default 0.5) and is not tied with another category, else
NA (residues claimed by two categories of one locus go to the canonically earlier one).
A single 'reference' locus is not
used: DUF3320 is absent from ~54% of Prometheus loci, so any one reference would
mislabel whole domains for the rest.

Residue -> column mapping walks each aligned row counting non-gap characters; it is
verified per locus (same residue count, same residues after the predictor's B/U/X/Z
conversion) so a stale or re-ordered alignment fails instead of being silently
mis-mapped. Genuinely empty inputs (no scored residues, no alignment, or the <3-sequence
pass-through that family_align writes instead of an MSA) give a header-only output.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

HEADER_COLUMNS = ["alignment_col", "n_loci", "mean_disorder", "frac_loci_disordered",
                  "domain", "occupancy", "domain_support", "mean_meta_score"]
HEADER = "\t".join(HEADER_COLUMNS)
GAP_CHARS = frozenset("-.")
CONVERSION = {"B": "N", "U": "C", "X": "G", "Z": "Q"}   # = predict_disorder.CONVERSION


class MapError(RuntimeError):
    """Inconsistent inputs: the rule must fail."""


def read_residue_table(path):
    """-> {locus: {"cons": [0/1 by position], "meta": [float|None], "aa": str|None}}."""
    tracks = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        need = ("locus_id", "position", "consensus")
        if any(c not in header for c in need):
            raise MapError(f"{path}: residue table lacks one of {need}")
        li, pi, ci = (header.index(c) for c in need)
        mi = header.index("meta_score") if "meta_score" in header else None
        ai = header.index("aa") if "aa" in header else None
        for ln, line in enumerate(fh, 2):
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            try:
                locus, pos, cons = f[li], int(f[pi]), int(float(f[ci]))
            except (IndexError, ValueError):
                raise MapError(f"{path}:{ln}: malformed residue row")
            if cons not in (0, 1):
                raise MapError(f"{path}:{ln}: consensus must be 0/1, got {cons}")
            t = tracks.setdefault(locus, {"cons": [], "meta": [], "aa": []})
            if pos != len(t["cons"]) + 1:
                raise MapError(f"{path}:{ln}: {locus} positions not contiguous from 1 "
                               f"(got {pos} after {len(t['cons'])})")
            t["cons"].append(cons)
            meta = f[mi] if mi is not None and mi < len(f) else ""
            try:
                t["meta"].append(float(meta) if meta not in ("", "NA", "nan") else None)
            except ValueError:
                raise MapError(f"{path}:{ln}: non-numeric meta_score {meta!r}")
            t["aa"].append(f[ai] if ai is not None and ai < len(f) else "")
    for t in tracks.values():
        t["aa"] = "".join(t["aa"]) if all(t["aa"]) else None
    return tracks


def residue_columns(aligned_seq):
    """0-based alignment column of each residue (1-based residue i -> cols[i-1])."""
    return [c for c, ch in enumerate(aligned_seq) if ch not in GAP_CHARS]


def ungapped(aligned_seq):
    return "".join(ch for ch in aligned_seq if ch not in GAP_CHARS)


def locus_residue_categories(hits, spec):
    """1-based residue -> category for one locus via utils.domain_intervals.
    Residues claimed by two categories (overlapping hit ends) go to the category that
    comes first in the canonical N->C order."""
    res = {}
    ivs = utils.domain_intervals(hits, spec)
    for cat in utils.canonical_order(spec):
        for s, e in ivs.get(cat, []):
            for r in range(s, e + 1):
                res.setdefault(r, cat)
    return res


def column_domains(aln, hits, spec, min_support):
    """-> ({col: category}, {col: support}, n_annotated_loci)."""
    order = utils.canonical_order(spec)
    votes, n_res = {}, {}
    n_annot = 0
    for locus, seq in aln:
        lh = hits.get(locus, [])
        rescat = locus_residue_categories(lh, spec) if lh else {}
        if not rescat:
            continue
        n_annot += 1
        for i, col in enumerate(residue_columns(seq), 1):
            n_res[col] = n_res.get(col, 0) + 1
            cat = rescat.get(i)
            if cat is not None:
                votes.setdefault(col, {}).setdefault(cat, 0)
                votes[col][cat] += 1
    dom, sup = {}, {}
    for col, n in n_res.items():
        v = votes.get(col)
        if not v:
            sup[col] = 0.0
            continue
        ranked = sorted(v.items(), key=lambda kv: (-kv[1], order.index(kv[0])))
        best, top = ranked[0]
        sup[col] = top / n
        tied = len(ranked) > 1 and ranked[1][1] == top     # two categories equally supported
        if sup[col] >= min_support and not tied:
            dom[col] = best
    return dom, sup, n_annot


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--residue", required=True, help="disorder_per_residue.tsv")
    ap.add_argument("--alignment", required=True,
                    help="UNtrimmed family alignment (evolution/align/family.aln)")
    ap.add_argument("--domtbl", default="", help="hmmscan --domtblout for the domain track")
    ap.add_argument("--domains", default="", help="prometheus.domains spec as JSON")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-domain-support", type=float, default=0.5,
                    help="min fraction of annotated loci agreeing on a column's domain")
    args = ap.parse_args(argv)
    if not (0.0 < args.min_domain_support <= 1.0):
        raise MapError("--min-domain-support must be in (0, 1]")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    for p in (args.residue, args.alignment):
        if not os.path.isfile(p):
            raise MapError(f"input {p!r} not found")
    tracks = read_residue_table(args.residue)
    aln = list(utils.read_fasta(args.alignment))
    if not tracks or not aln:
        utils.write_empty(args.out, HEADER)
        utils.log(f"{'no scored residues' if not tracks else 'empty alignment'}; "
                  "header-only disorder-by-column table.")
        return 0
    ids = [k for k, _ in aln]
    if len(ids) != len(set(ids)):
        raise MapError(f"{args.alignment}: duplicated sequence ids")
    widths = {len(s) for _, s in aln}
    if len(widths) != 1:
        if len(aln) < 3:
            # family_align copies the raw FASTA when it has < 3 sequences (no MSA).
            utils.write_empty(args.out, HEADER)
            utils.log(f"{args.alignment} holds {len(aln)} unaligned sequence(s) "
                      "(family_align pass-through for < 3 loci); header-only output.")
            return 0
        raise MapError(f"{args.alignment} is not an alignment (row lengths {sorted(widths)})")
    aln_len = widths.pop()

    aln_ids = set(ids)
    missing = sorted(set(tracks) - aln_ids)
    if missing:
        raise MapError(f"{len(missing)} scored locus/loci absent from the alignment "
                       f"(e.g. {missing[:3]}); alignment and disorder inputs disagree")
    for locus, seq in aln:
        t = tracks.get(locus)
        if t is None:
            continue
        raw = ungapped(seq).upper().rstrip("*")
        if len(raw) != len(t["cons"]):
            raise MapError(f"{locus}: {len(raw)} residues in the alignment but "
                           f"{len(t['cons'])} in the disorder table")
        if t["aa"] is not None:
            conv = "".join(CONVERSION.get(c, c) for c in raw)
            if conv != t["aa"]:
                i = next(k for k, (a, b) in enumerate(zip(conv, t["aa"])) if a != b)
                raise MapError(f"{locus}: alignment residue {i + 1} is {conv[i]!r} but the "
                               f"disorder table has {t['aa'][i]!r}")

    dom, sup, n_annot = {}, {}, 0
    if args.domtbl:
        if not os.path.isfile(args.domtbl):
            raise MapError(f"domtbl {args.domtbl!r} not found")
        spec = utils.load_domain_spec(args.domains)
        if not spec:
            raise MapError("--domtbl given but the --domains spec is empty")
        hits = utils.parse_domtbl(args.domtbl)
        dom, sup, n_annot = column_domains(aln, hits, spec, args.min_domain_support)
        if n_annot == 0:
            raise MapError(f"no aligned locus has a hit to any configured domain in "
                           f"{args.domtbl}; the domain track cannot be built")
    else:
        utils.log("WARNING: no --domtbl; every column's domain is NA.")

    n_scored = len(tracks)
    col_n = [0] * aln_len
    col_dis = [0] * aln_len
    col_meta = [0.0] * aln_len
    col_meta_n = [0] * aln_len
    for locus, seq in aln:
        t = tracks.get(locus)
        if t is None:
            continue
        for i, col in enumerate(residue_columns(seq)):
            if i >= len(t["cons"]):        # a trailing '*' kept by the aligner
                break
            col_n[col] += 1
            col_dis[col] += t["cons"][i]
            m = t["meta"][i]
            if m is not None:
                col_meta[col] += m
                col_meta_n[col] += 1

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for col in range(aln_len):
            n = col_n[col]
            if n == 0:
                continue
            frac = round(col_dis[col] / n, 4)
            mm = round(col_meta[col] / col_meta_n[col], 4) if col_meta_n[col] else ""
            ds = round(sup[col], 4) if col in sup else ""
            out.write(f"{col + 1}\t{n}\t{frac}\t{frac}\t{dom.get(col, 'NA')}\t"
                      f"{round(n / n_scored, 4)}\t{ds}\t{mm}\n")
    utils.log(f"disorder-by-column: {sum(1 for n in col_n if n)} columns, {n_scored} scored "
              f"loci, {len(aln)} aligned, {n_annot} with domain hits.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except MapError as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)
