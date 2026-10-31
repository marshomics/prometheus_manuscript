"""R-loop-forming-sequence (RLFS) prediction in viral sequences -- a faithful
re-implementation of QmRLFS-finder v2.0.

Model (sequence-feature prediction, NOT experimental R-loop mapping)
--------------------------------------------------------------------
QmRLFS (Wongsurawat et al., Nucleic Acids Res. 2012, 40:e16; Jenjaroenpun et al.,
Nucleic Acids Res. 2015, 43:W527-W534, doi:10.1093/nar/gkv344; reference code
github.com/piroonj/QmRLFS-finder, QmRLFS-finder.py v2.0, sha256 253b91bb...). This is the
predictor van den Berg et al. (Cell Host Microbe 2024, 32:1427) applied to phage genomes:
phages restricted by Prometheus carried predicted RLFS, unaffected phages had none.
An RLFS is RIZ + linker + REZ on one strand:

  RIZ  (R-loop initiation zone), two published models, matched with the reference regexes
       m1: G{3,}[ATCGU]{1,10}?G{3,}(?:[ATCGU]{1,10}?G{3,}){1,}?   (>=3 G3-runs)
       m2: G{4,}(?:[ATCGU]{1,10}?G{4,}){1,}?                        (>=2 G4-runs)
       (non-overlapping regex matches, left to right) with G content >= 50 %.
  linker 0-49 nt between the RIZ end and the REZ start (reference: range(50)).
  REZ  (R-loop elongation zone): the LONGEST window of 100..2000 nt starting at one of the
       linker offsets whose G content is >= 40 %, ending strictly before the end of the
       sequence; ties go to the smallest linker (the reference loop order).
Percentages are rounded to two decimals before the threshold test, exactly as the
reference `percent_g` does. Each model is searched on the given strand and on its reverse
complement ('-' calls are reported in plus-strand coordinates, as the reference does).
One RIZ yields at most one RLFS; m1 and m2 can report the same locus.

Equivalence with the reference implementation is tested call-for-call (coordinates,
lengths, G counts, %G) in tests/ecology/run_tests.py on real virus genomes (lambda: one
RLFS, SSV1: none, SH1: 25) and on synthetic constructs.

Per-sequence output (--out; one row per FASTA record)
  genome_id, provirus_id, length, n_rlfs, rloop_prone   (columns read by aggregate_viral.py)
  + n_rlfs_loci, n_rlfs_plus, n_rlfs_minus, n_rlfs_m1, n_rlfs_m2, gc, rlfs_loci_per_kb
  n_rlfs       = number of RLFS reported (all models, both strands) = reference row count
  n_rlfs_loci  = RLFS after merging overlapping calls on the same strand (m1/m2 redundancy)
  rloop_prone  = n_rlfs_loci >= --min-rlfs   (default 1: at least one predicted RLFS)
  gc           = (G+C)/(A+C+G+T), N and other symbols excluded ('' for an all-N record)
Per-RLFS calls (--calls, optional): the reference table columns without the sequences.

Parameters (CLI; defaults are the published model)
  --models     m1,m2 (reference default)
  --riz-min-g  minimum G-run of the m1 RIZ model [3]; m2 uses riz_min_g + 1. Only 3
               reproduces QmRLFS (a different value is logged as a non-QmRLFS variant).
  --window     minimum REZ length = REZ search window [100] (reference window_step).
  --require-rez True (QmRLFS: RLFS = RIZ + linker + REZ). False counts RIZs only
               (permissive, NOT QmRLFS; logged).
  --min-rlfs   RLFS loci needed to call a sequence R-loop-prone [1].

Pure Python + numpy (prefix sums make the REZ search O(50 x 1900) per RIZ instead of the
reference's O(50 x 1900 x length)). An empty FASTA (geNomad found no virus) -> header-only
outputs; duplicated record ids or a malformed FASTA -> error.
"""
import argparse
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

HEADER = ("genome_id\tprovirus_id\tlength\tn_rlfs\trloop_prone\tn_rlfs_loci\tn_rlfs_plus\t"
          "n_rlfs_minus\tn_rlfs_m1\tn_rlfs_m2\tgc\trlfs_loci_per_kb")
CALLS_HEADER = ("genome_id\tprovirus_id\tmodel\tlocation\tstart_RIZ\tend_RIZ\tlength_RIZ\tG_RIZ\t"
                "3Gs_RIZ\t4Gs_RIZ\tperc_G_RIZ\tLinker\tstart_REZ\tend_REZ\tlength_REZ\tG_REZ\t"
                "3Gs_REZ\t4Gs_REZ\tperc_G_REZ\tstrand")

# ---- QmRLFS-finder v2.0 parameters (reference `params` dict) --------------------------
MIN_PERC_G_RIZ = 50       # RIZ G content (%)
NUM_LINKER = 50           # linker offsets 0..49
WINDOW_STEP = 100         # minimum REZ length / search window
MAX_LENGTH_REZ = 2000     # maximum REZ length
MIN_PERC_G_REZ = 40       # REZ G content (%)
QMRLFS_RIZ_MIN_G = 3      # m1 G-run (m2 = 4) in the published models

# IUPAC complement (as Bio.Seq.reverse_complement for upper-case DNA; U pairs like T)
_COMP = str.maketrans("ACGTUNRYKMSWBDHV", "TGCAANYRMKSWVHDB")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def riz_regex(model: str, min_g: int = QMRLFS_RIZ_MIN_G):
    """The reference RIZ regexes (min_g = 3 reproduces QmRLFS exactly)."""
    if model == "m1":
        g = f"G{{{min_g},}}"
        return re.compile(rf"{g}[ATCGU]{{1,10}}?{g}(?:[ATCGU]{{1,10}}?{g}){{1,}}?")
    if model == "m2":
        g = f"G{{{min_g + 1},}}"
        return re.compile(rf"{g}(?:[ATCGU]{{1,10}}?{g}){{1,}}?")
    raise ValueError(f"unknown QmRLFS model {model!r} (use m1 and/or m2)")


def percent_g(seq: str) -> float:
    """Reference percent_g: G % rounded to two decimals (string-formatted, as QmRLFS)."""
    return float("%.2f" % ((seq.count("G") / float(len(seq))) * 100))


def _passes(pct: np.ndarray, threshold: float) -> np.ndarray:
    """pct >= threshold after the reference's two-decimal rounding ("%.2f")."""
    ok = pct >= threshold
    border = np.nonzero((pct >= threshold - 0.01) & (pct < threshold))[0]
    for k in border:            # only values that can round up to the threshold
        if float("%.2f" % pct[k]) >= threshold:
            ok[k] = True
    return ok


def find_rez(prefix_g: np.ndarray, riz_end: int, seqlength: int, window: int = WINDOW_STEP):
    """Reference rez_search (non-quick mode) on the strand sequence with G prefix sums.

    Returns (start, end) of the REZ in strand coordinates, or None. Reference semantics:
    offsets i in 0..49 (i < downstream length), ends j in [i+window, i+2000] with
    riz_end + j < seqlength; keep the longest window with %G >= 40; on equal length the
    first (smallest i) wins.
    """
    n_down = seqlength - riz_end
    best_len, best = 0, None
    for i in range(min(NUM_LINKER, n_down)):
        j_lo = i + window
        j_hi = min(i + MAX_LENGTH_REZ, n_down - 1)       # riz_end + j < seqlength
        if j_hi < j_lo or (j_hi - i) <= best_len:
            continue                                     # cannot beat the current best
        js = np.arange(j_lo, j_hi + 1)
        g = prefix_g[riz_end + js] - prefix_g[riz_end + i]
        lens = js - i
        pct = (g / lens) * 100.0                         # same IEEE ops as the reference
        ok = np.nonzero(_passes(pct, MIN_PERC_G_REZ))[0]
        if ok.size:
            k = ok[-1]                                   # largest j for this offset
            if lens[k] > best_len:
                best_len, best = int(lens[k]), (riz_end + i, riz_end + int(js[k]))
    return best


def search_strand(seq: str, model: str, rx, strand: str, require_rez: bool, window: int):
    """All RLFS of one model on one strand sequence (already upper-case)."""
    n = len(seq)
    arr = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    prefix_g = np.concatenate([[0], np.cumsum(arr == ord("G"), dtype=np.int64)])
    calls = []
    for m in rx.finditer(seq):
        riz_s, riz_e = m.start(), m.end()
        riz_seq = m.group(0)
        riz_pg = percent_g(riz_seq)
        if riz_pg < MIN_PERC_G_RIZ:
            continue
        rez = find_rez(prefix_g, riz_e, n, window) if require_rez else None
        if require_rez and rez is None:
            continue
        rec = {"model": model, "strand": strand, "riz": (riz_s, riz_e), "riz_seq": riz_seq,
               "riz_pg": riz_pg, "rez": rez}
        if rez is not None:
            rec["rez_seq"] = seq[rez[0]:rez[1]]
        calls.append(rec)
    return calls


def to_plus(call, n):
    """Convert a call to plus-strand coordinates (reference output_rlfs)."""
    rs, re_ = call["riz"]
    if call["rez"] is None:
        zs, ze = None, None
    else:
        zs, ze = call["rez"]
    linker = (zs - re_) if zs is not None else None
    if call["strand"] == "-":
        rs, re_ = n - re_, n - rs
        if zs is not None:
            zs, ze = n - ze, n - zs
    if zs is None:
        start, end = rs, re_
    elif call["strand"] == "+":
        start, end = rs, ze
    else:
        start, end = zs, re_
    return {"riz": (rs, re_), "rez": (zs, ze) if zs is not None else None,
            "start": start, "end": end, "linker": linker}


def merge_loci(intervals):
    """Number of connected components of overlapping [start, end) intervals."""
    n, cur_end = 0, None
    for s, e in sorted(intervals):
        if cur_end is None or s >= cur_end:
            n += 1
            cur_end = e
        else:
            cur_end = max(cur_end, e)
    return n


def gc_fraction(seq: str):
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    if acgt == 0:
        return None
    return (s.count("G") + s.count("C")) / acgt


def predict(seq: str, models=("m1", "m2"), riz_min_g=QMRLFS_RIZ_MIN_G, window=WINDOW_STEP,
            require_rez=True):
    """All RLFS of a sequence in reference output order (model, then '+' then '-')."""
    plus = seq.upper()
    minus = revcomp(plus)
    n = len(plus)
    out = []
    for model in models:
        rx = riz_regex(model, riz_min_g)
        for strand, s in (("+", plus), ("-", minus)):
            for c in search_strand(s, model, rx, strand, require_rez, window):
                c.update(to_plus(c, n))
                out.append(c)
    return out


def summarize(calls, length, min_rlfs=1):
    by_strand = {"+": [], "-": []}
    for c in calls:
        by_strand[c["strand"]].append((c["start"], c["end"]))
    loci = merge_loci(by_strand["+"]) + merge_loci(by_strand["-"])
    return {
        "n_rlfs": len(calls),
        "n_rlfs_loci": loci,
        "n_rlfs_plus": len(by_strand["+"]),
        "n_rlfs_minus": len(by_strand["-"]),
        "n_rlfs_m1": sum(c["model"] == "m1" for c in calls),
        "n_rlfs_m2": sum(c["model"] == "m2" for c in calls),
        "rloop_prone": loci >= min_rlfs,
        "rlfs_loci_per_kb": (loci / (length / 1000.0)) if length else 0.0,
    }


def _call_row(genome, name, c):
    rs, re_ = c["riz"]
    riz_seq = c["riz_seq"]
    cols = [genome, name, c["model"], f"{name}:{c['start'] + 1}-{c['end']}", rs, re_, re_ - rs,
            riz_seq.count("G"), riz_seq.count("GGG"), riz_seq.count("GGGG"), c["riz_pg"]]
    if c["rez"] is None:
        cols += [""] * 8
    else:
        zs, ze = c["rez"]
        rez_seq = c["rez_seq"]
        cols += [c["linker"], zs, ze, ze - zs, rez_seq.count("G"), rez_seq.count("GGG"),
                 rez_seq.count("GGGG"), percent_g(rez_seq)]
    cols.append(c["strand"])
    return "\t".join(map(str, cols))


def _parse_bool(v: str) -> bool:
    s = str(v).strip().lower()
    if s in ("true", "1", "yes", "t"):
        return True
    if s in ("false", "0", "no", "f"):
        return False
    raise SystemExit(f"[rloop_predict] --require-rez must be True/False, got {v!r}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--riz-min-g", type=int, default=QMRLFS_RIZ_MIN_G)
    ap.add_argument("--window", type=int, default=WINDOW_STEP)
    ap.add_argument("--require-rez", default="True")
    ap.add_argument("--models", default="m1,m2")
    ap.add_argument("--min-rlfs", type=int, default=1)
    ap.add_argument("--calls", default="", help="optional per-RLFS table")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    require_rez = _parse_bool(args.require_rez)
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    for m in models:
        riz_regex(m)                                    # validates the model names
    if not models:
        raise SystemExit("[rloop_predict] --models is empty")
    if args.riz_min_g < 1 or args.window < 1 or args.min_rlfs < 1:
        raise SystemExit("[rloop_predict] --riz-min-g, --window and --min-rlfs must be >= 1")
    variant = []
    if args.riz_min_g != QMRLFS_RIZ_MIN_G:
        variant.append(f"riz_min_g={args.riz_min_g} (QmRLFS: 3)")
    if args.window != WINDOW_STEP:
        variant.append(f"window={args.window} (QmRLFS: 100)")
    if not require_rez:
        variant.append("require_rez=False (RIZ-only, not QmRLFS)")
    if variant:
        utils.log("WARNING: non-QmRLFS parameter(s): " + "; ".join(variant))

    if not os.path.exists(args.fasta):
        raise SystemExit(f"[rloop_predict] FASTA not found: {args.fasta}")
    # whitespace (e.g. '\r' of CRLF files) is not sequence: it would break the regex RIZ
    records = [(n, re.sub(r"\s+", "", s)) for n, s in utils.read_fasta(args.fasta)]
    if os.path.getsize(args.fasta) > 0 and not records:
        raise SystemExit(f"[rloop_predict] {args.fasta} is not empty but holds no FASTA record")
    names = [n for n, _ in records]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise SystemExit(f"[rloop_predict] duplicated record id(s) in {args.fasta}: {dup[:5]}")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    calls_fh = None
    if args.calls:
        os.makedirs(os.path.dirname(os.path.abspath(args.calls)), exist_ok=True)
        calls_fh = open(args.calls, "w")
        calls_fh.write(CALLS_HEADER + "\n")
    n_prone = 0
    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for name, seq in records:
            calls = predict(seq, models, args.riz_min_g, args.window, require_rez)
            s = summarize(calls, len(seq), args.min_rlfs)
            n_prone += s["rloop_prone"]
            gc = gc_fraction(seq)
            out.write("\t".join(map(str, [
                args.genome, name, len(seq), s["n_rlfs"], s["rloop_prone"], s["n_rlfs_loci"],
                s["n_rlfs_plus"], s["n_rlfs_minus"], s["n_rlfs_m1"], s["n_rlfs_m2"],
                "" if gc is None else round(gc, 5), round(s["rlfs_loci_per_kb"], 5)])) + "\n")
            if calls_fh:
                for c in calls:
                    calls_fh.write(_call_row(args.genome, name, c) + "\n")
    if calls_fh:
        calls_fh.close()
    utils.log(f"{args.genome}: {len(records)} viral sequence(s), {n_prone} R-loop-prone "
              f"(models={','.join(models)}, riz_min_g={args.riz_min_g}, window={args.window}, "
              f"require_rez={require_rez}, min_rlfs={args.min_rlfs})")


if __name__ == "__main__":
    main()
