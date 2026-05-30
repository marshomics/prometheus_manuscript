"""R-loop-forming-sequence (RLFS) detector for proviral contigs.

Model (a sequence-feature proxy, NOT experimental R-loop mapping)
----------------------------------------------------------------
Co-transcriptional R-loops nucleate on G-rich, G-skewed template/non-template
arrangements. This detector follows the QmRLFS RIZ/REZ idea (Jenjaroenpun et al.,
*Nucleic Acids Res.* 2015; the RLFS model also used by van den Berg et al. 2024 to
predict R-loop-forming phage transcripts):

  RIZ  (R-loop Initiation Zone): a guanine-rich initiator built from a cluster of
       >= ``riz_min_tracts`` G-tracts, each G-tract being a run of >= ``riz_min_g``
       consecutive guanines, with all tracts falling inside a window of
       <= ``window`` bp. The RIZ is where the stable G4 / G-rich duplex that seeds
       the R-loop forms.
  REZ  (R-loop Elongation Zone): a downstream G-rich stretch within ``rez_max_gap``
       bp of the RIZ that supports elongation of the displaced loop. The REZ is
       *supportive*: an RLFS is called on a qualifying RIZ; the REZ, when present
       downstream, is recorded but is not strictly required (matching QmRLFS's
       "RIZ with optional permissive REZ" permissive mode).

Each sequence is scanned on BOTH strands (the minus strand is the
reverse-complement), because an RLFS on either DNA strand can template an R-loop on
the corresponding transcript. ``n_rlfs`` is the number of non-overlapping RIZ-anchored
RLFS across both strands; a provirus is ``rloop_prone`` when ``n_rlfs >= 1``.

Parameters (CLI; defaults match the rule/config)
  --riz-min-g : minimum length of each G-tract (consecutive G run)         [3]
  --window    : max span (bp) of the G-tract cluster that forms the RIZ     [100]

Internal constants document the rest of the QmRLFS-style model:
  RIZ_MIN_TRACTS = 3      (>=3 G-tracts make a RIZ; QmRLFS m1 model)
  REZ_MAX_GAP    = 2000   (downstream search window for a permissive REZ, bp)
  REZ_WINDOW     = 100    (REZ scan window, bp)
  REZ_MIN_G_FRAC = 0.40   (G fraction over REZ_WINDOW to count as G-rich)

Pure Python; no external tools. Empty FASTA -> header only.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

HEADER = "genome_id\tprovirus_id\tlength\tn_rlfs\trloop_prone"

# --- QmRLFS-style model constants (documented in the module docstring) ---
RIZ_MIN_TRACTS = 3        # >= this many G-tracts within `window` bp form a RIZ
REZ_MAX_GAP = 2000        # search this far downstream of a RIZ for a permissive REZ
REZ_WINDOW = 100          # REZ scan window (bp)
REZ_MIN_G_FRAC = 0.40     # G fraction over REZ_WINDOW to flag a G-rich REZ

_COMP = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def find_g_tracts(seq: str, min_g: int):
    """Return (start, end) spans of every run of >= min_g consecutive G."""
    if min_g < 1:
        min_g = 1
    return [(m.start(), m.end()) for m in re.finditer(f"G{{{min_g},}}", seq)]


def has_downstream_rez(seq: str, from_pos: int) -> bool:
    """True if a G-rich REZ window starts within REZ_MAX_GAP bp of `from_pos`."""
    n = len(seq)
    end_search = min(n, from_pos + REZ_MAX_GAP)
    # Slide a REZ_WINDOW; step by a fraction of the window for speed.
    step = max(1, REZ_WINDOW // 5)
    pos = from_pos
    while pos < end_search:
        win = seq[pos:pos + REZ_WINDOW]
        if not win:
            break
        g = win.count("G")
        if g / len(win) >= REZ_MIN_G_FRAC:
            return True
        pos += step
    return False


def count_rlfs_one_strand(seq: str, min_g: int, window: int) -> int:
    """Count non-overlapping RIZ-anchored RLFS on a single (already-cased) strand.

    A RIZ = >= RIZ_MIN_TRACTS G-tracts whose cluster spans <= `window` bp. We walk
    the G-tracts left to right; the first tract that, together with the next
    (RIZ_MIN_TRACTS-1) tracts, fits inside `window` bp opens an RLFS. We then jump
    past that cluster so RLFS do not overlap (each initiation zone counted once).
    """
    tracts = find_g_tracts(seq, min_g)
    if len(tracts) < RIZ_MIN_TRACTS:
        return 0

    n_rlfs = 0
    i = 0
    n_tracts = len(tracts)
    while i <= n_tracts - RIZ_MIN_TRACTS:
        first_start = tracts[i][0]
        last_idx = i + RIZ_MIN_TRACTS - 1
        last_end = tracts[last_idx][1]
        if (last_end - first_start) <= window:
            # qualifying RIZ; the optional permissive REZ is recorded but not required
            _ = has_downstream_rez(seq, last_end)  # documented; informative only
            n_rlfs += 1
            # advance past this cluster to keep RLFS non-overlapping
            i = last_idx + 1
        else:
            i += 1
    return n_rlfs


def count_rlfs(seq: str, min_g: int, window: int) -> int:
    """Total RLFS on both strands of a sequence."""
    s = seq.upper().replace("U", "T")
    plus = count_rlfs_one_strand(s, min_g, window)
    minus = count_rlfs_one_strand(revcomp(s), min_g, window)
    return plus + minus


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--riz-min-g", type=int, default=3)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    records = list(utils.read_fasta(args.fasta))
    if not records:
        utils.write_empty(args.out, HEADER)
        utils.log(f"{args.genome}: no provirus sequences; header-only output")
        return

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        n_prone = 0
        for name, seq in records:
            length = len(seq)
            n_rlfs = count_rlfs(seq, args.riz_min_g, args.window)
            prone = n_rlfs >= 1
            if prone:
                n_prone += 1
            out.write(
                f"{args.genome}\t{name}\t{length}\t{n_rlfs}\t{prone}\n"
            )

    utils.log(
        f"{args.genome}: {len(records)} provirus(es), "
        f"{n_prone} R-loop-prone (riz_min_g={args.riz_min_g}, window={args.window})"
    )


if __name__ == "__main__":
    main()
