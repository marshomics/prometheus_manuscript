"""Regenerate the synthetic R-loop test sequences (deterministic, seeded).

    python tests/ecology/make_rloop_fixtures.py

Writes two FASTA files under tests/ecology/fixtures/:
  rloop_synthetic.fna           constructs with HAND-DERIVED expected QmRLFS answers
                                (asserted in run_tests.py: test_rloop_constructs_known_answers)
  rloop_reference_input_v1.fna  the exact input on which the reference QmRLFS-finder v2.0
                                was run to produce qmrlfs_ref_synthetic_v1.table.tsv
                                (asserted call-for-call in test_rloop_matches_reference_qmrlfs)
Both files share the seeded A/T backgrounds and the random sequences; they differ only in
the hand-built constructs (v1 lacks the 'T' spacer after the RIZ, so its RIZ absorbs the
first G of the REZ -- the reference output records exactly that).
"""
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
FIX = os.path.join(HERE, "fixtures")

_COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(_COMP)[::-1]


def at(n, rng):
    """A/T-only background: no G on either strand, so it can never seed a RIZ/REZ."""
    return "".join(rng.choice("AT") for _ in range(n))


def build(version=2):
    rng = random.Random(20240814)
    bg5, bg3 = at(300, rng), at(300, rng)
    rez = "GA" * 100                         # 200 nt, 50 % G
    sp = "T" if version >= 2 else ""         # spacer that keeps the REZ's G out of the RIZ
    seqs = {}
    # m1 RIZ (3 G3-runs, 9/11 = 81.82 % G) + spacer + 50 %-G REZ: plus strand only.
    # v2 expected (QmRLFS): RIZ [300,311), linker 0, REZ [311,561) = 250 nt, 100 G (40.0 %):
    # the REZ is the LONGEST >=40 %-G window, so it runs 50 nt into the A/T tail.
    seqs["syn_plus_m1"] = bg5 + "GGGTGGGTGGG" + sp + rez + bg3
    # m2 RIZ (2 G4-runs) + spacer + REZ, written as its reverse complement -> minus strand only
    seqs["syn_minus_m2"] = rc(bg5 + "GGGGAGGGG" + sp + rez + bg3)
    # RIZ but no G-rich REZ downstream -> no RLFS (RIZ-only mode: 1)
    seqs["syn_riz_no_rez"] = bg5 + "GGGTGGGTGGG" + bg3
    if version >= 2:
        # 100 nt after the RIZ: a 100-nt REZ would end exactly at the sequence end, which
        # QmRLFS rejects (REZ end must be < sequence length) -> no RLFS
        seqs["syn_rez_at_end"] = bg5 + "GGGTGGGTGGG" + "T" + "GA" * 49 + "G"
        # ... one extra trailing base makes the 100-nt REZ [311,411) valid -> 1 RLFS
        seqs["syn_rez_at_end_plus1"] = bg5 + "GGGTGGGTGGG" + "T" + "GA" * 49 + "G" + "A"
    else:
        seqs["syn_rez_at_end"] = bg5 + "GGGTGGGTGGG" + "GA" * 50
        seqs["syn_rez_at_end_plus1"] = bg5 + "GGGTGGGTGGG" + "GA" * 50 + "A"
    # m1 regex match whose G content is 9/27 = 33 % < 50 % -> not a RIZ
    seqs["syn_low_g_riz"] = bg5 + "GGGAAAAAAAAAGGGAAAAAAAAAGGG" + sp + rez + bg3
    # an N inside the RIZ breaks the [ATCGU] linker -> no RIZ
    seqs["syn_n_break"] = bg5 + "GGGNGGGTGGG" + sp + rez + bg3
    # soft-masked (lower-case) copy of syn_plus_m1 -> same call as upper case
    seqs["syn_lowercase"] = seqs["syn_plus_m1"].lower()
    # random 50 % GC sequences and a G-skewed random sequence (many overlapping calls)
    for k in range(3):
        seqs[f"syn_random_gc50_{k + 1}"] = "".join(rng.choice("ACGT") for _ in range(4000))
    seqs["syn_random_gskew"] = "".join(rng.choices("ACGT", weights=[17, 18, 47, 18], k=4000))
    return seqs


def write(seqs, path):
    with open(path, "w") as fh:
        for name, s in seqs.items():
            fh.write(f">{name}\n")
            for i in range(0, len(s), 70):
                fh.write(s[i:i + 70] + "\n")


def main():
    os.makedirs(FIX, exist_ok=True)
    write(build(2), os.path.join(FIX, "rloop_synthetic.fna"))
    write(build(1), os.path.join(FIX, "rloop_reference_input_v1.fna"))
    print(FIX)


if __name__ == "__main__":
    main()
