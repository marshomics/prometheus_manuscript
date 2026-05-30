"""Codon (in-frame) nucleotide alignment from a protein alignment + CDS via PAL2NAL.

Reads a protein alignment and the combined CDS FASTA (same combined ids,
``genome__orig``). For every aligned protein record it fetches the matching CDS,
builds an un-gapped nucleotide FASTA subset in the SAME order/names as the
alignment, and runs ``pal2nal.pl`` to back-translate the protein alignment into a
codon (in-frame) nucleotide alignment.

Defensive contract:
  - PAL2NAL needs each CDS length to be approximately 3x the protein length
    (it tolerates a trailing stop). Records whose CDS is absent or whose length
    is incompatible are skipped and logged, never crash the run.
  - If pal2nal.pl is missing on PATH, or fewer than 3 usable records remain,
    an empty output file is written (touch) and the reason is logged. Downstream
    selection steps treat an empty codon alignment as "nothing to test".
"""
import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

# pal2nal needs CDS length within this many *codons* of 3x the (un-gapped)
# protein length; a single trailing stop codon is the usual, tolerated excess.
LEN_TOL_CODONS = 1


def ungapped_len(seq: str) -> int:
    """Residue count of an aligned sequence (strip the usual gap characters)."""
    return sum(1 for c in seq if c not in "-.*")


def cds_compatible(prot_aln_seq: str, cds_seq: str) -> bool:
    """True if a CDS can back-translate the (aligned) protein under PAL2NAL.

    The CDS must be a whole number of codons and span at least the protein's
    residues, allowing up to LEN_TOL_CODONS extra codons (e.g. a stop).
    """
    aa = ungapped_len(prot_aln_seq)
    nt = len(cds_seq)
    if aa == 0 or nt == 0:
        return False
    if nt % 3 != 0:
        return False
    n_codons = nt // 3
    # CDS must cover the residues and not exceed them by more than the tolerance.
    return aa <= n_codons <= aa + LEN_TOL_CODONS


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prot-aln", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    def write_empty_codon(reason: str):
        utils.log(f"Empty codon alignment: {reason}")
        open(args.out, "w").close()

    # --- read inputs ---
    prot_records = list(utils.read_fasta(args.prot_aln))
    if not prot_records:
        write_empty_codon("protein alignment is empty")
        return

    cds = dict(utils.read_fasta(args.cds))
    if not cds:
        write_empty_codon("combined CDS FASTA is empty/missing")
        return

    # --- build a CDS subset in the SAME order/names as the protein alignment ---
    subset = []          # list of (name, cds_seq) for pal2nal
    kept_prot = []       # corresponding protein-alignment records (same order)
    n_missing = 0
    n_mismatch = 0
    for name, aln_seq in prot_records:
        cds_seq = cds.get(name)
        if cds_seq is None:
            n_missing += 1
            utils.log(f"skip {name}: no CDS in combined ffn")
            continue
        cds_seq = cds_seq.upper().replace("U", "T")
        if not cds_compatible(aln_seq, cds_seq):
            n_mismatch += 1
            utils.log(
                f"skip {name}: CDS length {len(cds_seq)} nt incompatible with "
                f"{ungapped_len(aln_seq)} aa (need ~3x, whole codons)"
            )
            continue
        subset.append((name, cds_seq))
        kept_prot.append((name, aln_seq))

    if len(subset) < 3:
        write_empty_codon(
            f"only {len(subset)} usable record(s) "
            f"(missing CDS={n_missing}, length-mismatch={n_mismatch}); need >=3"
        )
        return

    # --- pal2nal availability ---
    pal2nal = shutil.which("pal2nal.pl") or shutil.which("pal2nal")
    if pal2nal is None:
        write_empty_codon("pal2nal.pl not found on PATH")
        return

    # PAL2NAL requires the protein alignment and the CDS to contain the same
    # records in the same order. Write a pruned protein alignment + the matching
    # CDS subset to a working directory next to the output.
    workdir = os.path.dirname(os.path.abspath(args.out))
    prot_subset_path = os.path.join(workdir, "codon_prot_subset.faa")
    cds_subset_path = os.path.join(workdir, "codon_cds_subset.fna")
    utils.write_fasta(kept_prot, prot_subset_path)
    utils.write_fasta(subset, cds_subset_path)

    cmd = [
        pal2nal, prot_subset_path, cds_subset_path,
        "-output", "fasta", "-nogap",
    ]
    utils.log("running: " + " ".join(cmd))
    try:
        with open(args.out, "w") as out_fh:
            proc = subprocess.run(
                cmd, stdout=out_fh, stderr=subprocess.PIPE, check=False
            )
        if proc.stderr:
            utils.log("pal2nal stderr: " + proc.stderr.decode(errors="replace").strip())
        if proc.returncode != 0:
            utils.log(f"pal2nal exited {proc.returncode}; writing empty output")
            open(args.out, "w").close()
            return
    except Exception as e:  # noqa: BLE001 - never crash the DAG
        utils.log(f"pal2nal failed: {e}; writing empty output")
        open(args.out, "w").close()
        return

    # Guard against pal2nal writing nothing useful (header-less / empty).
    n_out = len(utils.fasta_ids(args.out))
    if n_out == 0:
        utils.log("pal2nal produced no records; output is empty")
    else:
        utils.log(
            f"codon alignment: {n_out} records "
            f"(skipped missing CDS={n_missing}, mismatch={n_mismatch})"
        )


if __name__ == "__main__":
    main()
