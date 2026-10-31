"""Codon (in-frame) nucleotide alignment from a protein alignment + CDS via PAL2NAL.

Reads a protein alignment and one or more CDS FASTA files (same combined ids,
``genome__orig``). Every aligned protein is checked against its CDS BEFORE PAL2NAL
runs, the consistent subset is back-translated with ``pal2nal.pl -output fasta
-codontable <T>``, and the PAL2NAL output is verified record by record.

Why the pre-check (PAL2NAL v14.1 behaviour, verified by running it):
  * ONE inconsistent pep/CDS pair makes PAL2NAL abort with an "ERROR: inconsistency"
    block -- which contains '>' header lines -- and exit 0. On the command line v14.1
    prints the block to STDERR and writes NOTHING to stdout (pal2nal.pl:564; the
    STDOUT variant is its HTML/web mode), so the pristine caller silently wrote an
    EMPTY codon alignment for the whole family. A wrapper that merges the streams
    would instead see the block as two FASTA records. Either way it must never reach
    HyPhy: the output ids/order are verified against the input below.
  * The default codon table (1) does not accept GTG as an initiator Met, so every
    GTG-start gene (Prodigal translates the start codon as M) is "inconsistent";
    archaeal/bacterial CDS need ``-codontable 11`` (initiators ATG/GTG/TTG/CTG...).
  * A terminal stop codon in the CDS is tolerated (and dropped) when the protein
    has no trailing '*'; ``-nogap`` deletes every column with a gap in ANY sequence,
    which destroys the 1:1 codon-column <-> protein-column correspondence the
    site->domain mapping relies on, so it is NOT used.

Records are classified (``--report``): kept | no_cds | cds_not_codon_multiple |
length_mismatch | translation_mismatch | internal_stop | protein_mismatch.
Records listed in ``--proteins`` (the loci the caller expects to be codon-aligned;
several files may be given) that cannot be used count toward ``--max-drop-frac``:
above it the script FAILS (an id-scheme or CDS bug must not silently shrink the
selection data set). Aligned records not listed in ``--proteins`` (e.g. reference
sequences without CDS in a cross-domain alignment) are kept when their CDS is
consistent and otherwise dropped with a loud warning.

Exit status is non-zero on: malformed/unaligned input, PAL2NAL missing or failing,
PAL2NAL output inconsistent with the verified input, or too many expected loci
dropped. An empty codon alignment (0 bytes) is written only for a genuinely empty
input (no aligned records) or when fewer than ``--min-seqs`` consistent records exist.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

GAPS = set("-.")
REPORT_HEADER = "id\tstatus\tin_proteins\tn_aa\tn_codons\tnote"

# NCBI translation tables (bases TCAG order, standard NCBI string encoding).
_BASES = "TCAG"
_CODONS = [a + b + c for a in _BASES for b in _BASES for c in _BASES]
_NCBI = {
    1: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "---M------**--*----M---------------M----------------------------"),
    4: ("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "--MM------**-------M------------MMMM---------------M------------"),
    11: ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
         "---M------**--*----M------------MMMM---------------M------------"),
}


def codon_table(n: int):
    """(codon->aa dict, set of initiator codons, set of stop codons) for NCBI table n."""
    if n not in _NCBI:
        raise SystemExit(f"[codon_align] unsupported codon table {n} (supported: {sorted(_NCBI)})")
    aas, starts = _NCBI[n]
    table = dict(zip(_CODONS, aas))
    return table, {c for c, s in zip(_CODONS, starts) if s == "M"}, {c for c, a in table.items() if a == "*"}


def translate(cds: str, table: dict) -> str:
    """Codon-by-codon translation; codons with non-ACGT bases translate to X."""
    return "".join(table.get(cds[i:i + 3], "X") for i in range(0, len(cds) - len(cds) % 3, 3))


def ungap(seq: str) -> str:
    return "".join(c for c in seq if c not in GAPS)


def check_record(prot_aln_seq: str, cds: str, table: dict, starts: set, stops: set):
    """Validate one aligned protein against its CDS.

    Returns (status, cds_used, note). cds_used is the CDS trimmed of a terminal stop
    codon, exactly 3 x (number of residues) long, when status == 'kept'.
    """
    prot = ungap(prot_aln_seq).upper()
    if prot.endswith("*"):
        prot = prot[:-1]
    if "*" in prot:
        return "internal_stop", "", "protein has an internal stop ('*')"
    n_aa = len(prot)
    cds = cds.upper().replace("U", "T")
    if len(cds) % 3:
        return "cds_not_codon_multiple", "", f"CDS {len(cds)} nt is not a whole number of codons"
    n_cod = len(cds) // 3
    if n_cod == n_aa + 1 and cds[-3:] in stops:
        cds = cds[:-3]
        n_cod -= 1
    if n_cod != n_aa:
        return "length_mismatch", "", f"{n_aa} aa vs {n_cod} codons (after removing a terminal stop)"
    tr = translate(cds, table)
    bad = []
    for i, (p, t) in enumerate(zip(prot, tr)):
        if t == "*":                       # an in-frame stop must never reach HyPhy
            return "internal_stop", "", f"in-frame stop codon at codon {i + 1}"
        if p == t or p == "X" or t == "X":
            continue
        if i == 0 and p == "M" and cds[:3] in starts:      # alternative initiator (GTG/TTG...)
            continue
        bad.append(i + 1)
    if bad:
        return "translation_mismatch", "", f"{len(bad)} residue(s) differ from the CDS translation (first at {bad[0]})"
    return "kept", cds, ""


def read_cds_subset(paths, wanted):
    """Stream CDS FASTA file(s), keeping only `wanted` ids (memory stays O(family))."""
    out = {}
    for p in paths:
        if not os.path.exists(p):
            raise SystemExit(f"[codon_align] CDS file not found: {p}")
        for name, seq in utils.read_fasta(p):
            if name not in wanted:
                continue
            seq = seq.replace(" ", "").upper()
            if name in out and out[name] != seq:
                raise SystemExit(f"[codon_align] CDS id {name} occurs twice with different sequences")
            out[name] = seq
    return out


def read_proteins(paths):
    out = {}
    for p in paths:
        if not os.path.exists(p):
            raise SystemExit(f"[codon_align] --proteins file not found: {p}")
        for name, seq in utils.read_fasta(p):
            # same cleaning as family_align: terminal '*' dropped, internal '*' -> X
            s = seq.upper().rstrip("*").replace("*", "X")
            if name in out and out[name] != s:
                raise SystemExit(f"[codon_align] protein id {name} occurs twice with different sequences")
            out[name] = s
    return out


def run_pal2nal(pal2nal, kept, aln_len, codon_tab, workdir):
    """Run PAL2NAL on the verified subset and verify its output. Returns [(id, codon_seq)]."""
    prot_path = os.path.join(workdir, "pal2nal_input.aln")
    cds_path = os.path.join(workdir, "pal2nal_input.fna")
    utils.write_fasta([(n, a) for n, a, _ in kept], prot_path)
    utils.write_fasta([(n, c) for n, _, c in kept], cds_path)
    cmd = ["perl", pal2nal] if pal2nal.endswith(".pl") else [pal2nal]
    cmd += [prot_path, cds_path, "-output", "fasta", "-codontable", str(codon_tab)]
    utils.log("running: " + " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.stderr.strip():
        utils.log("pal2nal stderr: " + proc.stderr.strip()[:2000])
    if proc.returncode != 0:
        raise SystemExit(f"[codon_align] pal2nal exited {proc.returncode}: {proc.stderr[-1000:]}")
    if "ERROR" in proc.stdout or "inconsistency" in proc.stdout:
        raise SystemExit("[codon_align] pal2nal reported an inconsistency despite the pre-check:\n"
                         + proc.stdout[:2000])
    recs, name, buf = [], None, []
    for line in proc.stdout.splitlines():
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    exp_ids = [n for n, _, _ in kept]
    if [n for n, _ in recs] != exp_ids:
        raise SystemExit(f"[codon_align] pal2nal output ids/order differ from the input "
                         f"({len(recs)} records vs {len(exp_ids)} expected)")
    for (n, s), (_, _, cds) in zip(recs, kept):
        if len(s) != 3 * aln_len:
            raise SystemExit(f"[codon_align] pal2nal record {n} has {len(s)} nt, expected {3 * aln_len}")
        if ungap(s).upper() != cds:
            raise SystemExit(f"[codon_align] pal2nal record {n} does not reproduce its verified CDS")
    return recs


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--prot-aln", required=True, help="protein alignment (FASTA)")
    ap.add_argument("--cds", required=True, nargs="+", action="extend",
                    help="CDS FASTA file(s); ids must equal the alignment ids")
    ap.add_argument("--proteins", required=True, nargs="+", action="extend",
                    help="protein FASTA file(s) listing the loci that MUST be codon-aligned "
                         "(repeat the flag or give several files)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--codon-table", type=int, default=11,
                    help="NCBI translation table (11 = bacterial/archaeal; default)")
    ap.add_argument("--max-drop-frac", type=float, default=0.05,
                    help="fail when more than this fraction of the --proteins loci present "
                         "in the alignment cannot be codon-aligned (default 0.05)")
    ap.add_argument("--min-seqs", type=int, default=3)
    ap.add_argument("--report", default=None, help="optional per-record status TSV")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    table, starts, stops = codon_table(args.codon_table)

    def write_report(rows):
        if args.report:
            os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
            with open(args.report, "w") as fh:
                fh.write(REPORT_HEADER + "\n")
                for r in rows:
                    fh.write("\t".join(str(x) for x in r) + "\n")

    prot_records = list(utils.read_fasta(args.prot_aln))
    if not prot_records:
        utils.log("protein alignment is empty; writing an empty codon alignment")
        open(args.out, "w").close()  # genuinely empty input
        write_report([])
        return
    names = [n for n, _ in prot_records]
    if len(set(names)) != len(names):
        raise SystemExit("[codon_align] duplicated ids in the protein alignment")
    lengths = {len(s) for _, s in prot_records}
    if len(lengths) != 1:
        raise SystemExit(f"[codon_align] {args.prot_aln} is not an alignment (record lengths {sorted(lengths)[:5]})")
    aln_len = lengths.pop()

    proteins = read_proteins(args.proteins)
    expected = set(proteins)
    if expected and not expected & set(names):
        raise SystemExit("[codon_align] no aligned record matches an id in --proteins "
                         "(id scheme mismatch between the alignment and the protein FASTA?)")
    cds = read_cds_subset(args.cds, set(names))

    rows, kept = [], []
    for name, aln_seq in prot_records:
        in_p = name in expected
        n_aa = len(ungap(aln_seq).rstrip("*"))
        if in_p and ungap(aln_seq).upper().rstrip("*").replace("*", "X") != proteins[name]:
            rows.append((name, "protein_mismatch", 1, n_aa, "", "aligned sequence != --proteins sequence"))
            continue
        if name not in cds:
            rows.append((name, "no_cds", int(in_p), n_aa, "", "no CDS with this id"))
            continue
        status, cds_used, note = check_record(aln_seq, cds[name], table, starts, stops)
        rows.append((name, status, int(in_p), n_aa, len(cds[name]) // 3, note))
        if status == "kept":
            kept.append((name, aln_seq.upper().replace(".", "-").replace("*", "-"), cds_used))

    dropped = [r for r in rows if r[1] != "kept"]
    exp_in_aln = [r for r in rows if r[2] == 1]
    exp_dropped = [r for r in exp_in_aln if r[1] != "kept"]
    other_dropped = [r for r in dropped if r[2] == 0]
    for r in dropped:
        utils.log(f"drop {r[0]}: {r[1]} ({r[5]})")
    if other_dropped:
        utils.log(f"WARNING: {len(other_dropped)} aligned record(s) NOT listed in --proteins were "
                  f"dropped (e.g. reference sequences without CDS): "
                  f"{', '.join(r[0] for r in other_dropped[:10])}")
    write_report(rows)
    if exp_in_aln:
        frac = len(exp_dropped) / len(exp_in_aln)
        if frac > args.max_drop_frac:
            raise SystemExit(
                f"[codon_align] {len(exp_dropped)}/{len(exp_in_aln)} expected loci could not be "
                f"codon-aligned ({frac:.1%} > --max-drop-frac {args.max_drop_frac:.1%}); e.g. "
                + "; ".join(f"{r[0]}: {r[1]}" for r in exp_dropped[:10]))
    if len(kept) < args.min_seqs:
        utils.log(f"only {len(kept)} consistent record(s) (< {args.min_seqs}); writing an empty "
                  "codon alignment")
        open(args.out, "w").close()  # nonfatal: genuinely too few sequences to align
        return

    pal2nal = shutil.which("pal2nal.pl") or shutil.which("pal2nal")
    if pal2nal is None:
        raise SystemExit("[codon_align] pal2nal.pl not found on PATH")
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.out))) as td:
        recs = run_pal2nal(pal2nal, kept, aln_len, args.codon_table, td)
        tmp_out = os.path.join(td, "codon.fasta")
        utils.write_fasta(recs, tmp_out)
        os.replace(tmp_out, args.out)
    utils.log(f"codon alignment: {len(recs)} records x {aln_len} codons -> {args.out} "
              f"(dropped {len(dropped)}: {len(exp_dropped)} expected, {len(other_dropped)} other)")


if __name__ == "__main__":
    main()
