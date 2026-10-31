"""Rewrite Prodigal gene calls to Contract-v2 gene ids and write the per-gene table.

Input: the raw Prodigal 2.6.3 outputs (-a proteins, -d CDS nucleotides, -f gff). Formats
(verified on real Prodigal 2.6.3 runs, see tests/preprocess_qc/REVIEW.md):

  FASTA header  >{seqid}_{n} # {start} # {end} # {1|-1} # ID={k}_{n};partial={xy};start_type=..
  GFF comments  # Sequence Data: seqnum={k};seqlen={L};seqhdr="{full FASTA header}"
                # Model Data: version=Prodigal.v2.6.3;run_type=..;...;transl_table={T};...
  GFF CDS line  seqid  Prodigal_v2.6.3  CDS  start  end  score  +|-  0  ID={k}_{n};partial=..

  k = 1-based index of the record in the input FASTA (every record, also those without a
  gene); n = 1-based index of the gene on that record, in coordinate order. The ID=k_n
  attribute (not the first header token, which repeats when contig names repeat) is the
  key that links the .faa, .fna and .gff records of one gene.

Output (Contract v2): every gene gets id utils.make_gene_id(genome, k, n) =
`{genome}_{k}_{n}`, identical in the .faa, the .fna, the .gff (ID attribute) and
genes.tsv. FASTA headers carry only that id. genes.tsv columns:

  protein_id, contig, contig_index, gene_index, start, end, strand (+/-),
  partial (Prodigal code: 00 complete, 10/01/11 = runs off the left/right/both ends of
  the contig, left = lower coordinate), contig_length, length_aa (residues, excluding the
  terminal '*' stop symbol), transl_table (from the GFF Model Data line)

Every protein must have exactly one CDS record and one GFF feature with identical
coordinates, and the CDS length must be 3 x the protein length; any mismatch, a malformed
header, or zero predicted genes is an error (exit 1) -- never a silently empty proteome.
"""
from __future__ import annotations
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

GENES_COLS = ["protein_id", "contig", "contig_index", "gene_index", "start", "end", "strand",
              "partial", "contig_length", "length_aa", "transl_table"]
HDR_RE = re.compile(r"^(\S+) # (\d+) # (\d+) # (-?1) # (.*)$")
SEQDATA_RE = re.compile(r'^# Sequence Data: seqnum=(\d+);seqlen=(\d+);seqhdr="(.*)"\s*$')
TT_RE = re.compile(r"transl_table=(\d+)")
GFF_ID_RE = re.compile(r"(^|;)ID=([^;]*)")


class ProdigalFormatError(Exception):
    pass


def read_fasta_full(path):
    """[(full_header_without_>, sequence)] in file order."""
    recs, hdr, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if hdr is not None:
                    recs.append((hdr, "".join(seq)))
                hdr, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if hdr is not None:
        recs.append((hdr, "".join(seq)))
    return recs


def parse_attrs(s: str) -> dict:
    out = {}
    for kv in s.strip().strip(";").split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_prodigal_fasta(path, label):
    """{prodigal_id 'k_n': dict(start, end, strand, partial, seq)} preserving file order."""
    genes = {}
    for hdr, seq in read_fasta_full(path):
        m = HDR_RE.match(hdr)
        if not m:
            raise ProdigalFormatError(f"{label}: header is not Prodigal format: {hdr[:120]!r}")
        attrs = parse_attrs(m.group(5))
        pid = attrs.get("ID", "")
        if not re.fullmatch(r"\d+_\d+", pid):
            raise ProdigalFormatError(f"{label}: header has no ID=k_n attribute: {hdr[:120]!r}")
        partial = attrs.get("partial", "")
        if partial not in ("00", "01", "10", "11"):
            raise ProdigalFormatError(f"{label}: bad partial code {partial!r} in {hdr[:120]!r}")
        if pid in genes:
            raise ProdigalFormatError(f"{label}: duplicated Prodigal ID {pid}")
        genes[pid] = {"start": int(m.group(2)), "end": int(m.group(3)),
                      "strand": "+" if m.group(4) == "1" else "-", "partial": partial, "seq": seq}
    return genes


def parse_prodigal_gff(path):
    """Return (seqs {k: (name, length, transl_table)}, cds {k_n: (start, end, strand)}, lines)."""
    seqs, cds, lines = {}, {}, []
    last_k = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            lines.append(line)
            if line.startswith("# Sequence Data:"):
                m = SEQDATA_RE.match(line)
                if not m:
                    raise ProdigalFormatError(f"GFF: malformed Sequence Data line: {line[:160]!r}")
                k = int(m.group(1))
                name = m.group(3).split()[0] if m.group(3).split() else f"seq{k}"
                seqs[k] = [name, int(m.group(2)), ""]
                last_k = k
            elif line.startswith("# Model Data:"):
                m = TT_RE.search(line)
                if last_k is not None and m:
                    seqs[last_k][2] = m.group(1)
            elif line and not line.startswith("#"):
                f = line.split("\t")
                if len(f) != 9 or f[2] != "CDS":
                    raise ProdigalFormatError(f"GFF: unexpected feature line: {line[:160]!r}")
                pid = parse_attrs(f[8]).get("ID", "")
                if pid in cds:
                    raise ProdigalFormatError(f"GFF: duplicated ID {pid}")
                cds[pid] = (int(f[3]), int(f[4]), f[6])
    return {k: tuple(v) for k, v in seqs.items()}, cds, lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--faa", required=True, help="output proteins (input too, unless --in-faa)")
    ap.add_argument("--fna", required=True, help="output CDS (input too, unless --in-fna)")
    ap.add_argument("--gff", default="", help="output GFF (input too, unless --in-gff)")
    ap.add_argument("--in-faa", default="")
    ap.add_argument("--in-fna", default="")
    ap.add_argument("--in-gff", default="")
    ap.add_argument("--genes-tsv", default="", help="per-gene coordinate table (Contract v2)")
    args = ap.parse_args()
    try:
        _run(args)
    except (ProdigalFormatError, OSError) as exc:
        utils.log(f"ERROR: {exc}")
        sys.exit(1)


def _run(args) -> None:
    in_faa = args.in_faa or args.faa
    in_fna = args.in_fna or args.fna
    in_gff = args.in_gff or args.gff
    if not in_gff:
        raise ProdigalFormatError("a Prodigal GFF (--in-gff/--gff) is required for contig "
                                  "indices and lengths")
    g = args.genome
    if utils.sanitize(g) != g:
        raise ProdigalFormatError(f"genome id {g!r} contains characters outside [A-Za-z0-9_.-]")

    prot = parse_prodigal_fasta(in_faa, "faa")
    nucl = parse_prodigal_fasta(in_fna, "fna")
    seqs, cds, gff_lines = parse_prodigal_gff(in_gff)

    if not prot:
        raise ProdigalFormatError(
            f"Prodigal predicted 0 genes for genome {g} ({len(seqs)} sequence record(s), "
            f"{sum(s[1] for s in seqs.values())} bp); refusing to write an empty proteome "
            "(it would read as 'Prometheus absent')")
    if set(prot) != set(nucl) or set(prot) != set(cds):
        raise ProdigalFormatError(
            f"Prodigal outputs disagree: {len(prot)} proteins, {len(nucl)} CDS, {len(cds)} GFF "
            f"features (only in faa: {sorted(set(prot) - set(nucl))[:3]}, only in fna: "
            f"{sorted(set(nucl) - set(prot))[:3]}, gff-only: {sorted(set(cds) - set(prot))[:3]})")

    rows, new_faa, new_fna, new_id_of = [], [], [], {}
    for pid, p in prot.items():
        k, n = (int(x) for x in pid.split("_"))
        if k not in seqs:
            raise ProdigalFormatError(f"gene {pid}: no '# Sequence Data: seqnum={k}' line in the GFF")
        q = nucl[pid]
        if (q["start"], q["end"], q["strand"], q["partial"]) != (p["start"], p["end"], p["strand"], p["partial"]):
            raise ProdigalFormatError(f"gene {pid}: faa and fna headers disagree")
        if cds[pid] != (p["start"], p["end"], p["strand"]):
            raise ProdigalFormatError(f"gene {pid}: GFF coordinates {cds[pid]} != FASTA "
                                      f"{(p['start'], p['end'], p['strand'])}")
        nt_len = p["end"] - p["start"] + 1
        if len(q["seq"]) != nt_len or nt_len % 3 or len(p["seq"]) * 3 != nt_len:
            raise ProdigalFormatError(f"gene {pid}: CDS length {len(q['seq'])} / span {nt_len} / "
                                      f"protein length {len(p['seq'])} are inconsistent")
        name, clen, tt = seqs[k]
        if p["end"] > clen:
            raise ProdigalFormatError(f"gene {pid}: end {p['end']} beyond contig length {clen}")
        new = utils.make_gene_id(g, k, n)
        new_id_of[pid] = new
        new_faa.append((new, p["seq"]))
        new_fna.append((new, q["seq"]))
        rows.append([new, name, k, n, p["start"], p["end"], p["strand"], p["partial"], clen,
                     len(p["seq"].rstrip("*")), tt])
    if len(set(new_id_of.values())) != len(new_id_of):
        raise ProdigalFormatError("gene ids are not unique after renaming")

    utils.write_fasta(new_faa, args.faa)
    utils.write_fasta(new_fna, args.fna)
    if args.gff:
        with open(args.gff, "w") as out:
            for line in gff_lines:
                if line and not line.startswith("#"):
                    f = line.split("\t")
                    pid = parse_attrs(f[8]).get("ID", "")
                    f[8] = GFF_ID_RE.sub(lambda m: f"{m.group(1)}ID={new_id_of[pid]}", f[8], count=1)
                    line = "\t".join(f)
                out.write(line + "\n")
    if args.genes_tsv:
        os.makedirs(os.path.dirname(os.path.abspath(args.genes_tsv)), exist_ok=True)
        with open(args.genes_tsv, "w") as out:
            out.write("\t".join(GENES_COLS) + "\n")
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")
    n_partial = sum(1 for r in rows if r[7] != "00")
    tts = sorted({r[10] for r in rows})
    utils.log(f"{g}: {len(rows)} genes on {len({r[2] for r in rows})}/{len(seqs)} contig(s) "
              f"({n_partial} partial; translation table(s) {tts}) renamed to {g}_<contig>_<gene>")
    if any(t not in ("", "11") for t in tts):
        utils.log(f"WARNING: {g}: Prodigal used translation table(s) {tts}; archaea use table 11 "
                  "(meta mode ignores -g and may pick a table-4 model)")


if __name__ == "__main__":
    main()
