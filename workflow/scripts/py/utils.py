"""Small shared helpers: FASTA/TSV IO, newick-safe IDs, empty-input guards."""
from __future__ import annotations
import os
import re
import sys
from typing import Iterator, Tuple, Dict, List

ID_SEP = "__"  # genome__origProteinId  (newick-safe; no pipes/colons)


def log(msg: str) -> None:
    print(f"[{os.path.basename(sys.argv[0])}] {msg}", file=sys.stderr)


def sanitize(label: str) -> str:
    """Make a label safe for newick / filenames (keep alnum, _, ., -)."""
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", str(label))


def read_fasta(path: str) -> Iterator[Tuple[str, str]]:
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            elif line and not line.startswith("#"):
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def write_fasta(records: List[Tuple[str, str]], path: str, width: int = 60) -> int:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    n = 0
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")
            n += 1
    return n


def fasta_ids(path: str) -> List[str]:
    return [name for name, _ in read_fasta(path)]


def genome_of(protein_id: str) -> str:
    """Recover the genome id from a combined protein id (genome__orig)."""
    return protein_id.split(ID_SEP, 1)[0]


def make_protein_id(genome_id: str, orig_id: str) -> str:
    return f"{sanitize(genome_id)}{ID_SEP}{sanitize(orig_id)}"


def write_empty(path: str, header: str) -> None:
    """Write a header-only TSV so downstream steps and the DAG stay valid."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as out:
        out.write(header.rstrip("\n") + "\n")


def parse_domtbl(path: str) -> Dict[str, List[dict]]:
    """Parse hmmscan --domtblout into {query_protein: [hit dicts]} (by ali coords).

    hmmscan --domtblout is whitespace-separated with 22 FIXED columns followed by
    an OPTIONAL free-text description (columns 23+). A domain hit whose target has
    no description has EXACTLY 22 fields, so we require len(f) >= 22 (NOT >= 23):
    requiring 23 silently dropped every description-less domain hit. 0-indexed
    columns used here: target/domain name=f[0], target accession=f[1],
    query(protein) name=f[3], i-Evalue=f[12], ali_from=f[17], ali_to=f[18].
    Multiple spaces between columns are handled by str.split(); index access is
    guarded so a short/garbled line is skipped (with a warning) rather than
    crashing the parse or zeroing out a real protein's other valid hits.
    """
    hits: Dict[str, List[dict]] = {}
    if not path or not os.path.exists(path):
        return hits
    skipped = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            # Need the 22 fixed columns; the description (col 23+) is optional and
            # ignored. A description-less hit has exactly 22 fields and MUST pass.
            if len(f) < 22:
                skipped += 1
                continue
            # hmmscan domtblout: col1 target(domain) name, col2 acc, col4 query name,
            # ali coords at 18,19 (env 20,21). Use ali_from/ali_to (cols 18,19, 1-based).
            target, acc, query = f[0], f[1], f[3]
            try:
                ali_from, ali_to = int(f[17]), int(f[18])
                ievalue = float(f[12])
            except (ValueError, IndexError):
                skipped += 1
                continue
            hits.setdefault(query, []).append({
                "domain": target, "accession": acc,
                "start": ali_from, "end": ali_to, "ievalue": ievalue,
            })
    if skipped:
        log(f"WARNING: parse_domtbl skipped {skipped} malformed line(s) in {path}")
    return hits
