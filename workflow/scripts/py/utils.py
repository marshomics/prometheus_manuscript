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


def local_of(protein_id: str) -> str:
    """Per-genome gene id from a combined id ('G__G_3_17' -> 'G_3_17'); identity otherwise."""
    return protein_id.split(ID_SEP, 1)[1] if ID_SEP in protein_id else protein_id


def make_gene_id(genome_id: str, contig_index: int, gene_index: int) -> str:
    """Per-genome gene id '{genome}_{contig#}_{gene#}' (both 1-based).

    This is MacSyFinder/DefenseFinder 'gembase'-compatible: the replicon is everything
    before the LAST '_' ('{genome}_{contig#}'), so `defense-finder run --db-type gembase`
    treats each contig as its own replicon and never assembles a system across a
    contig junction. Gene adjacency = same contig index, consecutive gene index.
    """
    return f"{sanitize(genome_id)}_{int(contig_index)}_{int(gene_index)}"


def split_gene_id(gene_id: str):
    """'{genome}_{contig#}_{gene#}' (or a combined 'genome__...' id) ->
    (genome_prefix, contig_index, gene_index); None when not in that form."""
    loc = local_of(gene_id)
    parts = loc.rsplit("_", 2)
    if len(parts) != 3:
        return None
    try:
        return parts[0], int(parts[1]), int(parts[2])
    except ValueError:
        return None


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
            try:
                score = float(f[13])
                env_from, env_to = int(f[19]), int(f[20])
                hmm_from, hmm_to = int(f[15]), int(f[16])
                qlen = int(f[5])
            except (ValueError, IndexError):
                score, env_from, env_to, hmm_from, hmm_to, qlen = (
                    float("nan"), ali_from, ali_to, 0, 0, 0)
            hits.setdefault(query, []).append({
                "domain": target, "accession": acc,
                "start": ali_from, "end": ali_to, "ievalue": ievalue,
                "score": score, "env_from": env_from, "env_to": env_to,
                "hmm_from": hmm_from, "hmm_to": hmm_to, "qlen": qlen,
                "description": " ".join(f[22:]),
            })
    if skipped:
        log(f"WARNING: parse_domtbl skipped {skipped} malformed line(s) in {path}")
    return hits


# --------------------------------------------------------------------------- #
# Prometheus domain specification -- the SINGLE source of truth for every script
# that classifies hmmscan hits (validate_prometheus, architecture_typology,
# domain_trees, run_selection, map_disorder_to_alignment, plots).
#
# The spec is config `prometheus.domains`, passed to scripts as JSON (`--domains`).
# Its key ORDER is the canonical N->C order of ProA. Verified empirically
# (hmmscan --cut_ga of UniProt ProA homologues and of the DefenseFinder
# Prometheus__ProA HMM consensus against Pfam):
#     DUF4011 (PF13195) -> SF1 helicase AAA_11/AAA_12 (PF13086/PF13087)
#     -> REase_MTES_1575 (PF18741; Pfam ID "MTES_1575") -> DUF3320 (PF11784)
# DUF3320 is ABSENT in ~54% of Prometheus instances (van den Berg et al. 2024,
# Cell Host Microbe 32:1427), so it is configured `required: false`.
#
# Matching is EXACT (not substring): a hit belongs to a category when its Pfam
# accession, version suffix stripped (PF13086.13 -> PF13086), is listed in the
# category's `accessions`, or else when its Pfam ID (hmmscan target name,
# case-insensitive) equals one of the category's `names`. Accession matches are
# tried across all categories before any name match, and each hit is assigned to
# at most one category, so no hit is ever counted twice.
# --------------------------------------------------------------------------- #
def pfam_acc_base(acc: str) -> str:
    """'PF13086.13' -> 'PF13086' ('' for missing / '-')."""
    a = str(acc or "").strip()
    if a in ("", "-", "nan", "None"):
        return ""
    return a.split(".")[0].upper()


def load_domain_spec(spec) -> "Dict[str, dict]":
    """Normalise the --domains JSON (str) or a dict into an insertion-ordered spec.

    Each category becomes {"names": [lowercased exact IDs], "accessions": [base
    accessions], "required": bool}. `required` defaults to True when absent, so an
    old-style spec without the key keeps its every-domain-required semantics.
    """
    import json
    if spec is None or spec == "":
        return {}
    if isinstance(spec, str):
        spec = json.loads(spec)
    out: Dict[str, dict] = {}
    for cat, s in spec.items():
        s = s or {}
        out[str(cat)] = {
            "names": [str(n).strip().lower() for n in s.get("names", []) if str(n).strip()],
            "accessions": [pfam_acc_base(a) for a in s.get("accessions", []) if pfam_acc_base(a)],
            "required": bool(s.get("required", True)),
        }
    return out


def hit_category(hit: dict, spec: Dict[str, dict]):
    """Return the ONE category a domtbl hit belongs to, or None (exact matching)."""
    acc = pfam_acc_base(hit.get("accession", ""))
    if acc:
        for cat, s in spec.items():
            if acc in s["accessions"]:
                return cat
    name = str(hit.get("domain", "")).strip().lower()
    if name:
        for cat, s in spec.items():
            if name in s["names"]:
                return cat
    return None


def canonical_order(spec: Dict[str, dict]) -> List[str]:
    return list(spec.keys())


def required_categories(spec: Dict[str, dict]) -> List[str]:
    return [c for c, s in spec.items() if s.get("required", True)]


def architecture_tokens(hits: List[dict], spec: Dict[str, dict]) -> List[str]:
    """N->C list of category tokens for one protein.

    Hits are ordered by alignment start; consecutive hits of the same category
    (e.g. the AAA_11 / AAA_11 / AAA_12 pieces of the split SF1 helicase) collapse
    to a single token. A category that recurs after another category (a genuine
    duplication or rearrangement) is kept, so reorderings stay visible.
    """
    toks: List[str] = []
    for h in sorted(hits, key=lambda x: (x.get("start", 0), x.get("end", 0))):
        cat = hit_category(h, spec)
        if cat is None:
            continue
        if not toks or toks[-1] != cat:
            toks.append(cat)
    return toks


def is_canonical_order(tokens: List[str], spec: Dict[str, dict]) -> bool:
    """True when every category occurs once and in canonical (spec) order."""
    order = canonical_order(spec)
    idx = [order.index(t) for t in tokens if t in order]
    return len(idx) == len(set(idx)) and idx == sorted(idx)


def domain_intervals(hits: List[dict], spec: Dict[str, dict], gap: int = 0) -> Dict[str, List[Tuple[int, int]]]:
    """Per-category merged (start, end) alignment-coordinate intervals (1-based, inclusive).

    Overlapping or abutting hits of one category (within `gap` residues) are merged.
    """
    raw: Dict[str, List[Tuple[int, int]]] = {}
    for h in hits:
        cat = hit_category(h, spec)
        if cat is None:
            continue
        s, e = int(h.get("start", 0)), int(h.get("end", 0))
        if e < s:
            s, e = e, s
        raw.setdefault(cat, []).append((s, e))
    merged: Dict[str, List[Tuple[int, int]]] = {}
    for cat, ivs in raw.items():
        ivs.sort()
        out = [list(ivs[0])]
        for s, e in ivs[1:]:
            if s <= out[-1][1] + 1 + gap:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[cat] = [(a, b) for a, b in out]
    return merged


def other_domains(hits: List[dict], spec: Dict[str, dict]) -> List[str]:
    """N->C list of Pfam IDs of hits that match NO configured category (fusions,
    e.g. the ResIII/Helicase_C N-terminal extension of haloarchaeal ProA)."""
    seen: List[str] = []
    for h in sorted(hits, key=lambda x: (x.get("start", 0), x.get("end", 0))):
        if hit_category(h, spec) is None:
            nm = str(h.get("domain", ""))
            if nm and nm not in seen:
                seen.append(nm)
    return seen
