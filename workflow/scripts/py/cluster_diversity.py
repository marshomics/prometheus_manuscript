"""Cluster validated Prometheus proteins into subfamilies at graded identities
(design §2.1).

For each requested sequence-identity threshold we run MMseqs2 `easy-cluster` with
a coverage floor and record, per protein, which cluster it falls into. The 30%
level approximates subfamilies; higher levels resolve recent diversification. The
output column suffixes are the *actual* identity values passed in (e.g.
cluster_0.30 ... cluster_0.90), so the schema follows the configured thresholds.

Robustness: clustering needs >=2 sequences and the `mmseqs` binary on PATH. If
either is missing we emit a header-only table (with whatever identity columns the
caller requested) so the DAG stays valid; per-protein cluster ids are then simply
absent rather than the step crashing.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import shutil
import subprocess


def parse_identities(s: str):
    """'0.30,0.50,0.70,0.90' -> ['0.30','0.50',...] preserving the literal text.

    The literal token is used as the column suffix so the header matches exactly
    what was requested (and the SCHEMA's cluster_0.30 etc.).
    """
    out = []
    for tok in str(s).split(","):
        tok = tok.strip()
        if tok:
            out.append(tok)
    return out


def header_for(identities):
    cols = ["protein_id", "genome_id"] + [f"cluster_{i}" for i in identities]
    return "\t".join(cols)


def run_mmseqs(faa, ident, coverage, threads, workdir):
    """Run MMseqs2 easy-cluster at one identity; return {protein -> cluster_id}.

    cluster id = the representative protein id from {prefix}_cluster.tsv
    (two columns: representative<TAB>member). Returns None on any failure so the
    caller can leave that column blank without aborting other identities.
    """
    prefix = os.path.join(workdir, f"clu_{ident}")
    tmp = os.path.join(workdir, f"tmp_{ident}")
    os.makedirs(workdir, exist_ok=True)
    cmd = [
        "mmseqs", "easy-cluster", faa, prefix, tmp,
        "--min-seq-id", str(ident), "-c", str(coverage),
        "--threads", str(threads),
    ]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as exc:
        utils.log(f"WARNING: mmseqs failed to launch for id={ident}: {exc}")
        return None
    if res.returncode != 0:
        utils.log(f"WARNING: mmseqs returncode={res.returncode} for id={ident}; "
                  f"stderr tail: {res.stderr[-300:]}")
        return None

    cluster_tsv = f"{prefix}_cluster.tsv"
    if not os.path.exists(cluster_tsv):
        utils.log(f"WARNING: missing {cluster_tsv} for id={ident}")
        return None

    mapping = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            rep, member = parts[0].strip(), parts[1].strip()
            mapping[member] = rep
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--identities", required=True,
                    help="comma list, e.g. 0.30,0.50,0.70,0.90")
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    identities = parse_identities(args.identities)
    header = header_for(identities)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    proteins = [pid for pid, _ in utils.read_fasta(args.faa)]
    n_seq = len(proteins)

    if n_seq < 2:
        utils.write_empty(args.out, header)
        utils.log(f"Only {n_seq} sequence(s); wrote header-only subfamilies.")
        return
    if shutil.which("mmseqs") is None:
        utils.write_empty(args.out, header)
        utils.log("mmseqs not on PATH; wrote header-only subfamilies.")
        return

    # cluster assignments per identity
    per_ident = {}  # ident -> {protein -> cluster}
    for ident in identities:
        mapping = run_mmseqs(args.faa, ident, args.coverage, args.threads, args.workdir)
        per_ident[ident] = mapping or {}

    # If every identity failed, degrade gracefully to header-only.
    if all(len(m) == 0 for m in per_ident.values()):
        utils.write_empty(args.out, header)
        utils.log("All mmseqs runs failed/empty; wrote header-only subfamilies.")
        return

    # Build one row per protein.
    rows = []
    for pid in proteins:
        row = [pid, utils.genome_of(pid)]
        for ident in identities:
            row.append(per_ident[ident].get(pid, ""))
        rows.append(row)

    with open(args.out, "w") as out:
        out.write(header + "\n")
        for row in rows:
            out.write("\t".join(str(x) for x in row) + "\n")

    n_clusters = {i: len(set(v for v in per_ident[i].values() if v)) for i in identities}
    utils.log(f"Subfamilies: {n_seq} proteins; clusters per identity {n_clusters}.")


if __name__ == "__main__":
    main()
