"""Cluster validated Prometheus proteins into subfamilies at graded identities
(design §2.1).

For each identity threshold MMseqs2 `easy-cluster` is run with (MMseqs2 15-6f452
Parameters.cpp / Alignment.cpp, verified):
  --min-seq-id X         identity over the alignment length (--seq-id-mode 0)
  -c C --cov-mode 0      both query AND target must be >= C covered (bidirectional;
                         keeps fusion/fission length variants from being lumped by a
                         shared core alone)
  --alignment-mode 3     exact sequence identity (auto mode 0 also selects it when
                         -c > 0 and --min-seq-id > 0; explicit for robustness)
  --cluster-mode 0       greedy set cover (MMseqs2 default)
  --cluster-reassign 1   corrects cascaded-clustering members that violate the criteria
  -s 7.5                 sensitive prefilter (the 30% level is remote homology)
Output columns are cluster_<identity formatted %.2f> (cluster_0.30 ... cluster_0.90, as
in SCHEMA) whatever the literal form of --identities ('0.3' and '0.30' are the same
column). The cluster id is the MMseqs2 representative member. Every input protein must
appear exactly once per level; a crashed or missing mmseqs fails the step. Zero input
sequences -> header-only table; a single sequence is its own cluster at every level.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import shutil
import subprocess


def parse_identities(s: str):
    """'0.3,0.50,0.7' -> ['0.30', '0.50', '0.70'] (validated, unique, order kept)."""
    out = []
    for tok in str(s).split(","):
        tok = tok.strip()
        if not tok:
            continue
        x = float(tok)
        if not 0.0 < x <= 1.0:
            raise ValueError(f"identity {tok} outside (0, 1]")
        lab = f"{x:.2f}"
        if lab in out:
            raise ValueError(f"identity {tok} duplicates {lab}")
        out.append(lab)
    if not out:
        raise ValueError("no identity thresholds given")
    return out


def header_for(identities):
    return "\t".join(["protein_id", "genome_id"] + [f"cluster_{i}" for i in identities])


def mmseqs_cmd(faa, prefix, tmp, ident, coverage, threads, cov_mode=0, sensitivity=7.5,
               cluster_mode=0):
    return ["mmseqs", "easy-cluster", faa, prefix, tmp,
            "--min-seq-id", str(ident), "-c", str(coverage), "--cov-mode", str(cov_mode),
            "--alignment-mode", "3", "--seq-id-mode", "0",
            "--cluster-mode", str(cluster_mode), "--cluster-reassign", "1",
            "-s", str(sensitivity), "--threads", str(threads)]


def run_mmseqs(faa, ident, coverage, threads, workdir, cov_mode=0, sensitivity=7.5,
               cluster_mode=0):
    """One easy-cluster run -> {member: representative}. Raises on any failure."""
    prefix = os.path.join(workdir, f"clu_{ident}")
    tmp = os.path.join(workdir, f"tmp_{ident}")
    shutil.rmtree(tmp, ignore_errors=True)       # never reuse stale intermediate DBs
    os.makedirs(workdir, exist_ok=True)
    cmd = mmseqs_cmd(faa, prefix, tmp, ident, coverage, threads, cov_mode, sensitivity,
                     cluster_mode)
    utils.log("running: " + " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"mmseqs easy-cluster failed (exit {res.returncode}) at "
                           f"identity {ident}: {(res.stderr or res.stdout)[-1500:]}")
    cluster_tsv = f"{prefix}_cluster.tsv"
    if not os.path.exists(cluster_tsv):
        raise RuntimeError(f"mmseqs wrote no {cluster_tsv}")
    mapping = {}
    with open(cluster_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            rep, member = parts[0], parts[1]
            if member in mapping:
                raise RuntimeError(f"{member} assigned to two clusters at identity {ident}")
            mapping[member] = rep
    return mapping


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--identities", required=True, help="comma list, e.g. 0.30,0.50,0.70,0.90")
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cov-mode", type=int, default=0)
    ap.add_argument("--sensitivity", type=float, default=7.5)
    ap.add_argument("--cluster-mode", type=int, default=0)
    ap.add_argument("--summary", default="", help="optional per-level summary TSV")
    args = ap.parse_args()

    try:
        identities = parse_identities(args.identities)
        cov = float(args.coverage)
        if not 0.0 <= cov <= 1.0:
            raise ValueError(f"coverage {args.coverage} outside [0, 1]")
    except ValueError as exc:
        raise SystemExit(f"[cluster_diversity] FATAL: {exc}")
    header = header_for(identities)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    proteins = utils.fasta_ids(args.faa)
    if len(set(proteins)) != len(proteins):
        raise SystemExit(f"[cluster_diversity] FATAL: duplicated ids in {args.faa}")
    per_ident = {}
    if len(proteins) == 0:
        utils.write_empty(args.out, header)
        utils.log("No sequences: header-only subfamilies table.")
    else:
        if len(proteins) == 1:
            per_ident = {i: {proteins[0]: proteins[0]} for i in identities}
        else:
            if shutil.which("mmseqs") is None:
                raise SystemExit("[cluster_diversity] FATAL: mmseqs not on PATH")
            try:
                for ident in identities:
                    per_ident[ident] = run_mmseqs(args.faa, ident, cov, args.threads,
                                                  args.workdir, args.cov_mode,
                                                  args.sensitivity, args.cluster_mode)
            except RuntimeError as exc:
                raise SystemExit(f"[cluster_diversity] FATAL: {exc}")
        for ident, m in per_ident.items():
            missing = [p for p in proteins if p not in m]
            extra = sorted(set(m) - set(proteins))
            if missing or extra:
                raise SystemExit(f"[cluster_diversity] FATAL: identity {ident}: "
                                 f"{len(missing)} input protein(s) unclustered (e.g. "
                                 f"{missing[:3]}), {len(extra)} unknown member(s)")
        with open(args.out, "w") as out:
            out.write(header + "\n")
            for pid in proteins:
                out.write("\t".join([pid, utils.genome_of(pid)]
                                    + [per_ident[i][pid] for i in identities]) + "\n")

    if args.summary:
        with open(args.summary, "w") as fh:
            fh.write("identity\tn_proteins\tn_clusters\tn_singletons\tlargest_cluster\n")
            for ident in identities:
                m = per_ident.get(ident, {})
                sizes = {}
                for rep in m.values():
                    sizes[rep] = sizes.get(rep, 0) + 1
                fh.write(f"{ident}\t{len(m)}\t{len(sizes)}\t"
                         f"{sum(1 for v in sizes.values() if v == 1)}\t"
                         f"{max(sizes.values()) if sizes else 0}\n")
    utils.log(f"Subfamilies: {len(proteins)} proteins; clusters per identity "
              f"{ {i: len(set(per_ident.get(i, {}).values())) for i in identities} }")


if __name__ == "__main__":
    main()
