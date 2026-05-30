"""Detect the methyl-coenzyme M reductase genes mcrA/mcrB/mcrG by HMM search
(design §3.3 — methanogen status is defined from genes, not taxonomy).

For every profile HMM in --hmmdir we run hmmsearch against the combined proteome
and keep the best hit per gene per genome. McrA (encoded by mcrA / KO K00399) is
the diagnostic minimum for methanogenesis; mcrB (K00401) and mcrG (K00402) are
the other structural subunits. Detecting the complex from `mcr` rather than from
the order name keeps the trait independent of the GTDB taxonomy later used in the
association test.

Gene naming from the HMM filename stem:
  K00399 -> mcrA, K00401 -> mcrB, K00402 -> mcrG; if the stem already contains
  mcrA/mcrB/mcrG (case-insensitive) that wins; otherwise the bare stem is used.

Acceptance: if the HMM carries GA (gathering) thresholds we trust them via
--cut_ga (the curated, release-appropriate cutoff); otherwise we fall back to a
strict full-sequence E-value <= 1e-10. Both guard against spurious cross-hits to
other radical-SAM / nickel enzymes.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import glob
import shutil
import subprocess
import tempfile

HEADER = "genome_id\tgene\tprotein_id\tscore\tevalue"
EVALUE_CUTOFF = 1e-10
# KEGG KO -> mcr subunit (authoritative, rest.kegg.jp): K00399 = mcrA (alpha),
# K00401 = mcrB (beta), K00402 = mcrG (gamma).
KO_TO_GENE = {"K00399": "mcrA", "K00401": "mcrB", "K00402": "mcrG"}


def gene_from_stem(stem: str) -> str:
    """Infer the mcr subunit from a HMM filename stem."""
    s = stem.strip()
    low = s.lower()
    for g in ("mcra", "mcrb", "mcrg"):
        if g in low:
            return "mcr" + g[-1].upper()
    # KO accession (allow K00399 anywhere in the stem)
    for ko, gene in KO_TO_GENE.items():
        if ko.lower() in low:
            return gene
    return s


def hmm_has_ga(hmm_path: str) -> bool:
    """True if the HMM file declares GA (gathering) thresholds (line 'GA ...')."""
    try:
        with open(hmm_path) as fh:
            for line in fh:
                if line.startswith("GA "):
                    return True
                if line.startswith("HMM "):  # GA appears before the HMM matrix
                    break
    except Exception:
        return False
    return False


def run_hmmsearch(hmm, proteome, threads, use_ga):
    """Run hmmsearch; return list of (target, score, evalue) parsed from --tblout.

    tblout columns (per HMMER man page): col1 target name, col5 full-seq E-value,
    col6 full-seq score. Returns None on launch/exec failure.
    """
    fd, tbl = tempfile.mkstemp(suffix=".tblout")
    os.close(fd)
    cmd = ["hmmsearch", "--cpu", str(threads), "--tblout", tbl]
    if use_ga:
        cmd.append("--cut_ga")
    cmd += [hmm, proteome]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
    except Exception as exc:
        utils.log(f"WARNING: hmmsearch failed to launch for {hmm}: {exc}")
        try:
            os.remove(tbl)
        except OSError:
            pass
        return None
    if res.returncode != 0:
        utils.log(f"WARNING: hmmsearch returncode={res.returncode} for {hmm}; "
                  f"stderr tail: {res.stderr[-300:]}")

    hits = []
    try:
        with open(tbl) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split()
                if len(f) < 6:
                    continue
                target = f[0]
                try:
                    evalue = float(f[4])
                    score = float(f[5])
                except ValueError:
                    continue
                hits.append((target, score, evalue))
    finally:
        try:
            os.remove(tbl)
        except OSError:
            pass
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", required=True)
    ap.add_argument("--hmmdir", required=True)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    hmms = sorted(glob.glob(os.path.join(args.hmmdir, "*.hmm")))
    if not hmms:
        utils.write_empty(args.out, HEADER)
        utils.log("No HMMs in --hmmdir; wrote header-only mcr_hits.")
        return
    if shutil.which("hmmsearch") is None:
        utils.write_empty(args.out, HEADER)
        utils.log("hmmsearch not on PATH; wrote header-only mcr_hits.")
        return
    if not os.path.exists(args.proteome) or os.path.getsize(args.proteome) == 0:
        utils.write_empty(args.out, HEADER)
        utils.log("Empty proteome; wrote header-only mcr_hits.")
        return

    # best hit per (genome, gene): keep max score
    best = {}  # (genome, gene) -> (protein_id, score, evalue)
    for hmm in hmms:
        stem = os.path.splitext(os.path.basename(hmm))[0]
        gene = gene_from_stem(stem)
        use_ga = hmm_has_ga(hmm)
        hits = run_hmmsearch(hmm, args.proteome, args.threads, use_ga)
        if not hits:
            continue
        for target, score, evalue in hits:
            # When trusting GA, HMMER already filtered; otherwise enforce E-value.
            if not use_ga and evalue > EVALUE_CUTOFF:
                continue
            genome = utils.genome_of(target)
            key = (genome, gene)
            cur = best.get(key)
            if cur is None or score > cur[1]:
                best[key] = (target, score, evalue)

    rows = []
    for (genome, gene), (pid, score, evalue) in best.items():
        rows.append((genome, gene, pid, score, evalue))
    # deterministic order: genome, then gene
    rows.sort(key=lambda r: (r[0], r[1]))

    if not rows:
        utils.write_empty(args.out, HEADER)
        utils.log("No mcr hits passed thresholds; wrote header-only mcr_hits.")
        return

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for genome, gene, pid, score, evalue in rows:
            out.write(f"{genome}\t{gene}\t{pid}\t{score}\t{evalue}\n")

    utils.log(f"mcr_hits: {len(rows)} (genome,gene) hits across "
              f"{len(set(r[0] for r in rows))} genomes.")


if __name__ == "__main__":
    main()
