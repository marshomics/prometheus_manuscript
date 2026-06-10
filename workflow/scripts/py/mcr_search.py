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

Acceptance (uses methanogen.kofam_score_fraction): a hit is accepted when its
full-sequence bit score is >= score_fraction * THRESHOLD, where THRESHOLD is the
curated adaptive cutoff for that profile -- taken from the HMM's GA (gathering) line
if present, else from a KOfam ko_list (--ko-list) keyed by KO. When no curated
threshold is available we fall back to a strict full-sequence E-value <= 1e-10. This
makes the documented score_fraction parameter operative (the previous code ignored it
and only ever applied --cut_ga or the E-value floor) while staying anchored to the
curated per-profile thresholds. All guard against spurious cross-hits to other
radical-SAM / nickel enzymes; McrA is conserved enough that true hits clear the cutoff
by a wide margin.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import glob
import re
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


def hmm_ga_threshold(hmm_path: str):
    """Return the sequence GA (gathering) bit-score threshold from an HMM, or None.

    HMMER GA line: 'GA   <seq_score> <dom_score>;'. We use the sequence threshold.
    """
    try:
        with open(hmm_path) as fh:
            for line in fh:
                if line.startswith("GA "):
                    toks = line[3:].replace(";", " ").split()
                    if toks:
                        try:
                            return float(toks[0])
                        except ValueError:
                            return None
                if line.startswith("HMM "):  # GA appears before the HMM matrix
                    break
    except Exception:
        return None
    return None


def parse_ko_list(path: str) -> dict:
    """KofamKOALA ko_list -> {KO: adaptive_threshold(float)} (skips '-' thresholds)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    try:
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            try:
                ki = header.index("knum")
                ti = header.index("threshold")
            except ValueError:
                ki, ti = 0, 1  # KofamKOALA default column order
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(ki, ti):
                    continue
                ko = f[ki].strip()
                thr = f[ti].strip()
                if ko and thr and thr != "-":
                    try:
                        out[ko] = float(thr)
                    except ValueError:
                        continue
    except Exception as exc:
        utils.log(f"WARNING: could not parse ko_list {path}: {exc}")
    return out


def ko_from_stem(stem: str):
    """Extract a KEGG KO accession (K#####) from an HMM filename stem, or None."""
    m = re.search(r"K\d{5}", stem)
    return m.group(0) if m else None


def run_hmmsearch(hmm, proteome, threads):
    """Run hmmsearch (no built-in cutoff); return [(target, score, evalue)] from --tblout.

    tblout columns (per HMMER man page): col1 target name, col5 full-seq E-value,
    col6 full-seq score. Acceptance is applied by the caller against the curated
    threshold * score_fraction (or the E-value fallback). Returns None on failure.
    """
    fd, tbl = tempfile.mkstemp(suffix=".tblout")
    os.close(fd)
    cmd = ["hmmsearch", "--cpu", str(threads), "--tblout", tbl]
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
    ap.add_argument("--score-fraction", type=float, default=0.5,
                    help="accept hits with score >= fraction * curated threshold (GA/ko_list)")
    ap.add_argument("--ko-list", default="",
                    help="optional KofamKOALA ko_list giving per-KO adaptive thresholds")
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

    ko_thresholds = parse_ko_list(args.ko_list)
    fraction = float(args.score_fraction)

    # best hit per (genome, gene): keep max score
    best = {}  # (genome, gene) -> (protein_id, score, evalue)
    for hmm in hmms:
        stem = os.path.splitext(os.path.basename(hmm))[0]
        gene = gene_from_stem(stem)
        # curated adaptive threshold: HMM GA line, else ko_list by KO accession.
        threshold = hmm_ga_threshold(hmm)
        if threshold is None:
            ko = ko_from_stem(stem)
            if ko:
                threshold = ko_thresholds.get(ko)
        hits = run_hmmsearch(hmm, args.proteome, args.threads)
        if not hits:
            continue
        for target, score, evalue in hits:
            # Curated threshold available -> score >= fraction * threshold;
            # otherwise fall back to the strict E-value floor.
            if threshold is not None:
                if score < fraction * threshold:
                    continue
            elif evalue > EVALUE_CUTOFF:
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
