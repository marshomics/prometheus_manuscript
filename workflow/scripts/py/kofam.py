"""KofamScan-exact profile-HMM acceptance, shared by mcr_search.py and methanogenesis_search.py.

One implementation of the curated-threshold rule keeps every KO-based call in the pipeline
on the same definition (takaram/kofam_scan lib/kofam_scan/result/parser.rb and hit.rb):

* KOfam ko_list: a hit is significant iff score >= threshold, where score is the FULL-sequence
  bit score (tblout column 6) when the KO's score_type is 'full' and the best-1-domain bit score
  (tblout column 9) when it is 'domain'. A threshold scale (KofamScan -T / --threshold-scale)
  multiplies the threshold, and the scaled value is rounded to 5 decimals as KofamScan does.
* GA lines (Pfam / TIGRFAM style HMMs): HMMER --cut_ga semantics, full-sequence score >= GA1 and
  best-domain score >= GA2 (each scaled).
* A model with neither (KO absent from, or '-' in, the ko_list and no GA line) has no curated
  threshold: callers must treat that as an error; there is no E-value fallback.

hmmsearch is run with `-T 0 --noali --tblout`, as KofamScan does. A crashed hmmsearch raises.
Gene coordinates and Prodigal partial flags come from genes/{genome}/{genome}.genes.tsv (SCHEMA
Contract v2 item 2).
"""
import csv
import glob
import os
import re
import subprocess
import tempfile

import utils

AUDIT_EVALUE = 1e-3


def ko_from_stem(stem: str):
    m = re.search(r"K\d{5}", stem or "")
    return m.group(0) if m else None


def parse_hmm_models(path: str) -> list:
    """HMMER3 header fields per model: name, acc, desc, leng, ga=(GA1, GA2) or None."""
    models, cur = [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("HMMER3"):
                cur = {"name": "", "acc": "", "desc": "", "leng": None, "ga": None}
            elif cur is None:
                continue
            elif line.startswith("NAME "):
                cur["name"] = line[5:].strip()
            elif line.startswith("ACC "):
                cur["acc"] = line[4:].strip()
            elif line.startswith("DESC "):
                cur["desc"] = line[5:].strip()
            elif line.startswith("LENG "):
                cur["leng"] = int(line[5:].split()[0])
            elif line.startswith("GA "):
                vals = line[3:].replace(";", " ").split()
                if len(vals) >= 2:
                    cur["ga"] = (float(vals[0]), float(vals[1]))
                elif len(vals) == 1:
                    cur["ga"] = (float(vals[0]), float(vals[0]))
            elif line.startswith("HMM "):
                models.append(cur)
                cur = None
    return models


def parse_ko_list(path: str) -> dict:
    """KOfam ko_list -> {KO: {'threshold': float|None, 'score_type': 'full'|'domain'|None}}.

    The first line is a header (knum, threshold, score_type, profile_type, ...); '-' means not
    available (KofamScan KO.parse). Columns are found by header name.
    """
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            ki, ti, si = header.index("knum"), header.index("threshold"), header.index("score_type")
        except ValueError:
            raise ValueError(f"ko_list {path}: header lacks knum/threshold/score_type")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ki, ti, si) or not f[ki].strip():
                continue
            thr = f[ti].strip()
            st = f[si].strip()
            out[f[ki].strip()] = {
                "threshold": None if thr in ("", "-") else float(thr),
                "score_type": None if st in ("", "-") else st,
            }
    return out


def parse_tblout(path: str) -> list:
    """hmmsearch --tblout rows (HMMER 3.4 user guide; KofamScan parser positions)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 18:
                raise ValueError(f"malformed tblout line in {path}: {line[:120]!r}")
            rows.append({"target": f[0], "query": f[2], "query_acc": f[3],
                         "evalue": float(f[4]), "score": float(f[5]),
                         "dom_evalue": float(f[7]), "dom_score": float(f[8])})
    return rows


def run_hmmsearch(hmm: str, proteome: str, threads) -> list:
    """Run hmmsearch (-T 0, as KofamScan); raise on failure; return parsed tblout."""
    fd, tbl = tempfile.mkstemp(suffix=".tblout")
    os.close(fd)
    try:
        cmd = ["hmmsearch", "--cpu", str(threads), "-T", "0", "--noali",
               "-o", os.devnull, "--tblout", tbl, hmm, proteome]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"hmmsearch failed (exit {res.returncode}) for {hmm}: "
                               f"{res.stderr.strip()[-500:]}")
        return parse_tblout(tbl)
    finally:
        try:
            os.remove(tbl)
        except OSError:
            pass


def read_pmap(path: str) -> dict:
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None or not {"protein_id", "genome_id"} <= set(rd.fieldnames):
            raise ValueError(f"protein map {path} lacks protein_id/genome_id columns")
        for r in rd:
            out[r["protein_id"].strip()] = r["genome_id"].strip()
    return out


GENE_COLS = ("protein_id", "contig", "contig_index", "gene_index", "start", "end", "strand",
             "partial", "contig_length")


class GeneTable:
    """Per-gene records from genes/{genome}/{genome}.genes.tsv (read lazily, cached per genome).

    get(genome, protein_id) accepts the combined id (genome__local) or the local id and returns a
    dict with partial (two-character Prodigal code, '00' = complete), contig, contig_index,
    gene_index, start, end (1-based, inclusive, stop codon included as Prodigal reports it),
    strand and contig_length. A missing table or protein is an error (Contract v2 item 2).
    """

    def __init__(self, genes_dir: str):
        self.genes_dir = genes_dir
        self.cache = {}

    def table(self, genome: str) -> dict:
        if genome not in self.cache:
            path = os.path.join(self.genes_dir, genome, f"{genome}.genes.tsv")
            if not os.path.exists(path):
                raise FileNotFoundError(f"genes table missing for {genome}: {path} "
                                        "(Contract v2 item 2)")
            recs = {}
            with open(path) as fh:
                rd = csv.DictReader(fh, delimiter="\t")
                if rd.fieldnames is None or not {"protein_id", "partial"} <= set(rd.fieldnames):
                    raise ValueError(f"{path} lacks protein_id/partial columns")
                full = set(GENE_COLS) <= set(rd.fieldnames)
                for r in rd:
                    flag = r["partial"].strip()
                    rec = {"partial": flag.zfill(2) if flag.isdigit() else flag}
                    if full:
                        rec.update(contig=r["contig"].strip(),
                                   contig_index=int(r["contig_index"]),
                                   gene_index=int(r["gene_index"]),
                                   start=int(r["start"]), end=int(r["end"]),
                                   strand=r["strand"].strip(),
                                   contig_length=int(r["contig_length"]))
                    recs[r["protein_id"].strip()] = rec
            self.cache[genome] = recs
        return self.cache[genome]

    def get(self, genome: str, protein_id: str) -> dict:
        recs = self.table(genome)
        for key in (protein_id, utils.local_of(protein_id)):
            if key in recs:
                return recs[key]
        raise KeyError(f"protein {protein_id} not found in genes table of {genome}")


class PartialFlags(GeneTable):
    """Backward-compatible view: get() returns only the Prodigal partial flag ('NA' without a
    genes directory)."""

    def get(self, genome: str, protein_id: str) -> str:
        if not self.genes_dir:
            return "NA"
        return GeneTable.get(self, genome, protein_id)["partial"]


def scaled(t: float, f: float) -> float:
    return round(t * f, 5)  # KofamScan rounds the scaled threshold to 5 decimals


def evaluate(row: dict, thr: dict, f: float):
    """(passes, criterion_score, ratio, scaled_threshold) for one tblout hit."""
    if thr["source"] == "ko_list":
        s = row["dom_score"] if thr["score_type"] == "domain" else row["score"]
        st = scaled(thr["T"], f)
        return s >= st, s, s / thr["T"], st
    ga1, ga2 = thr["T"], thr["T2"]
    ok = row["score"] >= scaled(ga1, f) and row["dom_score"] >= scaled(ga2, f)
    ratio = min(row["score"] / ga1 if ga1 else float("inf"),
                row["dom_score"] / ga2 if ga2 else float("inf"))
    return ok, row["score"], ratio, scaled(ga1, f)


def curated_threshold(ko, ga, ko_table):
    """Threshold record for one model, or None when no curated threshold exists.

    ko_list (when the KO has a threshold) wins over GA, as in mcr_search; an unknown score_type
    is an error.
    """
    if ko and ko_table is not None:
        rec = ko_table.get(ko)
        if rec and rec["threshold"] is not None:
            if rec["score_type"] not in ("full", "domain"):
                raise ValueError(f"ko_list score_type for {ko} is {rec['score_type']!r} "
                                 "(expected full/domain)")
            return {"source": "ko_list", "score_type": rec["score_type"], "T": rec["threshold"]}
    if ga:
        return {"source": "GA", "score_type": "ga", "T": ga[0], "T2": ga[1]}
    return None


def hmm_files(hmmdir: str) -> list:
    if not os.path.isdir(hmmdir):
        raise FileNotFoundError(f"HMM directory {hmmdir} does not exist")
    return sorted(set(glob.glob(os.path.join(hmmdir, "*.hmm")) +
                      glob.glob(os.path.join(hmmdir, "*.HMM"))))


def fasta_has_sequences(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    with open(path) as fh:
        return any(line.startswith(">") for line in fh)
