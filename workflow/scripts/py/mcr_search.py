"""Detect the methyl-coenzyme M reductase genes by profile-HMM search (design §3.3:
methanogen status is defined from genes, not taxonomy).

Profiles
--------
Every HMM file in --hmmdir (*.hmm / *.HMM; one or more HMMER3 models per file) is
searched with `hmmsearch -T 0 --tblout` against the combined proteome (one run per
file; KofamScan runs hmmsearch with the same -T 0). Each model is mapped to a subunit
from its KEGG KO (model NAME / ACC or the file stem: K00399 mcrA, K00401 mcrB, K00402
mcrG, K03421 mcrC, K03422 mcrD), its TIGRFAM/NCBIfam accession (TIGR03256 mcrA,
TIGR03257 mcrB, TIGR03259 mcrG, TIGR03260 mcrD; checked against InterPro) or an
explicit mcrA..mcrD in the model name / file stem. --kos (config methanogen.mcr_kos)
lists the subunits that MUST be searched: a listed KO without a matching model is an
error, so an incomplete profile set can never silently zero the subunit counts. With
--kos only those subunits are searched (the directory may be the full KOfam profile
set); unmapped models are skipped and counted in the log.

Acceptance: the curated adaptive threshold, applied exactly as KofamScan does
--------------------------------------------------------------------------------
* KOfam ko_list (takaram/kofam_scan lib/kofam_scan/result/parser.rb and result/hit.rb):
  a hit is significant iff score >= threshold, where score is the FULL-sequence bit
  score (tblout column 6) when the KO's score_type is 'full' and the best-1-domain bit
  score (tblout column 9) when it is 'domain'. --score-fraction multiplies the
  threshold like KofamScan's -T/--threshold-scale (scaled threshold rounded to 5
  decimals, as KofamScan does); 1.0 (default) is KofamScan's own definition. ko_list
  values (genome.jp, 2026): K00399 McrA 741.47 full; K00401 McrB 506.80 full;
  K00402 McrG 359.37 full. A fraction of 0.5 would accept McrA-like hits at ~371 bits,
  where divergent McrA homologues (e.g. alkyl-coenzyme M reductases of alkane
  oxidisers) can pass.
* GA lines (TIGRFAM / Pfam style HMMs): HMMER --cut_ga semantics, i.e. full-sequence
  score >= GA1 and best-domain score >= GA2 (each scaled by --score-fraction).
* A model with no curated threshold (KO absent from, or '-' in, the ko_list and no GA
  line) is an error; there is no E-value fallback.

Truncated genes (--partial-score-fraction with --genes-dir)
-------------------------------------------------------------
A hit below the threshold is still accepted, as a FRAGMENT, when its gene runs off a
contig end (Prodigal partial flag != '00' in genes/{genome}/{genome}.genes.tsv,
Contract v2) and it scores >= partial_score_fraction x threshold: truncation, not
divergence, explains the lower score. A complete gene (partial '00') below the
threshold is never accepted.

Assignment
----------
A protein counts for at most one subunit (the model with the highest score/threshold
ratio). Per genome and subunit the reported protein is a full-threshold hit when one
exists (else the best fragment), highest score first, protein id as tie-break;
n_copies counts every accepted protein of that subunit.

Outputs
-------
--out      mcr_hits.tsv: SCHEMA columns genome_id, gene, protein_id, score (full-sequence
           bits), evalue, followed by additive audit columns (see HEADER).
--out-all  every hit with full-sequence E <= 1e-3 or accepted, with its decision.
Failure: a missing or empty HMM directory, a crashed hmmsearch, a ko_list that a profile
needs but that is missing, or a proteome / protein-map / genes.tsv inconsistency is an
error (exit != 0). Only an empty proteome yields header-only outputs.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import glob
import re
import shutil

# KofamScan-exact acceptance lives in kofam.py (shared with methanogenesis_search.py); the names
# are re-exported here so mcr_search keeps its public helpers.
from kofam import (AUDIT_EVALUE, PartialFlags, evaluate, fasta_has_sequences, ko_from_stem,  # noqa: F401
                   parse_hmm_models, parse_ko_list, parse_tblout, read_pmap, run_hmmsearch, scaled)

HEADER_COLS = ["genome_id", "gene", "protein_id", "score", "evalue",
               "best_domain_score", "criterion_score", "threshold", "score_type",
               "threshold_source", "model", "call", "partial", "n_copies"]
HEADER = "\t".join(HEADER_COLS)
ALL_COLS = ["genome_id", "protein_id", "model", "gene", "score", "evalue",
            "best_domain_score", "best_domain_evalue", "criterion_score", "threshold",
            "scaled_threshold", "score_type", "threshold_source", "partial", "decision"]

# KEGG KO -> mcr subunit (KEGG: K00399 alpha, K00401 beta, K00402 gamma, K03421 C, K03422 D)
KO_TO_GENE = {"K00399": "mcrA", "K00401": "mcrB", "K00402": "mcrG",
              "K03421": "mcrC", "K03422": "mcrD"}
# TIGRFAM / NCBIfam (InterPro: met_CoM_red_alp/bet/gam, met_CoM_red_D)
TIGR_TO_GENE = {"TIGR03256": "mcrA", "TIGR03257": "mcrB", "TIGR03259": "mcrG",
                "TIGR03260": "mcrD"}
GENES = ("mcrA", "mcrB", "mcrG", "mcrC", "mcrD")


def gene_from_stem(stem: str) -> str:
    """Backward-compatible helper: subunit from a name/stem, or the stem itself."""
    g, _ = map_model(stem, "", stem)
    return g or stem.strip()


def map_model(name: str, acc: str, stem: str):
    """(subunit or None, KO or None) for one HMM model."""
    toks = [t for t in (name, (acc or "").split(".")[0], stem) if t]
    for t in toks:
        for ko in re.findall(r"K\d{5}", t):
            if ko in KO_TO_GENE:
                return KO_TO_GENE[ko], ko
    for t in toks:
        for tg in re.findall(r"TIGR\d{5}", t):
            if tg in TIGR_TO_GENE:
                return TIGR_TO_GENE[tg], None
    for t in (name, stem):
        low = (t or "").lower()
        for g in GENES:
            if re.search(rf"(^|[^a-z]){g.lower()}([^a-z]|$)", low):
                return g, None
    return None, None


def hmm_ga_threshold(hmm_path: str):
    """Backward-compatible helper: sequence GA (GA1) of the first model, or None."""
    try:
        ms = parse_hmm_models(hmm_path)
    except OSError:
        return None
    return ms[0]["ga"][0] if ms and ms[0]["ga"] else None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", required=True)
    ap.add_argument("--hmmdir", required=True)
    ap.add_argument("--threads", default="1")
    ap.add_argument("--score-fraction", type=float, default=1.0,
                    help="threshold scale (KofamScan -T); 1.0 = KofamScan's definition")
    ap.add_argument("--partial-score-fraction", type=float, default=None,
                    help="accept contig-edge-truncated genes scoring >= this x threshold "
                         "(needs --genes-dir); omit to disable")
    ap.add_argument("--genes-dir", default="",
                    help="RESULTS/genes (genes/{g}/{g}.genes.tsv with Prodigal partial flags)")
    ap.add_argument("--ko-list", default="", help="KOfam ko_list (per-KO thresholds)")
    ap.add_argument("--kos", default="",
                    help="comma list of KOs whose subunits must be searched (methanogen.mcr_kos)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-all", default="", help="audit table of every homologous hit")
    args = ap.parse_args(argv)

    for p in (args.out, args.out_all):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    frac = float(args.score_fraction)
    pfrac = args.partial_score_fraction
    if not (0 < frac <= 10):
        raise ValueError("--score-fraction must be > 0")
    if pfrac is not None and not (0 < pfrac <= frac):
        raise ValueError("--partial-score-fraction must be in (0, --score-fraction]")
    if pfrac is not None and not args.genes_dir:
        raise ValueError("--partial-score-fraction needs --genes-dir (Prodigal partial flags)")

    # ---- profile set ------------------------------------------------------------
    if not os.path.isdir(args.hmmdir):
        raise FileNotFoundError(f"--hmmdir {args.hmmdir} does not exist")
    want_kos = [k.strip() for k in args.kos.split(",") if k.strip()]
    bad = [k for k in want_kos if k not in KO_TO_GENE]
    if bad:
        raise ValueError(f"--kos contains KOs that are not mcr subunits: {bad}")
    want_genes = {KO_TO_GENE[k] for k in want_kos}
    files = sorted(set(glob.glob(os.path.join(args.hmmdir, "*.hmm")) +
                       glob.glob(os.path.join(args.hmmdir, "*.HMM"))))
    if not files:
        raise FileNotFoundError(f"no *.hmm profiles in {args.hmmdir}")
    ko_table = None
    selected = []      # (file, [model dicts with gene/ko])
    n_unmapped = 0
    for path in files:
        stem = os.path.splitext(os.path.basename(path))[0]
        sk = ko_from_stem(stem)
        if want_genes and sk and KO_TO_GENE.get(sk) not in want_genes:
            continue                       # fast skip in a full KOfam profile directory
        models = parse_hmm_models(path)
        if not models:
            raise ValueError(f"{path} contains no HMMER3 model")
        keep = []
        for m in models:
            gene, ko = map_model(m["name"], m["acc"], stem)
            if gene is None or (want_genes and gene not in want_genes):
                n_unmapped += gene is None
                continue
            m.update(gene=gene, ko=ko, file=path)
            keep.append(m)
        if keep:
            selected.append((path, keep))
    if n_unmapped:
        utils.log(f"{n_unmapped} HMM model(s) map to no mcr subunit and were skipped.")
    found = {m["gene"] for _, ms in selected for m in ms}
    missing = sorted(want_genes - found)
    if missing:
        raise ValueError(f"no HMM in {args.hmmdir} for configured subunit(s) {missing} "
                         "(methanogen.mcr_kos)")
    if not selected:
        raise ValueError(f"no HMM in {args.hmmdir} maps to an mcr subunit")
    if "mcrA" not in found:
        raise ValueError("no McrA profile (K00399 / TIGR03256) found: McrA defines methanogens")

    # thresholds (curated only)
    for _, ms in selected:
        for m in ms:
            thr = None
            if m["ko"]:
                if ko_table is None and args.ko_list:
                    if not os.path.exists(args.ko_list):
                        raise FileNotFoundError(f"--ko-list {args.ko_list} does not exist")
                    ko_table = parse_ko_list(args.ko_list)
                rec = (ko_table or {}).get(m["ko"])
                if rec and rec["threshold"] is not None:
                    if rec["score_type"] not in ("full", "domain"):
                        raise ValueError(f"ko_list score_type for {m['ko']} is "
                                         f"{rec['score_type']!r} (expected full/domain)")
                    thr = {"source": "ko_list", "score_type": rec["score_type"],
                           "T": rec["threshold"]}
            if thr is None and m["ga"]:
                thr = {"source": "GA", "score_type": "ga", "T": m["ga"][0], "T2": m["ga"][1]}
            if thr is None:
                raise ValueError(
                    f"model {m['name'] or m['file']} ({m['gene']}) has no curated threshold: "
                    f"no GA line and KO {m['ko']} has no threshold in --ko-list "
                    f"{args.ko_list or '(not given)'}")
            m["thr"] = thr

    # ---- empty proteome: the only header-only case -------------------------------
    if not fasta_has_sequences(args.proteome):
        utils.write_empty(args.out, HEADER)  # nonfatal: empty proteome (no genomes / genes)
        if args.out_all:
            utils.write_empty(args.out_all, "\t".join(ALL_COLS))
        utils.log("Empty proteome; wrote header-only mcr_hits.")
        return 0
    if shutil.which("hmmsearch") is None:
        raise FileNotFoundError("hmmsearch not on PATH")
    pmap = read_pmap(args.pmap)
    flags = PartialFlags(args.genes_dir)

    # ---- search + per-hit decisions ----------------------------------------------
    hits = []
    for path, ms in selected:
        by_name = {m["name"]: m for m in ms}
        for row in run_hmmsearch(path, args.proteome, args.threads):
            m = by_name.get(row["query"])   # tblout query name == model NAME
            if m is None:
                continue                      # a model of this file that was not selected
            tid = row["target"]
            if tid not in pmap:
                raise KeyError(f"hmmsearch target {tid} is not in the protein map {args.pmap}")
            genome = pmap[tid]
            ok, crit, ratio, st = evaluate(row, m["thr"], frac)
            decision = "accepted_full" if ok else "below_threshold"
            partial = "NA"
            if not ok and pfrac is not None:
                ok_p, _, _, _ = evaluate(row, m["thr"], pfrac)
                if ok_p:
                    partial = flags.get(genome, tid)
                    if partial not in ("00", "NA"):
                        decision = "accepted_partial"
            hits.append({"genome_id": genome, "protein_id": tid, "model": m["name"] or
                         os.path.basename(path), "gene": m["gene"], "score": row["score"],
                         "evalue": row["evalue"], "best_domain_score": row["dom_score"],
                         "best_domain_evalue": row["dom_evalue"], "criterion_score": crit,
                         "threshold": m["thr"]["T"], "scaled_threshold": st,
                         "score_type": m["thr"]["score_type"],
                         "threshold_source": m["thr"]["source"], "partial": partial,
                         "decision": decision, "ratio": ratio})

    # partial flag for every accepted hit (audit), one subunit per protein
    best_for_protein = {}
    for h in hits:
        if h["decision"].startswith("accepted"):
            if h["partial"] == "NA" and args.genes_dir:
                h["partial"] = flags.get(h["genome_id"], h["protein_id"])
            cur = best_for_protein.get(h["protein_id"])
            key = (h["decision"] == "accepted_full", h["ratio"], h["score"])
            if cur is None or key > (cur["decision"] == "accepted_full", cur["ratio"], cur["score"]):
                best_for_protein[h["protein_id"]] = h
    for h in hits:
        win = best_for_protein.get(h["protein_id"])
        if h["decision"].startswith("accepted") and win is not h:
            h["decision"] = f"superseded_by_{win['gene']}"

    groups = {}
    for h in best_for_protein.values():
        groups.setdefault((h["genome_id"], h["gene"]), []).append(h)
    rows = []
    for (genome, gene), hs in sorted(groups.items()):
        hs.sort(key=lambda h: (h["decision"] != "accepted_full", -h["score"], h["protein_id"]))
        b = hs[0]
        rows.append([genome, gene, b["protein_id"], b["score"], b["evalue"],
                     b["best_domain_score"], b["criterion_score"], b["threshold"],
                     b["score_type"], b["threshold_source"], b["model"],
                     "full_threshold" if b["decision"] == "accepted_full" else "partial_rescue",
                     b["partial"], len(hs)])

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    if args.out_all:
        with open(args.out_all, "w") as out:
            out.write("\t".join(ALL_COLS) + "\n")
            for h in sorted(hits, key=lambda h: (h["genome_id"], h["protein_id"], h["model"])):
                if h["evalue"] <= AUDIT_EVALUE or h["decision"] != "below_threshold":
                    out.write("\t".join(str(h[c]) for c in ALL_COLS) + "\n")
    n_frag = sum(1 for r in rows if r[11] == "partial_rescue")
    utils.log(f"mcr_hits: {len(rows)} (genome, subunit) calls across "
              f"{len({r[0] for r in rows})} genomes ({n_frag} contig-edge fragments); "
              f"score_fraction={frac}, partial_score_fraction={pfrac}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
