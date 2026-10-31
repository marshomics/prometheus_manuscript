"""Detect the methanogenesis pathway KOs (design §3.3) with KofamScan-exact thresholds.

KO set
------
Every KO referenced by a component of the rules file (--rules,
workflow/resources/methanogenesis/rules/methanogenesis_rules.yaml) is searched. Each KO needs a
profile in --hmmdir: a file <KO>.hmm (KOfam layout) or any *.hmm whose model NAME / ACC is the KO.
A KO without a profile, or without a curated threshold (ko_list or GA line), is an ERROR: an
incomplete profile set must never silently turn into "pathway absent".

Genomes searched
----------------
Pathway capacity is only defined for methanogens. With --methanogen-calls (methanogen_calls.tsv)
only the proteins of genomes with methanogen == True are searched; --all-genomes searches the
whole proteome (used by the benchmark and tests).

Acceptance (kofam.py, identical to mcr_search.py)
-------------------------------------------------
score >= --score-fraction x KOfam threshold on the KO's score_type (full / best-domain), GA
semantics for GA-only profiles; contig-edge fragments (Prodigal partial != '00') are accepted at
>= --partial-score-fraction x threshold (call = partial_rescue). A protein counts for at most one
KO: the accepted KO with the highest score / threshold ratio (ties: lower KO id); the others are
recorded as superseded in the audit table.

Stop-codon read-through (call = readthrough_TAG / readthrough_TGA)
-----------------------------------------------------------------
Methylamine methyltransferases (mtmB, mtbB, mttB) carry pyrrolysine encoded by an in-frame UAG,
and several Methanococcales enzymes (e.g. FwdB, FdhA, HdrA, FruA) carry selenocysteine encoded by
UGA. Prodigal stops at both codons, so these genes arrive as N-terminal fragments that fail
full-length thresholds. For every homologous hit (full-sequence E <= --readthrough-evalue) of a KO
listed for codon C in the rules file's `readthrough` block whose gene is not accepted at the full
threshold, is complete at its 3' end (not a contig-edge fragment) and ends in C, the gene is
re-translated in frame from its start through C (C -> the rules' residue letter: O = pyrrolysine,
U = selenocysteine; Easel scores O as K and U as C) to the next stop codon, using the genome
sequence (--genome-dir/{genome}.{ext}). The read-through protein is re-scored against the same
profile and accepted only at the FULL threshold. Downstream Prodigal genes lying inside the
read-through span on the same strand are recorded as merged_into (they are part of the same
gene and are not counted again). Read-through is never applied to other KOs.

Near-threshold hits (call = near_threshold)
-------------------------------------------
A KO whose best hit reaches >= --near-threshold-fraction x threshold but not the threshold is
written with call = near_threshold. These are NEVER presence: the pathway caller uses them only
to refuse a confident ABSENCE, because a curated KOfam threshold optimises F-measure over all
taxa and can sit just above a divergent archaeal orthologue (e.g. KOfam K00625 phosphotrans-
acetylase: Methanosarcina acetivorans C2A scores 379.0 against the threshold 387.57, a 2%
shortfall, although the organism is a textbook acetoclastic methanogen).

Outputs
-------
--out      methanogenesis_ko_hits.tsv (PLAN contract): one row per accepted protein: genome_id,
           protein_id, ko, score, best_domain_score, criterion_score, threshold, score_type,
           threshold_source, call {full_threshold, partial_rescue, readthrough_TAG,
           readthrough_TGA}, partial, contig_index, gene_index, start, end, strand, then
           (additive) readthrough_length_aa, merged_proteins.
--out-all  audit: every hit with full-sequence E <= 1e-3 or accepted, with its decision.
--out-readthrough  FASTA of every attempted read-through translation (audit).
Failure: missing profile / threshold / genes table / genome FASTA needed for read-through, a
crashed hmmsearch or an inconsistent protein map is an error. An empty proteome or zero
methanogens gives header-only outputs (nonfatal: empty input).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils
import kofam
import kegg_module

import argparse
import csv
import gzip
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor

HIT_COLS = ["genome_id", "protein_id", "ko", "score", "best_domain_score", "criterion_score",
            "threshold", "score_type", "threshold_source", "call", "partial", "contig_index",
            "gene_index", "start", "end", "strand", "readthrough_length_aa", "merged_proteins"]
ALL_COLS = ["genome_id", "protein_id", "ko", "model", "score", "evalue", "best_domain_score",
            "best_domain_evalue", "criterion_score", "threshold", "scaled_threshold", "score_type",
            "threshold_source", "partial", "decision"]
STOPS = {"TAA", "TAG", "TGA"}
CODON = {}
_bases = "TCAG"
_aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for _i, _a in enumerate(_aas):
    CODON[_bases[_i // 16] + _bases[(_i // 4) % 4] + _bases[_i % 4]] = _a
COMP = str.maketrans("ACGTRYKMBDHVNacgtrykmbdhvn", "TGCAYRMKVHDBNtgcayrmkvhdbn")


def load_yaml(path):
    try:
        import yaml
    except ImportError:
        raise ImportError("PyYAML is required to read the methanogenesis rules file")
    with open(path) as fh:
        return yaml.safe_load(fh)


def rules_kos(rules: dict) -> list:
    """Every KO referenced by the rules' components."""
    out = set()
    comps = rules.get("components") or {}
    if not comps:
        raise ValueError("rules file defines no components")
    for name, c in comps.items():
        if "expr" not in c:
            raise ValueError(f"component {name} has no expr")
        out |= kegg_module.kos(kegg_module.parse(c["expr"]))
    return sorted(out)


def readthrough_table(rules: dict) -> dict:
    """{codon: {'residue': letter, 'kos': set}} from the rules' readthrough block."""
    out = {}
    for codon, spec in (rules.get("readthrough") or {}).items():
        codon = str(codon).upper()
        if codon not in ("TAG", "TGA"):
            raise ValueError(f"readthrough codon {codon} must be TAG or TGA")
        res = str(spec.get("residue", "")).upper()
        if res not in ("O", "U"):
            raise ValueError(f"readthrough residue for {codon} must be O or U")
        out[codon] = {"residue": res, "kos": set(spec.get("kos") or [])}
    return out


def map_profiles(hmmdir: str, want) -> dict:
    """{KO: (file, model)} for every wanted KO; an unmatched KO is an error."""
    found = {}
    for path in kofam.hmm_files(hmmdir):
        stem = os.path.splitext(os.path.basename(path))[0]
        sk = kofam.ko_from_stem(stem)
        if sk and sk not in want and stem == sk:
            continue                                   # fast skip in a full KOfam directory
        models = kofam.parse_hmm_models(path)
        if not models:
            raise ValueError(f"{path} contains no HMMER3 model")
        for m in models:
            ko = None
            for tok in (m["name"], (m["acc"] or "").split(".")[0]):
                if tok in want:
                    ko = tok
                    break
            if ko is None and sk in want and len(models) == 1:
                ko = sk
            if ko and ko not in found:
                m["file"] = path
                found[ko] = m
    missing = [k for k in want if k not in found]
    if missing:
        raise ValueError(f"no profile in {hmmdir} for {len(missing)} rules KO(s): "
                         f"{', '.join(missing)} (fetch them with fetch_kofam_subset.py)")
    return found


def read_fasta(path: str) -> dict:
    seqs, name, buf = {}, None, []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def translate_readthrough(contig: str, start: int, end: int, strand: str, codon: str,
                          residue: str, max_extra_codons: int = 1500):
    """Read-through translation of a Prodigal gene that ends in `codon`.

    start/end are 1-based inclusive with the stop codon included (Prodigal / genes.tsv). Returns
    (protein, nt_start, nt_end) for the extended ORF (stop codon included), or None when the gene
    does not end in `codon` or no later in-frame stop exists before the contig end."""
    if strand == "+":
        cds = contig[start - 1:end].upper()
        tail = contig[end:end + 3 * max_extra_codons].upper()
    else:
        cds = contig[start - 1:end].translate(COMP)[::-1].upper()
        lo = max(0, start - 1 - 3 * max_extra_codons)
        tail = contig[lo:start - 1].translate(COMP)[::-1].upper()
    if len(cds) % 3 or cds[-3:] != codon:
        return None
    aa = [CODON.get(cds[i:i + 3], "X") for i in range(0, len(cds) - 3, 3)]
    if not aa or "*" in aa:
        return None                                    # internal stop: not a single read-through
    aa[0] = "M"                                        # initiator, as Prodigal translates GTG/TTG
    aa.append(residue)
    for i in range(0, len(tail) - 2, 3):
        c = tail[i:i + 3]
        if c in STOPS:
            n_nt = len(cds) + i + 3
            if strand == "+":
                return "".join(aa), start, start + n_nt - 1
            return "".join(aa), end - n_nt + 1, end
        aa.append(CODON.get(c, "X"))
    return None


def write_header_only(args, note):
    utils.write_empty(args.out, "\t".join(HIT_COLS))
    if args.out_all:
        utils.write_empty(args.out_all, "\t".join(ALL_COLS))
    if args.out_readthrough:
        utils.write_empty(args.out_readthrough, "")
    utils.log(note)


def methanogen_genomes(path: str) -> set:
    out = set()
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None or not {"genome_id", "methanogen"} <= set(rd.fieldnames):
            raise ValueError(f"{path} lacks genome_id/methanogen columns")
        for r in rd:
            if str(r["methanogen"]).strip().lower() in ("true", "1"):
                out.add(r["genome_id"].strip())
    return out


def subset_proteome(proteome: str, pmap: dict, keep: set, dest: str) -> int:
    n, write = 0, False
    with open(proteome) as fh, open(dest, "w") as out:
        for line in fh:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                if pid not in pmap:
                    raise KeyError(f"proteome record {pid} is not in the protein map")
                write = pmap[pid] in keep
                n += write
            if write:
                out.write(line)
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", required=True)
    ap.add_argument("--genes-dir", required=True)
    ap.add_argument("--hmmdir", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--ko-list", default="", help="KOfam ko_list (per-KO thresholds)")
    sel = ap.add_mutually_exclusive_group(required=True)
    sel.add_argument("--methanogen-calls", default="", help="methanogen_calls.tsv")
    sel.add_argument("--all-genomes", action="store_true")
    ap.add_argument("--genome-dir", default="", help="genome FASTA directory (read-through)")
    ap.add_argument("--genome-ext", default="fna")
    ap.add_argument("--score-fraction", type=float, default=1.0)
    ap.add_argument("--partial-score-fraction", type=float, default=0.5)
    ap.add_argument("--near-threshold-fraction", type=float, default=0.9,
                    help="record (never accept) hits >= this fraction of the threshold as "
                         "call=near_threshold; 0 disables")
    ap.add_argument("--readthrough", default="True", help="True/False: stop-codon read-through")
    ap.add_argument("--readthrough-evalue", type=float, default=1e-5)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-all", default="")
    ap.add_argument("--out-readthrough", default="")
    args = ap.parse_args(argv)

    for p in (args.out, args.out_all, args.out_readthrough):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    frac, pfrac = float(args.score_fraction), args.partial_score_fraction
    if not (0 < frac <= 10):
        raise ValueError("--score-fraction must be > 0")
    if pfrac is not None and not (0 < pfrac <= frac):
        raise ValueError("--partial-score-fraction must be in (0, --score-fraction]")
    nfrac = float(args.near_threshold_fraction)
    if not (0 <= nfrac < 1):
        raise ValueError("--near-threshold-fraction must be in [0, 1)")
    flag = str(args.readthrough).strip().lower()
    if flag not in ("true", "false", "1", "0", "yes", "no"):
        raise ValueError(f"--readthrough must be True or False (got {args.readthrough!r})")
    do_rt = flag in ("true", "1", "yes")

    rules = load_yaml(args.rules)
    want = rules_kos(rules)
    rt = readthrough_table(rules) if do_rt else {}
    stray = sorted(set().union(*(v["kos"] for v in rt.values())) - set(want)) if rt else []
    if stray:
        raise ValueError(f"readthrough KOs not used by any component: {stray}")
    profiles = map_profiles(args.hmmdir, set(want))
    ko_table = None
    if args.ko_list:
        if not os.path.exists(args.ko_list):
            raise FileNotFoundError(f"--ko-list {args.ko_list} does not exist")
        ko_table = kofam.parse_ko_list(args.ko_list)
    for ko, m in profiles.items():
        m["thr"] = kofam.curated_threshold(ko, m["ga"], ko_table)
        if m["thr"] is None:
            raise ValueError(f"profile {m['file']} ({ko}) has no curated threshold: no GA line and "
                             f"no ko_list threshold ({args.ko_list or '--ko-list not given'})")

    if not kofam.fasta_has_sequences(args.proteome):
        write_header_only(args, "Empty proteome; wrote header-only methanogenesis_ko_hits.")
        return 0                                       # nonfatal: empty input
    if shutil.which("hmmsearch") is None:
        raise FileNotFoundError("hmmsearch not on PATH")
    pmap = kofam.read_pmap(args.pmap)
    tmpd = tempfile.mkdtemp(prefix="mgsearch_")
    try:
        target = args.proteome
        if not args.all_genomes:
            keep = methanogen_genomes(args.methanogen_calls)
            if not keep:
                write_header_only(args, "No methanogen genome; wrote header-only "
                                        "methanogenesis_ko_hits.")
                return 0                               # nonfatal: empty input (no methanogens)
            target = os.path.join(tmpd, "methanogens.faa")
            n = subset_proteome(args.proteome, pmap, keep, target)
            utils.log(f"searching {n} proteins of {len(keep)} methanogen genome(s)")
        genes = kofam.GeneTable(args.genes_dir)

        # ---- search: one hmmsearch per profile file, run concurrently ------------------
        by_file = {}
        for ko, m in profiles.items():
            by_file.setdefault(m["file"], {})[m["name"]] = ko
        files = sorted(by_file)
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            results = list(ex.map(lambda f: kofam.run_hmmsearch(f, target, 1), files))
        hits = []
        for f, rows in zip(files, results):
            names = by_file[f]
            for row in rows:
                ko = names.get(row["query"])
                if ko is None:
                    continue
                tid = row["target"]
                if tid not in pmap:
                    raise KeyError(f"hmmsearch target {tid} is not in the protein map")
                genome = pmap[tid]
                m = profiles[ko]
                ok, crit, ratio, st = kofam.evaluate(row, m["thr"], frac)
                g = genes.get(genome, tid)
                decision = "accepted_full" if ok else "below_threshold"
                if not ok and pfrac is not None and g["partial"] not in ("00", "NA"):
                    if kofam.evaluate(row, m["thr"], pfrac)[0]:
                        decision = "accepted_partial"
                if decision == "below_threshold" and nfrac and \
                        kofam.evaluate(row, m["thr"], nfrac * frac)[0]:
                    decision = "near_threshold"
                hits.append({"genome_id": genome, "protein_id": tid, "ko": ko, "model": m["name"],
                             "score": row["score"], "evalue": row["evalue"],
                             "best_domain_score": row["dom_score"],
                             "best_domain_evalue": row["dom_evalue"], "criterion_score": crit,
                             "threshold": m["thr"]["T"], "scaled_threshold": st,
                             "score_type": m["thr"]["score_type"],
                             "threshold_source": m["thr"]["source"], "partial": g["partial"],
                             "decision": decision, "ratio": ratio, "gene": g,
                             "rt_len": "", "merged": ""})

        # ---- stop-codon read-through ------------------------------------------------------
        rt_records, merged = [], {}
        if rt:
            accepted_full = {(h["protein_id"], h["ko"]) for h in hits if h["decision"] == "accepted_full"}
            cands = []
            for h in hits:
                if h["decision"] == "accepted_full" or h["evalue"] > args.readthrough_evalue:
                    continue
                g = h["gene"]
                if "start" not in g:
                    raise ValueError(f"genes table of {h['genome_id']} lacks coordinates "
                                     "(contig/start/end/strand) needed for read-through")
                three_prime = g["partial"][1] if g["strand"] == "+" else g["partial"][0]
                if three_prime != "0":                 # 3' end runs off the contig: no stop codon
                    continue
                for codon, spec in rt.items():
                    if h["ko"] in spec["kos"] and (h["protein_id"], h["ko"]) not in accepted_full:
                        cands.append((h, codon, spec["residue"]))
            if cands:
                if not args.genome_dir:
                    raise ValueError("read-through candidates exist but --genome-dir is not set")
                genome_seqs = {}
                rt_faa = os.path.join(tmpd, "readthrough.faa")
                index = {}
                with open(rt_faa, "w") as out:
                    for h, codon, res in cands:
                        gid = h["genome_id"]
                        if gid not in genome_seqs:
                            path = os.path.join(args.genome_dir, f"{gid}.{args.genome_ext}")
                            if not os.path.exists(path):
                                raise FileNotFoundError(f"genome FASTA for read-through missing: {path}")
                            genome_seqs[gid] = read_fasta(path)
                        g = h["gene"]
                        contig = genome_seqs[gid].get(g["contig"])
                        if contig is None:
                            raise KeyError(f"contig {g['contig']} of {h['protein_id']} not in genome FASTA")
                        tr = translate_readthrough(contig, g["start"], g["end"], g["strand"], codon, res)
                        if tr is None:
                            continue
                        rid = f"{h['protein_id']}|rt{codon}|{h['ko']}"
                        if rid in index:
                            continue
                        index[rid] = (h, codon, tr)
                        out.write(f">{rid}\n{tr[0]}\n")
                        rt_records.append((rid, tr[0]))
                if index:
                    by_ko = {}
                    for rid, (h, codon, tr) in index.items():
                        by_ko.setdefault(h["ko"], []).append(rid)
                    for ko, rids in sorted(by_ko.items()):
                        m = profiles[ko]
                        sub = os.path.join(tmpd, f"rt_{ko}.faa")
                        with open(sub, "w") as out:
                            for rid in rids:
                                out.write(f">{rid}\n{index[rid][2][0]}\n")
                        for row in kofam.run_hmmsearch(m["file"], sub, 1):
                            if row["query"] != m["name"] or row["target"] not in index:
                                continue
                            h, codon, tr = index[row["target"]]
                            ok, crit, ratio, st = kofam.evaluate(row, m["thr"], frac)
                            if not ok:
                                continue
                            g = h["gene"]
                            lo, hi = tr[1], tr[2]
                            mates = []
                            for pid, rec in genes.table(h["genome_id"]).items():
                                in_frame = ((rec.get("start", 0) - lo) % 3 == 0 if g["strand"] == "+"
                                            else (hi - rec.get("end", 0)) % 3 == 0)
                                if (rec.get("contig") == g["contig"] and rec.get("strand") == g["strand"]
                                        and pid != utils.local_of(h["protein_id"]) and in_frame
                                        and rec["start"] >= lo and rec["end"] <= hi):
                                    mates.append(f"{h['genome_id']}{utils.ID_SEP}{pid}")
                            h.update(decision=f"accepted_readthrough_{codon}", score=row["score"],
                                     evalue=row["evalue"], best_domain_score=row["dom_score"],
                                     best_domain_evalue=row["dom_evalue"], criterion_score=crit,
                                     scaled_threshold=st, ratio=ratio, rt_len=len(tr[0]),
                                     merged=",".join(sorted(mates)))
                            for mp in mates:
                                merged[mp] = h["protein_id"]
        # proteins merged into a read-through gene are not counted on their own
        for h in hits:
            if h["protein_id"] in merged and h["decision"] in ("accepted_full", "accepted_partial"):
                h["decision"] = f"merged_into_{merged[h['protein_id']]}"

        # ---- one KO per protein ---------------------------------------------------------
        best = {}
        for h in hits:
            if h["decision"].startswith("accepted"):
                key = (h["decision"] != "accepted_partial", h["ratio"], h["score"], -int(h["ko"][1:]))
                cur = best.get(h["protein_id"])
                if cur is None or key > cur[0]:
                    best[h["protein_id"]] = (key, h)
        for h in hits:
            if h["decision"].startswith("accepted"):
                win = best[h["protein_id"]][1]
                if win is not h:
                    h["decision"] = f"superseded_by_{win['ko']}"

        call_of = {"accepted_full": "full_threshold", "accepted_partial": "partial_rescue",
                   "accepted_readthrough_TAG": "readthrough_TAG",
                   "accepted_readthrough_TGA": "readthrough_TGA",
                   "near_threshold": "near_threshold"}
        rows = []
        for pid, (_, h) in best.items():
            g = h["gene"]
            rows.append([h["genome_id"], pid, h["ko"], h["score"], h["best_domain_score"],
                         h["criterion_score"], h["threshold"], h["score_type"],
                         h["threshold_source"], call_of[h["decision"]], h["partial"],
                         g.get("contig_index", ""), g.get("gene_index", ""), g.get("start", ""),
                         g.get("end", ""), g.get("strand", ""), h["rt_len"], h["merged"]])
        # near-threshold evidence: the best such hit per (genome, KO) that was NOT accepted
        accepted_pairs = {(r[0], r[2]) for r in rows}
        near = {}
        for h in hits:
            if h["decision"] != "near_threshold" or (h["genome_id"], h["ko"]) in accepted_pairs:
                continue
            key = (h["genome_id"], h["ko"])
            if key not in near or h["ratio"] > near[key]["ratio"]:
                near[key] = h
        for (gid, ko), h in near.items():
            if h["protein_id"] in best:
                continue                       # that protein already counts for another KO
            g = h["gene"]
            rows.append([gid, h["protein_id"], ko, h["score"], h["best_domain_score"],
                         h["criterion_score"], h["threshold"], h["score_type"],
                         h["threshold_source"], "near_threshold", h["partial"],
                         g.get("contig_index", ""), g.get("gene_index", ""), g.get("start", ""),
                         g.get("end", ""), g.get("strand", ""), "", ""])
        rows.sort(key=lambda r: (r[0], r[2], r[1]))
        with open(args.out, "w") as out:
            out.write("\t".join(HIT_COLS) + "\n")
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")
        if args.out_all:
            with open(args.out_all, "w") as out:
                out.write("\t".join(ALL_COLS) + "\n")
                for h in sorted(hits, key=lambda h: (h["genome_id"], h["protein_id"], h["ko"])):
                    if h["evalue"] <= kofam.AUDIT_EVALUE or h["decision"] != "below_threshold":
                        out.write("\t".join(str(h[c]) for c in ALL_COLS) + "\n")
        if args.out_readthrough:
            with open(args.out_readthrough, "w") as out:
                for rid, seq in rt_records:
                    out.write(f">{rid}\n{seq}\n")
        n_rt = sum(1 for r in rows if r[9].startswith("readthrough"))
        n_near = sum(1 for r in rows if r[9] == "near_threshold")
        utils.log(f"methanogenesis_ko_hits: {len(rows) - n_near} accepted proteins in "
                  f"{len({r[0] for r in rows})} genomes; {n_rt} via stop-codon read-through "
                  f"({len(rt_records)} attempted); {n_near} near-threshold (evidence only); "
                  f"{len(profiles)} KO profiles.")
        return 0
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
