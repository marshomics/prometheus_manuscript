"""Build a small, fully known project for the end-to-end integration test.

World (17 input genomes, 13 species representatives after Phase 1):
  * outgroup phylum p__Synth_O: O1-O3                      (no ProA, no mcr)
  * methanogen phylum p__Synth_M: M1-M5 (+ M1b, same species as M1, lower quality)
                                                          (mcrA/B/G; ProA in all five)
  * non-methanogen phylum p__Synth_N: N1-N5              (ProA in N1-N3)
  * Q1 fails QC, U1 is unclassified by GTDB-Tk, B1 is a bacterium (all dropped in Phase 1)
Genome ids are NCBI-style accessions, so the GTDB-Tk reference-id guard (user_ prefix)
is exercised in both phases (the stand-in refuses bare accessions, as GTDB-Tk 2.4.0 does).

Planted genes (real sequences, back-translated): six UniProt ProA homologues with all four
domains, a DUF3320-less copy (A0A0F0LN71 cut after residue 1960, i.e. after REase_MTES),
and a haloarchaeal ResIII/Helicase_C fusion. Gene layout follows tests/preprocess_qc/synth.py
(random sense codons without start codons, AGGAGG RBS, stop codons in all frames between
genes), so real Prodigal calls exactly the planted genes.
Synthetic McrA/McrB/McrG families give real hmmsearch hits against HMMs built here, with
KOfam-format ko_list thresholds set well above background.
"""
from __future__ import annotations

import os
import random
import shutil
import stat
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(HERE.parent / "preprocess_qc"))
import synth  # noqa: E402  (format-faithful GTDB-Tk/CheckM2/BBMap stand-ins)
sys.path.insert(0, str(HERE.parents[1] / "workflow" / "scripts" / "py"))
import methanogenesis_search as _msearch  # noqa: E402  (the rules' KO universe)

RULES_YAML = (HERE.parents[1] / "workflow" / "resources" / "methanogenesis" / "rules" /
              "methanogenesis_rules.yaml")

FIX = REPO / "tests" / "fixtures"
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_W = [8.3, 1.4, 5.4, 6.7, 3.9, 7.1, 2.3, 5.9, 5.8, 9.6, 2.4, 4.1, 4.7, 3.9, 5.4, 6.6, 5.4, 6.9,
        1.1, 2.9]
# Non-start sense codons per amino acid (Met only has ATG; GTG/TTG are avoided so the
# planted ATG is the only alternative-free start).
SYN = {}
for _c, _a in synth.CODE11.items():
    if _a != "*" and _c not in ("GTG", "TTG"):
        SYN.setdefault(_a, []).append(_c)

TAX = {
    "O1": "d__Archaea;p__Synth_O;c__Synth_Oc;o__Synth_Oo;f__Synth_Of;g__Synth_Og;s__Synth_Og alpha",
    "O2": "d__Archaea;p__Synth_O;c__Synth_Oc;o__Synth_Oo;f__Synth_Of;g__Synth_Og;s__Synth_Og beta",
    "O3": "d__Archaea;p__Synth_O;c__Synth_Oc;o__Synth_Oo;f__Synth_Of;g__Synth_Og;s__Synth_Og gamma",
    "M1": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo;f__Synth_Mf;g__Synth_Mg1;s__Synth_Mg1 one",
    "M1b": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo;f__Synth_Mf;g__Synth_Mg1;s__Synth_Mg1 one",
    "M2": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo;f__Synth_Mf;g__Synth_Mg1;s__Synth_Mg1 two",
    "M3": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo;f__Synth_Mf;g__Synth_Mg2;s__Synth_Mg2 three",
    "M4": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo2;f__Synth_Mf2;g__Synth_Mg3;s__Synth_Mg3 four",
    "M5": "d__Archaea;p__Synth_M;c__Synth_Mc;o__Synth_Mo2;f__Synth_Mf2;g__Synth_Mg3;s__Synth_Mg3 five",
    "N1": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No;f__Synth_Nf;g__Synth_Ng1;s__Synth_Ng1 one",
    "N2": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No;f__Synth_Nf;g__Synth_Ng1;s__Synth_Ng1 two",
    "N3": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No;f__Synth_Nf2;g__Synth_Ng2;s__Synth_Ng2 three",
    "N4": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No2;f__Synth_Nf3;g__Synth_Ng3;s__Synth_Ng3 four",
    "N5": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No2;f__Synth_Nf3;g__Synth_Ng3;s__Synth_Ng3 five",
    "Q1": "d__Archaea;p__Synth_N;c__Synth_Nc;o__Synth_No2;f__Synth_Nf3;g__Synth_Ng3;s__Synth_Ng3 six",
    "U1": "Unclassified",
    "B1": "d__Bacteria;p__Synth_B;c__Synth_Bc;o__Synth_Bo;f__Synth_Bf;g__Synth_Bg;s__Synth_Bg one",
}
LABELS = list(TAX)
# NCBI-style accessions (deliberately identical in form to GTDB reference ids).
ACC = {lab: f"GCA_9100000{10 + i:02d}.1" for i, lab in enumerate(LABELS)}
ACC["B1"] = "GCF_910000099.1"
QC = {lab: (99.0, 0.5) for lab in LABELS}
QC["M1b"] = (93.0, 1.2)          # same species as M1, lower quality score -> collapsed
QC["Q1"] = (70.0, 1.0)           # fails completeness >= 90
PROA = {"M1": "A0A011V4U2", "M2": "A0A060RCK4", "M3": "A0A073IP09", "M4": "A0A0A6X990",
        "M5": "A0A0A7PKR5", "M1b": "A0A011V4U2", "N1": "A0A0F0LN71", "N2": "A0A0F0LN71/DUF3320-",
        "N3": "A0A1H8VVA9"}
# Ground truth for DefenseFinder: A0A0A6X990 is a genuine ProA homologue (full score 1024.8
# vs Prometheus__ProA, GA 600) whose alignment splits into two domains (441.8 + 574.1), each
# below the per-domain GA2 of 600, so `hmmsearch --cut_ga` (MacSyFinder) reports no hit and
# DefenseFinder correctly does not call it: M4 must come out as a non-carrier.
DF_MISSED = {"A0A0A6X990"}
METHANOGENS = {"M1", "M1b", "M2", "M3", "M4", "M5"}
# Planted methanogenesis pathway capacity (one KO per component; the rules' relaxations make a
# single detected subunit enough for the multi-subunit components, so the world stays small).
#   M1, M1b CO2-reducing      methyl branch + Moco + Mtr + Hdr + hydrogenase
#   M2      methyl-reducing   methanol MTs + Hdr + hydrogenase, no methyl branch
#   M3      acetoclastic      CODH/ACS + Mtr + HdrDE + AckA/Pta
#   M4      dismutating       methyl branch + Moco + Mtr + HdrDE + methanol MTs
#   M5      (none)            Mcr only: every pathway absent
# Fwd needs >= 60% of its 5 mandatory subunit slots (rules relax.min_fraction) plus the
# FwdA/FwdC signature, and Mtr >= 50% of its 6: planting one subunit each would leave the
# component missing, which is the rules' intended behaviour, not a bug.
MG_METHYL_BRANCH = ["K00200", "K00201", "K00202", "K00203", "K00672", "K01499", "K00319",
                    "K00320"]
MG_MOCO = ["K03639", "K03637", "K03635", "K03638", "K03750"]
MG_MTR = ["K00577", "K00578", "K00579", "K00584"]
MG_HDR_ABC = ["K03388", "K03389"]
MG_HDR_DE = ["K08264", "K08265"]
MG_FRH = ["K00440", "K00441"]
MG_MTA = ["K14080", "K04480", "K14081"]
MG_CODH = ["K00193", "K00197"]
MG_ACKPTA = ["K00925"]
MG_PLANT = {
    "M1":  MG_METHYL_BRANCH + MG_MOCO + MG_MTR + MG_HDR_ABC + MG_FRH,
    "M1b": MG_METHYL_BRANCH + MG_MOCO + MG_MTR + MG_HDR_ABC + MG_FRH,
    "M2":  MG_MTA + MG_HDR_ABC + MG_FRH,
    "M3":  MG_CODH + MG_MTR + MG_HDR_DE + MG_ACKPTA,
    "M4":  MG_METHYL_BRANCH + MG_MOCO + MG_MTR + MG_HDR_DE + MG_MTA,
    "M5":  [],
}
# Expected pathway calls (methanogenesis_pathways.py) for the species representatives.
MG_TRUTH = {
    "M1": {"co2_reducing": "present", "acetoclastic": "absent", "methylotrophic": "absent",
           "methyl_dismutating": "absent", "methyl_reducing": "absent"},
    "M2": {"co2_reducing": "absent", "acetoclastic": "absent", "methylotrophic": "present",
           "methyl_dismutating": "absent", "methyl_reducing": "present"},
    "M3": {"co2_reducing": "absent", "acetoclastic": "present", "methylotrophic": "absent",
           "methyl_dismutating": "absent", "methyl_reducing": "absent"},
    "M4": {"co2_reducing": "absent", "acetoclastic": "absent", "methylotrophic": "present",
           "methyl_dismutating": "present", "methyl_reducing": "absent"},
    "M5": {"co2_reducing": "absent", "acetoclastic": "absent", "methylotrophic": "absent",
           "methyl_dismutating": "absent", "methyl_reducing": "absent"},
}
MG_SUBTYPE = {"M1": "CO2-reducing", "M2": "methyl-reducing", "M3": "acetoclastic",
              "M4": "methyl-dismutating", "M5": "none_detected"}
REPS = ["O1", "O2", "O3", "M1", "M2", "M3", "M4", "M5", "N1", "N2", "N3", "N4", "N5"]
# Prometheus references for the validation tree: canonical, full-length ProA homologues not used
# as candidates (the HMM consensus, the four-domain A0A0A6X990 and the haloarchaeal fusion
# A0A2A2F1Z1). The sub-GA Pkinase fusions A0AAV5N545 / A0ABR6IN60 are NOT used: as references
# they grouped with the paralogs, left the paralog root edge at ~77 UFBoot and the candidates'
# Prometheus clade at 91-95 (Unplaced); with this set it is 98-100 across seeds.
REFS_PROMETHEUS = ["Prometheus__ProA-consensus", "A0A0A6X990", "A0A2A2F1Z1"]
# Species-level niche columns ("" = unknown -> NA, never 0)
NICHE = {lab: {"group_Marine": "0", "niche_anoxic": "0", "niche_thermophile": ""} for lab in LABELS}
for lab in ("M1", "M2", "M3", "N1", "O1"):
    NICHE[lab]["group_Marine"] = "1"
for lab in ("M1", "M2", "M3", "M4", "M5", "N3"):
    NICHE[lab]["niche_anoxic"] = "1"
for lab in ("M4", "M5", "N4"):
    NICHE[lab]["niche_thermophile"] = "1"
for lab in ("O1", "O2", "N1"):
    NICHE[lab]["niche_thermophile"] = "0"
NICHE["M1b"]["niche_thermophile"] = "1"      # unioned into species M1 ('all' aggregation)
# True species tree (outgroup p__Synth_O); the simulated marker alignment follows it.
TRUE_TREE = ("((O1:0.12,(O2:0.10,O3:0.11):0.08):0.9,"
             "(((M1:0.08,M1b:0.01):0.02,(M2:0.09,M3:0.10):0.07):0.12,(M4:0.09,M5:0.08):0.14):0.25,"
             "(((N1:0.10,N2:0.09):0.08,N3:0.12):0.10,(N4:0.08,(N5:0.07,Q1:0.02):0.05):0.12):0.25,"
             "U1:0.4,B1:1.5);")


def read_fasta(path):
    recs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        recs.append((name, "".join(buf)))
    return recs


def acc_of(header):
    p = header.split("|")
    return p[1] if len(p) >= 3 else header


def back_translate(prot, rnd):
    prot = prot.rstrip("*").upper()
    body = prot[1:] if prot.startswith("M") else prot
    cod = []
    for a in body:
        a = a if (a in SYN or a == "M") else "A"
        cod.append("ATG" if a == "M" else rnd.choice(SYN[a]))
    return "ATG" + "".join(cod)


def random_protein(rnd, n):
    return "M" + "".join(rnd.choices(AA, weights=AA_W, k=n - 1))


def mutate(prot, rnd, frac):
    s = list(prot)
    for i in range(1, len(s)):
        if rnd.random() < frac:
            s[i] = rnd.choice([a for a in AA if a != s[i]])
    return "".join(s)


def build_genome(rnd, planted, n_background=62):
    """planted: protein sequences (exact). Returns two contig sequences."""
    genes = [("bg", None) for _ in range(n_background)]
    for p in planted:
        genes.insert(rnd.randint(3, len(genes) - 3), ("p", p))
    split = len(genes) * 2 // 3
    contigs = []
    for block in (genes[:split], genes[split:]):
        seq = synth.spacer(30)
        for gi, (kind, prot) in enumerate(block):
            if kind == "bg":
                ncod = rnd.randint(180, 420)
                body = "ATG" + "".join(rnd.choice(synth.CODONS) for _ in range(ncod - 1))
            else:
                body = back_translate(prot, rnd)
            unit = synth.spacer(40) + synth.RBS + "ATAACT" + body + "TAA"
            seq += synth.rc(unit) if gi % 3 == 2 else unit
        seq += synth.spacer(60)
        contigs.append(seq)
    return contigs


def parse_newick(s):
    """Tiny newick reader for TRUE_TREE -> nested (name, branch_length, children)."""
    pos = 0

    def node():
        nonlocal pos
        kids = []
        if s[pos] == "(":
            pos += 1
            while True:
                kids.append(node())
                if s[pos] == ",":
                    pos += 1
                    continue
                pos += 1   # ')'
                break
        start = pos
        while s[pos] not in ",();":
            pos += 1
        name, _, bl = s[start:pos].partition(":")
        return (name, float(bl or 0), kids)
    return node()


def simulate_msa(rnd, n_sites=500):
    root = "".join(rnd.choices(AA, weights=AA_W, k=n_sites))
    out = {}

    def walk(nd, seq):
        name, bl, kids = nd
        s = "".join(rnd.choice(AA) if rnd.random() < bl else c for c in seq)
        if kids:
            for k in kids:
                walk(k, s)
        else:
            out[name] = s
    walk(parse_newick(TRUE_TREE), root)
    return out


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def make_exe(path, body, interp):
    write(path, f"#!{interp}\n{body}")
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


GTDBTK_PHASE2 = r'''
import gzip, os, re, sys
REF = re.compile(r"^(RS_|GB_)?GC[AF]_[0-9]{9}\.[0-9]+$")
sub, a = sys.argv[1], sys.argv[2:]
opt, i = {}, 0
while i < len(a):
    if a[i].startswith("--") and i + 1 < len(a) and not a[i + 1].startswith("--"):
        opt[a[i]] = a[i + 1]; i += 2
    else:
        opt[a[i]] = True; i += 1
out = opt["--out_dir"]
if sub == "identify":
    ids = []
    for line in open(opt["--batchfile"]):
        path, gid = line.rstrip("\n").split("\t")[:2]
        if not os.path.isfile(path): sys.exit(f"stub gtdbtk identify: missing {path}")
        ids.append(gid)
    bad = [g for g in ids if REF.match(g)]
    if bad:   # GTDB-Tk 2.4.0 main._genomes_to_process
        sys.exit(f"You have {len(bad)} genomes with the same id as GTDB-Tk reference genomes: {bad[:3]}")
    os.makedirs(os.path.join(out, "identify"), exist_ok=True)
    open(os.path.join(out, "identify", "stub_ids.txt"), "w").write("\n".join(ids) + "\n")
    sys.exit(0)
if sub == "align":
    ids = open(os.path.join(opt["--identify_dir"], "identify", "stub_ids.txt")).read().split()
    msa, name = {}, None
    for line in open(os.environ["STUB_SPECIES_MSA"]):
        line = line.strip()
        if line.startswith(">"): name = line[1:]; msa[name] = ""
        elif line: msa[name] += line
    os.makedirs(os.path.join(out, "align"), exist_ok=True)
    with gzip.open(os.path.join(out, "align", "gtdbtk.ar53.user_msa.fasta.gz"), "wt") as fh:
        for g in ids:   # GTDB-Tk names MSA records by the batchfile id
            fh.write(f">{g}\n{msa[g[len('user_'):] if g.startswith('user_') else g]}\n")
    sys.exit(0)
sys.exit(f"stub gtdbtk: unsupported subcommand {sub}")
'''

GTDBTK_DISPATCH = r'''
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sub = sys.argv[1] if len(sys.argv) > 1 else ""
script = "gtdbtk_classify.py" if sub == "classify_wf" else "gtdbtk_phase2.py"
os.execv(sys.executable, [sys.executable, os.path.join(here, script)] + sys.argv[1:])
'''

GENOMAD = r'''
import os, sys, zlib
args = sys.argv[1:]
assert args[0] == "end-to-end" and {"--cleanup", "--restart", "--threads", "--splits"} <= set(args), args
pos = [a for a in args[1:] if not a.startswith("--") and not a.isdigit()]
inp, out, db = pos[0], pos[1], pos[2]
stem = os.path.basename(inp).rsplit(".", 1)[0]
if not os.path.isfile(os.path.join(db, "version.txt")): sys.exit("database missing")
seqs, name = {}, None
for line in open(inp):
    line = line.strip()
    if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = []
    elif line: seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
sd = os.path.join(out, f"{stem}_summary"); os.makedirs(sd, exist_ok=True)
hdr = ("seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\t"
       "n_hallmarks\tmarker_enrichment\ttaxonomy\n")
k = zlib.crc32(stem.encode()) % 3          # 0, 1 or 2 proviruses per genome (deterministic)
c1 = list(seqs)[0]
with open(os.path.join(sd, f"{stem}_virus_summary.tsv"), "w") as s, \
        open(os.path.join(sd, f"{stem}_virus.fna"), "w") as f:
    s.write(hdr)
    for j in range(k):
        a = 1000 + j * 9000; b = a + 6000 - 1
        prov = seqs[c1][a - 1:b]
        pn = f"{c1}|provirus_{a}_{b}"
        s.write(f"{pn}\t{len(prov)}\tProvirus\t{a}-{b}\t7\t11\t0.9500\tNA\t2\t3.0\tViruses\n")
        f.write(f">{pn}\n{prov}\n")
'''

CHECKV = r'''
import os, sys
a = sys.argv[1:]
assert a[0] == "end_to_end" and {"--restart", "--remove_tmp", "-d", "-t"} <= set(a), a
inp, out = a[1], a[2]
db = a[a.index("-d") + 1]
for f in ("genome_db", "hmm_db"):
    if not os.path.isdir(os.path.join(db, f)): sys.exit(f"Error: database dir not found '{db}'")
ids = [l[1:].split()[0] for l in open(inp) if l.startswith(">")]
if not ids:   # CheckV 1.0.3: bare sys.exit() on an empty FASTA (exit 0, no output)
    sys.stderr.write("You input FASTA file is empty or not properly formatted."); sys.exit()
os.makedirs(out, exist_ok=True)
cols = ["contig_id", "contig_length", "provirus", "proviral_length", "gene_count", "viral_genes",
        "host_genes", "checkv_quality", "miuvig_quality", "completeness", "completeness_method",
        "contamination", "kmer_freq", "warnings"]
with open(os.path.join(out, "quality_summary.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for c in ids:
        fh.write("\t".join([c, "6000", "No", "NA", "7", "3", "0", "Medium-quality", "Genome-fragment",
                            "61.0", "AAI-based", "0.0", "1.0", ""]) + "\n")
'''


def build(root: Path, python_exe: str, seed: int = 20261030):
    """Write the whole project under root; return a dict describing the planted truth."""
    rnd = random.Random(seed)
    root = Path(root)
    homs = {acc_of(n): s for n, s in read_fasta(FIX / "proA_homologs.faa")}
    proa = dict(homs)
    proa["A0A0F0LN71/DUF3320-"] = homs["A0A0F0LN71"][:1960]
    # ---- synthetic mcr families + HMMs + ko_list --------------------------------
    mcr = {"K00399": ("mcrA", 560), "K00401": ("mcrB", 430), "K00402": ("mcrG", 250)}
    # Every KO the rules reference needs a profile: methanogenesis_search.py fails on a missing
    # one (a KO that is not searched must never read as "absent"). Only MG_PLANT is planted.
    mgen_kos = sorted(set(_msearch.rules_kos(_msearch.load_yaml(RULES_YAML))))
    planted_set = {k for v in MG_PLANT.values() for k in v}
    assert planted_set <= set(mgen_kos), sorted(planted_set - set(mgen_kos))
    fams = dict(mcr)
    for i, ko in enumerate(mgen_kos):
        fams[ko] = (ko, 240 + 7 * i)              # distinct lengths keep the families apart
    hmmdir = root / "refs" / "mcr_hmm"
    hmmdir.mkdir(parents=True)
    mgdir = root / "db" / "kofam" / "methanogenesis_profiles"
    mgdir.mkdir(parents=True)
    seeds, planted_mcr = {}, {}
    for ko, (_gene, n) in fams.items():
        seeds[ko] = random_protein(rnd, n)
        sto = hmmdir / f"{ko}.sto"
        with open(sto, "w") as fh:
            fh.write("# STOCKHOLM 1.0\n")
            for i in range(24):
                fh.write(f"s{i} {mutate(seeds[ko], rnd, 0.12)}\n")
            fh.write("//\n")
        # mcr profiles go to BOTH directories: refs/mcr_hmm for mcr_search and the KOfam
        # subset for methanogenesis_search (fetch_kofam_subset.py stages them the same way,
        # because the rules file has an `mcr` component).
        dests = [mgdir] + ([hmmdir] if ko in mcr else [])
        subprocess.run(["hmmbuild", "--amino", "-n", ko, str(dests[0] / f"{ko}.hmm"), str(sto)],
                       check=True, capture_output=True)
        for extra in dests[1:]:
            shutil.copy(dests[0] / f"{ko}.hmm", extra / f"{ko}.hmm")
        sto.unlink()
    for lab in sorted(METHANOGENS):
        planted_mcr[lab] = {ko: mutate(seeds[ko], rnd, 0.12) for ko in mcr}
        for ko in MG_PLANT.get(lab, []):
            planted_mcr[lab][ko] = mutate(seeds[ko], rnd, 0.12)
    probe = root / "probe.faa"
    with open(probe, "w") as fh:
        for lab, d in planted_mcr.items():
            for ko, s in d.items():
                fh.write(f">{lab}_{ko}\n{s}\n")
        for ko, s in seeds.items():     # every family needs a score: a KO of the rules that is
            fh.write(f">SEED_{ko}\n{s}\n")   # planted in no genome still needs a threshold
    thr = {}
    for ko in fams:   # threshold = 40% of the lowest planted full-sequence score
        tbl = root / f"probe_{ko}.tbl"
        hmm = mgdir / f"{ko}.hmm"
        subprocess.run(["hmmsearch", "--tblout", str(tbl), "-o", os.devnull,
                        str(hmm), str(probe)], check=True, capture_output=True)
        sc = [float(l.split()[5]) for l in open(tbl)
              if not l.startswith("#") and l.split()[0].endswith(ko)]
        assert sc, f"no probe hit for family {ko}"
        thr[ko] = round(0.4 * min(sc), 2)
        tbl.unlink()
    probe.unlink()
    write(root / "db" / "kofam" / "ko_list",
          "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\talen\tmlen\t"
          "eff_nseq\tre/pos\tdefinition\n"
          + "".join(f"{ko}\t{thr[ko]}\tfull\tall\t1.0\t24\t24\t{fams[ko][1]}\t{fams[ko][1]}\t"
                    f"5.0\t0.6\tsynthetic family {fams[ko][0]}\n" for ko in fams))
    # ---- genomes + input metadata ---------------------------------------------------
    gdir = root / "data" / "genomes"
    gdir.mkdir(parents=True)
    truth = {"acc": dict(ACC), "proa": {}, "methanogens": sorted(ACC[l] for l in METHANOGENS),
             "reps": [ACC[l] for l in REPS], "mcr_thresholds": thr,
             "mgen_pathways": {ACC[l]: v for l, v in MG_TRUTH.items()},
             "mgen_subtype": {ACC[l]: v for l, v in MG_SUBTYPE.items()},
             "mgen_planted_kos": {ACC[l]: v for l, v in MG_PLANT.items()}}
    for lab in LABELS:
        planted = []
        if lab in PROA:
            planted.append(proa[PROA[lab]])
            truth["proa"][ACC[lab]] = PROA[lab]
            if PROA[lab] not in DF_MISSED:
                truth.setdefault("df_called", []).append(ACC[lab])
        if lab in planted_mcr:
            planted += list(planted_mcr[lab].values())
        contigs = build_genome(rnd, planted)
        synth.write_fasta(gdir / f"{ACC[lab]}.fasta", ["contig_1", "contig_2"], contigs)
    cols = ["group_Marine", "niche_anoxic", "niche_thermophile"]
    lines = ["genome_id\tfilename\tsource\t" + "\t".join(cols), "# synthetic integration set"]
    for lab in LABELS:
        lines.append(f"{ACC[lab]}\t{ACC[lab]}.fasta\tNCBI\t" + "\t".join(NICHE[lab][c] for c in cols))
    write(root / "data" / "input_metadata.tsv", "\n".join(lines) + "\n")
    write(root / "data" / "fake_tax.tsv", "".join(f"{ACC[l]}\t{TAX[l]}\n" for l in LABELS))
    write(root / "data" / "fake_checkm2.tsv",
          "".join(f"{ACC[l]} {QC[l][0]:.2f} {QC[l][1]:.2f}\n" for l in LABELS))
    msa = simulate_msa(rnd)
    write(root / "data" / "species_msa.faa", "".join(f">{ACC[l]}\n{msa[l]}\n" for l in LABELS))
    # ---- references and databases ----------------------------------------------------
    write(root / "refs" / "prometheus_reference.faa",
          "".join(f">{n}\n{homs[n]}\n" for n in REFS_PROMETHEUS))
    shutil.copy(FIX / "paralogs.faa", root / "refs" / "paralog_reference.faa")
    (root / "db" / "pfam").mkdir(parents=True)
    shutil.copy(FIX / "pfam_proA_subset.hmm", root / "db" / "pfam" / "Pfam-A.hmm")
    subprocess.run(["hmmpress", "-f", str(root / "db" / "pfam" / "Pfam-A.hmm")], check=True,
                   capture_output=True)
    (root / "db" / "gtdbtk").mkdir(parents=True)
    (root / "db" / "checkm2.dmnd").write_bytes(b"\0" * 64)
    write(root / "db" / "genomad_db" / "version.txt", "1.7\n")
    for sub in ("genome_db", "hmm_db"):
        (root / "db" / "checkv-db" / sub).mkdir(parents=True)
    # ---- stand-ins (only for tools whose reference data cannot be installed here) ----
    fb = root / "fakebin"
    synth.write_fake_tools(fb, python_exe)                 # gtdbtk (classify_wf), checkm2, statswrapper.sh
    os.replace(fb / "gtdbtk", fb / "gtdbtk_classify.py")
    make_exe(fb / "gtdbtk_phase2.py", GTDBTK_PHASE2, python_exe)
    make_exe(fb / "gtdbtk", GTDBTK_DISPATCH, python_exe)
    make_exe(fb / "genomad", GENOMAD, python_exe)
    make_exe(fb / "checkv", CHECKV, python_exe)
    synth.write_fake_gtdbtk_module(root / "fakepy")
    return truth
