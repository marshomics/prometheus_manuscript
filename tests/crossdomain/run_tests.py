#!/usr/bin/env python3
"""Synthetic test suite for the crossdomain stage (leaf 1.1.8).

    bash tests/run_leaf.sh crossdomain        # prints CROSSDOMAIN TESTS PASSED

Three layers:
  1. unit tests of every owned script on hand-built inputs with known answers
     (empty-input and failure behaviour included);
  2. an end-to-end Snakemake run of workflow/rules/crossdomain.smk on fixture-derived
     data (real hmmscan, MAFFT, ClipKIT, IQ-TREE, PAL2NAL, Rscript + the reused R/Python
     scripts; the archaeal producer rules are stubs whose outputs the test writes;
     geNomad/CheckV are not installed here and are replaced by fake executables that
     write their documented output layout), plus an all-empty run (zero bacterial
     Prometheus) that exercises every `# nonfatal:` branch;
  3. failure injection: a crashing hmmscan / IQ-TREE or a missing geNomad DB must make
     the rule fail (no placeholder output).
Scratch space is an absolute temp dir outside the repo (CROSSDOMAIN_TMP or /tmp); set
CROSSDOMAIN_KEEP=1 to keep it. At most 2 threads per tool. No bytecode is written.
"""
from __future__ import annotations

import os
import sys

sys.dont_write_bytecode = True
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import gzip  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import py_compile  # noqa: E402
import random  # noqa: E402
import re  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, run, FIXTURES, have, read_tsv, REPO, PYTHON, RSCRIPT,  # noqa: E402
                  SCRIPTS_PY, TOOLS_BIN)


def _scratch_root() -> Path:
    for cand in (os.environ.get("CROSSDOMAIN_TMP"), "/tmp", tempfile.gettempdir()):
        if not cand:
            continue
        try:
            base = Path(cand).resolve()
            base.mkdir(parents=True, exist_ok=True)
            return Path(tempfile.mkdtemp(prefix="cdtest_", dir=str(base)))
        except OSError:
            continue
    raise RuntimeError("no writable scratch directory")


ROOT = _scratch_root()
for _d in ("tmp", "home", "mpl", "pyc"):
    (ROOT / _d).mkdir()
os.environ["TMPDIR"] = str(ROOT / "tmp")           # MAFFT's mktemp needs an absolute TMPDIR
tempfile.tempdir = str(ROOT / "tmp")
os.environ["HOME"] = str(ROOT / "home")            # snakemake / matplotlib caches
os.environ["XDG_CACHE_HOME"] = str(ROOT / "home" / ".cache")
os.environ["MPLCONFIGDIR"] = str(ROOT / "mpl")
sys.path.insert(0, str(SCRIPTS_PY))

S = Suite("crossdomain")
OWNED_PY = ["extract_bacterial_prometheus.py", "build_bacterial_table.py",
            "build_bacterial_niche_master.py", "crossdomain_analysis.py",
            "plot_crossdomain_tree.py", "plot_crossdomain_extra.py",
            "plot_viral_crossdomain.py", "prune_tree_to_taxa.py"]
SPEC = {"DUF4011": {"names": ["DUF4011"], "accessions": ["PF13195"], "required": True},
        "DNA2": {"names": ["AAA_11", "AAA_12", "AAA_19", "AAA_30"],
                 "accessions": ["PF13086", "PF13087", "PF13245", "PF13604"], "required": True},
        "REase": {"names": ["MTES_1575"], "accessions": ["PF18741"], "required": True},
        "DUF3320": {"names": ["DUF3320"], "accessions": ["PF11784"], "required": False}}


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #
def mk(path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write(path, text):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


def tsv(path, header, rows):
    return write(path, "\t".join(header) + "\n" + "".join(
        "\t".join("" if v is None else str(v) for v in r) + "\n" for r in rows))


def fasta(path, recs):
    return write(path, "".join(f">{n}\n{s}\n" for n, s in recs))


def read_fa(path):
    out, name, buf = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def summary(path):
    return {r["metric"]: r["value"] for _, r in read_tsv(path).iterrows()}


def tips_of(newick_text):
    return sorted(re.findall(r"[(,]([A-Za-z0-9_.]+):", newick_text))


def fails(script, *args):
    p = py(script, *args, check=False)
    assert p.returncode != 0, f"{script} should have failed:\n{p.stdout[-800:]}\n{p.stderr[-800:]}"
    return p


BASES = "TCAG"
AA_TABLE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_AA = {a + b + c: AA_TABLE[16 * i + 4 * j + k]
            for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
SYN = {}
for _c, _a in CODON_AA.items():
    SYN.setdefault(_a, []).append(_c)


def back_translate(prot, rng):
    return "".join(rng.choice(sorted(SYN[a])) for a in prot) + "TAA"


def rand_codons(n, rng):
    sense = sorted(c for c, a in CODON_AA.items() if a != "*")
    return "".join(rng.choice(sense) for _ in range(n))


def find_snakemake():
    cands = [os.environ.get("PROM_SNAKEMAKE"), shutil.which("snakemake"),
             str(Path(PYTHON).resolve().parents[2] / "prom-smk" / "bin" / "snakemake")]
    for c in cands:
        if c and Path(c).exists():
            return c
    raise AssertionError("snakemake not found (set PROM_SNAKEMAKE); the end-to-end rule "
                         "tests cannot run")


# --------------------------------------------------------------------------- #
# 0. static checks
# --------------------------------------------------------------------------- #
@S.test
def py_compile_owned():
    for f in OWNED_PY:
        py_compile.compile(str(SCRIPTS_PY / f), cfile=str(ROOT / "pyc" / (f + "c")), doraise=True)


@S.test
def rules_hygiene():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"),
             str(REPO / "workflow" / "rules" / "crossdomain.smk")])
    assert "RULES HYGIENE OK" in p.stdout


@S.test
def rules_no_swallowed_failures_or_self_copy():
    txt = (REPO / "workflow" / "rules" / "crossdomain.smk").read_text()
    code = "\n".join(ln.split("#", 1)[0] for ln in txt.splitlines())
    assert "|| true" not in code and "|| cp" not in code and "|| echo" not in code
    assert "2>/dev/null" not in code
    # regression: the old rule cp'ed IQ-TREE's own output onto itself, failed, then
    # replaced the real tree with '();'
    assert "cp {params.workdir}/crossdomain.treefile" not in txt
    # regression: the codon alignment must get the COMBINED archaeal + bacterial list
    blk = txt.split("rule crossdomain_codon_align:")[1].split("\nrule ")[0]
    assert "--proteins {output.proteins}" in blk and "cat {input.arch} {input.bact}" in blk
    # regression: python is never run inside the R env
    blk = txt.split("rule bact_niche:")[1].split("\nrule ")[0]
    assert "python " not in blk and "rstats.yaml" in blk


def _rule_block(smk, rule):
    txt = (REPO / "workflow" / "rules" / smk).read_text()
    m = re.search(rf"^rule {rule}:\n(.*?)(?=^rule |\Z)", txt, re.S | re.M)
    assert m, f"rule {rule} not found in {smk}"
    return m.group(1)


def _call_flags(block, tool):
    """Option tokens of every `tool` invocation in a rule (backslash continuations)."""
    flags, lines, i = set(), block.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if re.search(rf"(^|[\s/]){re.escape(tool)}(\s|$)", ln) and not ln.strip().startswith("#"):
            cmd = [ln]
            while cmd[-1].rstrip().endswith("\\") and i + 1 < len(lines):
                i += 1
                cmd.append(lines[i])
            text = " ".join(c.split("#", 1)[0] for c in cmd)
            text = text[re.search(rf"(^|[\s/]){re.escape(tool)}", text).end():]
            flags |= set(re.findall(r"(?<=\s)(--?[A-Za-z][A-Za-z0-9_-]*)", " " + text))
        i += 1
    return flags


# (archaeal rule file, archaeal rule, bacterial rule, tool, flags the bacterial call may
# omit -- with the reason). The bacterial branch REUSES the archaeal scripts; the same
# evidence and thresholds on both sides is what makes the cross-domain comparison valid.
PARITY = [
    ("prometheus.smk", "prometheus_domain_scan", "bact_domain_scan", "hmmscan", {}),
    ("prometheus.smk", "prometheus_validation_tree", "bact_validation_tree", "combine_proteomes.py", {}),
    ("prometheus.smk", "prometheus_validation_tree", "bact_validation_tree", "mafft", {}),
    ("prometheus.smk", "prometheus_validation_tree", "bact_validation_tree", "clipkit", {}),
    ("prometheus.smk", "prometheus_validation_tree", "bact_validation_tree", "iqtree2", {}),
    ("prometheus.smk", "validate_prometheus", "bact_validate", "validate_prometheus.py", {}),
    ("distribution.smk", "compute_prevalence", "bact_prevalence", "compute_prevalence.py", {}),
    ("distribution.smk", "defense_matrix", "bact_defense_matrix", "build_defense_matrix.py", {}),
    ("distribution.smk", "cooccurrence", "bact_cooccurrence", "cooccurrence.R", {}),
    ("diversity.smk", "architecture_typology", "bact_architecture", "architecture_typology.py", {}),
    ("ecology.smk", "rloop_predict", "bact_rloop", "rloop_predict.py", {}),
    ("ecology.smk", "ecology_stats", "bact_ecology_stats", "ecology_stats.R",
     {"--predictors": "bacteria have no group_* columns; 'auto' = the niche_*/viral terms"}),
    ("ecology.smk", "ecology_stats", "bact_niche", "ecology_stats.R",
     {"--predictors": "as above", "--viral": "the niche control has no viral table"}),
    ("evolution.smk", "codon_align", "crossdomain_codon_align", "codon_align.py", {}),
    ("evolution.smk", "selection", "crossdomain_selection", "run_selection.py",
     {"--relax-labelling": "no RELAX Test/Reference partition across domains"}),
    ("ecology.smk", "genomad", "bact_genomad", "genomad", {}),
    ("ecology.smk", "checkv", "bact_checkv", "checkv", {}),
]


@S.test
def rules_parity_with_archaeal_calls():
    """Every option the archaeal rule passes to a shared script/tool is also passed by the
    bacterial counterpart (or exempted above with a reason). Fails when an archaeal call
    gains an option -- update crossdomain.smk (or add a reasoned exemption)."""
    missing = []
    for smk, arule, brule, tool, exempt in PARITY:
        af = _call_flags(_rule_block(smk, arule), tool)
        bf = _call_flags(_rule_block("crossdomain.smk", brule), tool)
        assert af, f"{tool} not called in {smk}:{arule} (renamed?)"
        assert bf, f"{tool} not called in crossdomain.smk:{brule}"
        gap = sorted(af - bf - set(exempt))
        if gap:
            missing.append(f"{brule} vs {smk}:{arule} ({tool}): {gap}")
    assert not missing, "bacterial calls lack archaeal options:\n  " + "\n  ".join(missing)
    iq = _call_flags(_rule_block("crossdomain.smk", "bact_validation_tree"), "iqtree2")
    assert "-fast" not in iq and {"-B", "--alrt", "--seed"} <= iq, iq     # supports required
    # the bacterial CheckV header (empty-input table + post-run check) == the archaeal one
    eco = (REPO / "workflow" / "rules" / "ecology.smk").read_text()
    m = re.search(r"CHECKV_QS_HEADER = \"\\t\"\.join\(\[(.*?)\]\)", eco, re.S)
    assert m, "CHECKV_QS_HEADER not found in ecology.smk"
    arch_hdr = re.findall(r'"([^"]+)"', m.group(1))
    cd = (REPO / "workflow" / "rules" / "crossdomain.smk").read_text()
    m2 = re.search(r"CHECKV_HEADER = \((.*?)\)\n", cd, re.S)
    bact_hdr = "".join(re.findall(r'"([^"]*)"', m2.group(1))).split("\\\\t")
    assert len(arch_hdr) == 14 and bact_hdr == arch_hdr, (bact_hdr, arch_hdr)


# --------------------------------------------------------------------------- #
# 1. extract_bacterial_prometheus.py
# --------------------------------------------------------------------------- #
def extract_inputs(d, *, drop_faa=False, universe=("g1", "g2", "g3", "g4"), no_prom=False):
    sys_rows = [
        ["g1", "g1_Prometheus_1", "Prometheus", "Prometheus", "Defense"],
        ["g1", "g1_RM_1", "RM", "RM_Type_I", "Defense"],
        ["g2", "g2_Prom_1", "PROMETHEUS", "Prometheus", "Defense"],     # case-insensitive
        ["g3", "g3_PII_1", "Prometheus_II", "Prometheus_II", "Defense"],  # not an exact match
        ["g3", "g3_ADF_1", "Anti_CBASS", "Anti_CBASS", "Antidefense"],  # dropped
        ["g3", "g3_Gab_1", "Gabija", "Gabija", "Defense"],
    ]
    gene_rows = [
        ["g1", "g1_Prometheus_1", "p1", "812.5", "Defense"],
        ["g1", "g1_RM_1", "p9", "100", "Defense"],
        ["g2", "g2_Prom_1", "WP_000002.1", "700", "Defense"],
        ["g3", "g3_PII_1", "p3", "650", "Defense"],
        ["g3", "g3_ADF_1", "p4", "50", "Antidefense"],
        ["g3", "g3_Gab_1", "p5", "90", "Defense"],
    ]
    if no_prom:
        sys_rows = [r for r in sys_rows if r[2].lower() != "prometheus"]
        gene_rows = [r for r in gene_rows if r[1] not in ("g1_Prometheus_1", "g2_Prom_1")]
    tsv(d / "systems.tsv", ["Genome", "sys_id", "type", "subtype", "activity"], sys_rows)
    tsv(d / "genes.tsv", ["genome_id", "sys_id", "hit_id", "hit_score", "activity"], gene_rows)
    faa = mk(d / "faa")
    with gzip.open(faa / "g1.faa.gz", "wt") as fh:                    # gzip-aware lookup
        fh.write(">p1 ProA-like protein\nMKKLLAAAEEE*\n>p9 RM\nMSSS\n")
    if not drop_faa:
        write(faa / "g2.faa", ">WP_000002.1 hypothetical\nMRRRAAAKK\n")
    write(faa / "g3.faa", ">p3 x\nMAAA\n>p4 y\nMCC\n>p5 z\nMDD\n")
    cds = mk(d / "cds")
    write(cds / "g1.fna", ">p1\nATGAAAAAACTGCTGGCGGCGGCGGAAGAAGAATAA\n>p9\nATGTCCTCCTCC\n")
    write(cds / "g2.fna", ">lcl|NZ_X.1_cds_WP_000002.1_1 [locus_tag=T_2] [protein_id=WP_000002.1]\n"
                          "ATGCGTCGTCGTGCGGCGGCGAAAAAATAA\n")
    tsv(d / "universe.tsv", ["genome_id"], [[g] for g in universe])


def run_extract(d, *extra, check=True):
    args = ["--df-systems", d / "systems.tsv", "--df-genes", d / "genes.tsv",
            "--faa-dir", d / "faa", "--cds-dir", d / "cds", "--systems-genome-col", "genome",
            "--genomes", d / "universe.tsv", "--out-defense-dir", d / "defense",
            "--out-candidates", d / "cand.faa", "--out-sources", d / "sources.tsv",
            "--out-raw", d / "raw.tsv", "--out-cds", d / "cds.ffn",
            "--out-report", d / "report.tsv", *extra]
    return py("extract_bacterial_prometheus.py", *args, check=check)


@S.test
def extract_exact_matching_split_and_ids():
    d = mk(ROOT / "extract1")
    extract_inputs(d)
    run_extract(d)
    cand = dict(read_fa(d / "cand.faa"))
    assert set(cand) == {"g1__p1", "g2__WP_000002.1"}, cand          # g3 Prometheus_II excluded
    assert cand["g1__p1"] == "MKKLLAAAEEE"                           # trailing '*' stripped
    cds = dict(read_fa(d / "cds.ffn"))
    assert set(cds) == set(cand)                                     # same ids as proteins
    assert cds["g2__WP_000002.1"].startswith("ATGCGT")               # via [protein_id=] tag
    src = read_tsv(d / "sources.tsv")
    assert list(src.columns) == ["protein_id", "genome_id", "detect_defensefinder", "df_score",
                                 "df_subtype"]
    assert dict(zip(src.protein_id, src.df_score)) == {"g1__p1": "812.5", "g2__WP_000002.1": "700.0"}
    for g in ("g1", "g2", "g3", "g4"):           # DefenseFinder layout for every genome
        for kind in ("systems", "genes"):
            assert (d / "defense" / g / f"{g}_defense_finder_{kind}.tsv").exists(), (g, kind)
    g3 = read_tsv(d / "defense" / "g3" / "g3_defense_finder_systems.tsv")
    assert "Genome" not in g3.columns and set(g3["type"]) == {"Prometheus_II", "Gabija"}
    g4 = read_tsv(d / "defense" / "g4" / "g4_defense_finder_systems.tsv")
    assert len(g4) == 0 and list(g4.columns) == ["sys_id", "type", "subtype", "activity"]
    rep = read_tsv(d / "report.tsv").set_index("genome_id")
    assert rep.loc["g3", "n_antidefense_dropped"] == "1" and rep.loc["g1", "n_cds_found"] == "1"


@S.test
def extract_stale_tables_removed_on_rerun():
    d = mk(ROOT / "extract_stale")
    extract_inputs(d)
    run_extract(d)
    tsv(d / "universe.tsv", ["genome_id"], [["g1"], ["g2"], ["g3"]])   # g4 leaves the universe
    run_extract(d)
    assert not (d / "defense" / "g4" / "g4_defense_finder_systems.tsv").exists()


@S.test
def extract_prodigal_genes_tables():
    """SCHEMA item-2 genes.tsv from Prodigal headers (direct fragmentation evidence for
    validate_prometheus --genes-dir), none for plain or mixed proteomes, stale removed."""
    d = mk(ROOT / "extract_genes")
    tsv(d / "systems.tsv", ["genome_id", "sys_id", "type", "subtype", "activity"],
        [[g, f"{g}_P_1", "Prometheus", "Prometheus", "Defense"] for g in ("g5", "g6", "g7")])
    tsv(d / "genes.tsv", ["genome_id", "sys_id", "hit_id", "hit_score", "activity"],
        [["g5", "g5_P_1", "ctg5_2", "900", "Defense"], ["g6", "g6_P_1", "q1", "800", "Defense"],
         ["g7", "g7_P_1", "ctg7_1", "700", "Defense"]])
    faa = mk(d / "faa")
    tail = ";start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.588"
    write(faa / "g5.faa", f">ctg5_1 # 2 # 1105 # -1 # ID=1_1;partial=00{tail}\nMKV\n"
                          f">ctg5_2 # 1200 # 1220 # 1 # ID=1_2;partial=01{tail}\nMKLLLL*\n"
                          f">ctg5b_1 # 1 # 99 # 1 # ID=2_1;partial=10{tail}\nMAA\n")
    write(faa / "g6.faa", ">q1 hypothetical protein\nMKKK\n")                 # not Prodigal
    write(faa / "g7.faa", f">ctg7_1 # 1 # 30 # 1 # ID=1_1;partial=00{tail}\nMKV\n"
                          ">plain2 some other header\nMAA\n")                     # mixed
    tsv(d / "universe.tsv", ["genome_id"], [["g5"], ["g6"], ["g7"], ["g8"]])
    write(d / "genes" / "OLD" / "OLD.genes.tsv", "protein_id\n")               # previous run
    py("extract_bacterial_prometheus.py", "--df-systems", d / "systems.tsv",
       "--df-genes", d / "genes.tsv", "--faa-dir", faa, "--genomes", d / "universe.tsv",
       "--out-defense-dir", d / "defense", "--out-genes-dir", d / "genes",
       "--out-candidates", d / "cand.faa", "--out-sources", d / "src.tsv",
       "--out-raw", d / "raw.tsv", "--out-report", d / "report.tsv")
    import rename_prodigal
    gt = read_tsv(d / "genes" / "g5" / "g5.genes.tsv")
    assert list(gt.columns) == rename_prodigal.GENES_COLS               # archaeal layout
    g = gt.set_index("protein_id")
    assert len(g) == 3 and tuple(g.loc["ctg5_2", ["contig", "contig_index", "gene_index", "start",
                                                  "end", "strand", "partial", "length_aa"]]) == \
        ("ctg5", "1", "2", "1200", "1220", "+", "01", "6")
    assert g.loc["ctg5_1", "strand"] == "-" and tuple(g.loc["ctg5b_1", ["contig", "contig_index"]]) == ("ctg5b", "2")
    assert not (d / "genes" / "g6").exists() and not (d / "genes" / "g7").exists()
    assert not (d / "genes" / "OLD" / "OLD.genes.tsv").exists()           # stale table removed
    rep = read_tsv(d / "report.tsv").set_index("genome_id")
    assert dict(rep.genes_table) == {"g5": "prodigal", "g6": "none", "g7": "none", "g8": ""}
    import validate_prometheus as vp                                     # the consumer reads it
    tbl = vp.load_genes_table(str(d / "genes"), "g5")
    assert tbl["ctg5_2"]["partial"] == "01" and tbl["ctg5b_1"]["contig_index"] == 2
    assert vp.load_genes_table(str(d / "genes"), "g6") is None           # -> length proxy


@S.test
def extract_missing_proteome_fails_unless_allowed():
    d = mk(ROOT / "extract2")
    extract_inputs(d, drop_faa=True)
    p = run_extract(d, check=False)
    assert p.returncode != 0 and "missing" in p.stderr
    run_extract(d, "--allow-missing", "True")
    rep = read_tsv(d / "report.tsv").set_index("genome_id")
    assert rep.loc["g2", "status"] == "no_proteome"
    assert set(dict(read_fa(d / "cand.faa"))) == {"g1__p1"}


@S.test
def extract_empty_and_mismatch():
    d = mk(ROOT / "extract3")
    extract_inputs(d, no_prom=True)
    run_extract(d)
    assert (d / "cand.faa").read_text() == "" and (d / "cds.ffn").read_text() == ""
    assert len(read_tsv(d / "sources.tsv")) == 0 and len(read_tsv(d / "raw.tsv")) == 0
    assert (d / "defense" / "g1" / "g1_defense_finder_systems.tsv").exists()
    d2 = mk(ROOT / "extract4")
    extract_inputs(d2, universe=("GCF_1", "GCF_2"))                 # zero id overlap
    p = run_extract(d2, check=False)
    assert p.returncode != 0 and "identical" in p.stderr
    d3 = mk(ROOT / "extract5")
    extract_inputs(d3)
    write(d3 / "systems.tsv", (d3 / "systems.tsv").read_text().replace("g1\tg1_RM_1", "g1__x\tg1_RM_1"))
    p = run_extract(d3, check=False)
    assert p.returncode != 0 and "__" in p.stderr


# --------------------------------------------------------------------------- #
# 2. build_bacterial_table.py
# --------------------------------------------------------------------------- #
CAT_HDR = ["locus_id", "genome_id", "protein_id", "length_aa", "detect_defensefinder",
           "architecture", "n_domains", "catalytic_status", "defense_context",
           "tree_assignment", "fragmented", "confidence_tier", "reason"]


def table_inputs(d):
    tsv(d / "qc.tsv", ["genome_id", "Completeness", "Contamination", "Genome_Size"], [
        ["g1", 95, 1.0, 4000000], ["g2", 96, 7.0, 3000000], ["g3", 80, 0.5, 2000000],
        ["g4", "", "", ""]])
    tsv(d / "tax.tsv", ["genome_id", "gtdb_domain", "gtdb_phylum", "gtdb_species"], [
        ["g1", "d__Bacteria", "p__Pseudomonadota", "s__X a"], ["g2", "Bacteria", "Bacillota", ""],
        ["g3", "Bacteria", "Bacillota", ""], ["g5", "Bacteria", "Bacillota", ""]])
    tsv(d / "cat.tsv", CAT_HDR, [
        ["g1__a", "g1", "g1__a", 2000, 1, "DUF4011-DNA2-REase", 3, "intact", True, "Prometheus",
         False, "A", ""],
        ["g1__b", "g1", "g1__b", 1900, 1, "DUF4011-DNA2", 2, "partial", True, "Unplaced", False,
         "C", "x"],
        ["g2__a", "g2", "g2__a", 2000, 1, "DUF4011-DNA2-REase", 3, "intact", True, "Unplaced",
         False, "B", "x"],
        ["g3__a", "g3", "g3__a", 300, 1, "DNA2", 1, "degenerate", False, "Paralog", True,
         "Rejected", "placed_with_cellular_paralog"]])
    tsv(d / "defense" / "g1" / "g1_defense_finder_systems.tsv", ["sys_id", "type", "activity"],
        [["s1", "Prometheus", "Defense"], ["s2", "RM", "Defense"], ["s3", "Gabija", "Defense"],
         ["s4", "Prometheus_II", "Defense"], ["s5", "Anti_CBASS", "Antidefense"]])
    for g in ("g2", "g3", "g4"):                                     # header-only, as extract writes
        tsv(d / "defense" / g / f"{g}_defense_finder_systems.tsv", ["sys_id", "type", "activity"], [])
    tsv(d / "samples.tsv", ["genome_id"], [["g1"], ["g2"], ["g3"], ["g4"]])


@S.test
def table_known_values():
    d = mk(ROOT / "table1")
    table_inputs(d)
    py("build_bacterial_table.py", "--taxonomy", d / "tax.tsv", "--qc", d / "qc.tsv",
       "--catalog", d / "cat.tsv", "--defense-dir", d / "defense", "--samples", d / "samples.tsv",
       "--min-completeness", 90, "--max-contamination", 5, "--out", d / "master.tsv")
    m = read_tsv(d / "master.tsv").set_index("genome_id")
    assert list(m.index) == ["g1", "g2", "g3", "g4"]                  # universe only (no g5)
    assert dict(m.pass_primary) == {"g1": "True", "g2": "False", "g3": "False", "g4": "False"}
    assert m.loc["g1", "prometheus_tier"] == "A" and m.loc["g1", "prometheus_n_loci"] == "2"
    assert m.loc["g2", "prometheus_tier"] == "B" and m.loc["g2", "prometheus_present"] == "True"
    assert m.loc["g3", "prometheus_present"] == "False" and m.loc["g3", "prometheus_tier"] == ""
    assert dict(m.detect_defensefinder) == {"g1": "1", "g2": "1", "g3": "1", "g4": "0"}
    assert m.loc["g1", "gtdb_phylum"] == "Pseudomonadota" and m.loc["g1", "gtdb_species"] == "X a"
    assert m.loc["g1", "genome_size"] == "4000000"
    # anti-defense row not counted; only the EXACT focal type is excluded (Prometheus_II stays)
    assert (m.loc["g1", "n_defense_systems"], m.loc["g1", "n_defense_systems_excl_prometheus"]) == ("4", "3")
    assert m.loc["g1", "tree_tip"] == "g1"
    # Contract v2 presence status: the B-only genome is EXCLUDED at tier A, PRESENT at A+B
    assert dict(m.prometheus_status_primary) == {"g1": "present", "g2": "excluded",
                                                 "g3": "absent", "g4": "absent"}
    assert dict(m.prometheus_status_sensitivity)["g2"] == "present"
    tsv(d / "missing_def.tsv", ["genome_id"], [["g1"], ["g2"], ["g3"], ["g4"], ["g9"]])
    write(d / "qc9.tsv", (d / "qc.tsv").read_text() + "g9\t99\t0\t1\n")
    p = fails("build_bacterial_table.py", "--taxonomy", d / "tax.tsv", "--qc", d / "qc9.tsv",
              "--catalog", d / "cat.tsv", "--defense-dir", d / "defense",
              "--samples", d / "missing_def.tsv", "--out", d / "m9.tsv")
    assert "not found" in p.stderr                                   # no silent zero defense load


@S.test
def table_genome_size_from_assembly():
    """genome_size fallback (no QC size column): streamed residue count of the assembly,
    plain or gzip-compressed (bacteria.genome_extension 'fna.gz')."""
    import build_bacterial_table as bbt
    d = mk(ROOT / "gsize")
    fasta(d / "plain.fasta", [("c1", "ACGT" * 25), ("c2", "A" * 50)])
    with gzip.open(d / "comp.fna.gz", "wt") as fh:
        fh.write(">c1 x\n" + "ACGT" * 25 + "\n>c2\n" + "G" * 30 + "\nG\n")
    assert bbt.assembly_size(str(d), "plain", "fasta") == 150
    assert bbt.assembly_size(str(d), "comp", "fna.gz") == 131
    assert bbt.assembly_size(str(d), "missing", "fasta") == ""
    assert bbt.assembly_size("", "plain", "fasta") == ""


@S.test
def table_duplicates_and_orphans_fail():
    d = mk(ROOT / "table2")
    table_inputs(d)
    write(d / "qc_dup.tsv", (d / "qc.tsv").read_text() + "g1\t90\t1\t1\n")
    p = fails("build_bacterial_table.py", "--taxonomy", d / "tax.tsv", "--qc", d / "qc_dup.tsv",
              "--catalog", d / "cat.tsv", "--defense-dir", d / "defense", "--out", d / "m.tsv")
    assert "duplicated" in p.stderr
    tsv(d / "samples2.tsv", ["genome_id"], [["g2"], ["g4"]])           # g1 catalog orphan
    p = fails("build_bacterial_table.py", "--taxonomy", d / "tax.tsv", "--qc", d / "qc.tsv",
              "--catalog", d / "cat.tsv", "--defense-dir", d / "defense",
              "--samples", d / "samples2.tsv", "--out", d / "m.tsv")
    assert "catalog genome" in p.stderr


# --------------------------------------------------------------------------- #
# 3. build_bacterial_niche_master.py
# --------------------------------------------------------------------------- #
@S.test
def niche_master_missing_is_not_absent():
    d = mk(ROOT / "niche1")
    tsv(d / "master.tsv", ["genome_id", "prometheus_present"],
        [[f"g{i}", i % 2 == 0] for i in range(1, 7)])
    tsv(d / "niche.tsv", ["genome_id", "Oxygen", "Host", "niche_soil"], [
        ["g1", "Anaerobe", "True", "1"], ["g2", "aerobe", "FALSE", "0"],
        ["g3", "anaerobe", "", ""], ["g4", "", "yes", "1"], ["g5", "", "", ""],
        ["g9", "aerobe", "no", "0"], ["g1", "facultative", "", ""]])   # g1 duplicated -> union
    py("build_bacterial_niche_master.py", "--master", d / "master.tsv", "--niche", d / "niche.tsv",
       "--out", d / "out.tsv")
    o = read_tsv(d / "out.tsv").set_index("genome_id")
    assert list(o.index) == ["g1", "g2", "g3", "g4"]                 # g5 unlabelled, g6 absent
    assert set(c for c in o.columns if c.startswith("niche_")) == {
        "niche_oxygen_aerobe", "niche_oxygen_anaerobe", "niche_oxygen_facultative",
        "niche_host", "niche_soil"}                                  # binary case-insensitive
    assert o.loc["g1", "niche_oxygen_anaerobe"] == "1" and o.loc["g1", "niche_oxygen_facultative"] == "1"
    assert o.loc["g4", "niche_oxygen_anaerobe"] == ""                 # missing stays missing
    assert o.loc["g3", "niche_host"] == "" and o.loc["g4", "niche_host"] == "1"
    assert o.loc["g2", "niche_host"] == "0"


@S.test
def niche_master_failures():
    d = mk(ROOT / "niche2")
    tsv(d / "master.tsv", ["genome_id"], [["g1"], ["g2"]])
    p = fails("build_bacterial_niche_master.py", "--master", d / "master.tsv",
              "--niche", d / "missing.tsv", "--out", d / "o.tsv")
    assert "does not exist" in p.stderr
    tsv(d / "niche.tsv", ["genome_id", "oxygen"], [["X1", "aerobe"], ["X2", "anaerobe"]])
    p = fails("build_bacterial_niche_master.py", "--master", d / "master.tsv",
              "--niche", d / "niche.tsv", "--out", d / "o.tsv")
    assert "identical" in p.stderr


# --------------------------------------------------------------------------- #
# 4. prune_tree_to_taxa.py
# --------------------------------------------------------------------------- #
GTDB_TREE = ("(((RS_GCF_000000001.2:1,GB_GCA_000000002.1:1)'100.0:g__X':0.5,"
             "(RS_GCF_000000003.1:1,RS_GCF_000000004.1:1)'99.0:g__Y; s__Y z':0.5)'100.0:p__A':0.2,"
             "((RS_GCF_000000005.1:1,RS_GCF_000000006.1:1)'98:g__Z':0.4,"
             "RS_GCF_000000007.1:1.5)'100.0:p__B':0.2)'100.0:d__Bacteria';\n")


@S.test
def prune_normalisation_collapse_and_rooting():
    d = mk(ROOT / "prune1")
    write(d / "gtdb.tree", GTDB_TREE)
    tsv(d / "keep.tsv", ["genome_id", "completeness", "contamination"], [
        ["GCF_000000001.2", 99, 0], ["GCF_000000002.1", 98, 0], ["GCF_000000003", 97, 0],
        ["NR_a", 95, 2], ["NR_b", 99, 0], ["GCF_000999999.1", 99, 0], ["RS_GCF_000000007.1", 90, 1]])
    tsv(d / "tax.tsv", ["genome_id", "classification"], [
        ["NR_a", "d__Bacteria;p__B;g__Z;s__Zeta one"], ["NR_b", "d__Bacteria;p__B;g__Z;s__Zeta one"]])
    write(d / "reftax.tsv", "RS_GCF_000000005.1\td__Bacteria;p__B;g__Z;s__Zeta one\n"
                            "RS_GCF_000000006.1\td__Bacteria;p__B;g__Z;s__Zeta two\n")
    py("prune_tree_to_taxa.py", "--in-tree", d / "gtdb.tree", "--keep", d / "keep.tsv",
       "--quality", d / "keep.tsv", "--taxonomy", d / "tax.tsv", "--ref-taxonomy", d / "reftax.tsv",
       "--min-match-fraction", 0.5, "--out-tree", d / "out.nwk", "--out-report", d / "rep.tsv")
    rep = read_tsv(d / "rep.tsv").set_index("genome_id")
    assert dict(rep["match"]) == {
        "GCF_000000001.2": "prefix_stripped", "GCF_000000002.1": "paired_accession",
        "GCF_000000003": "version_stripped", "NR_a": "collapsed", "NR_b": "species_representative",
        "GCF_000999999.1": "unmatched", "RS_GCF_000000007.1": "exact"}
    assert rep.loc["NR_b", "tree_tip"] == "RS_GCF_000000005.1"
    txt = (d / "out.nwk").read_text()
    assert "'" not in txt and "g__" not in txt                       # internal labels dropped
    assert tips_of(txt) == sorted(["GCF_000000001.2", "GCF_000000002.1", "GCF_000000003", "NR_b",
                                   "RS_GCF_000000007.1"])            # relabelled to genome_id
    r = run([RSCRIPT, "-e", f"t <- ape::read.tree('{d / 'out.nwk'}'); "
                            "cat(ape::Ntip(t), ape::is.rooted(t), '\\n')"])
    assert r.stdout.split() == ["5", "TRUE"], r.stdout                # rooting preserved


@S.test
def prune_low_match_fails_and_empty_is_trivial():
    d = mk(ROOT / "prune2")
    write(d / "gtdb.tree", GTDB_TREE)
    tsv(d / "keep.tsv", ["genome_id"], [["GCF_000000001.2"], ["foo"], ["bar"]])
    p = fails("prune_tree_to_taxa.py", "--in-tree", d / "gtdb.tree", "--keep", d / "keep.tsv",
              "--out-tree", d / "o.nwk", "--out-report", d / "r.tsv")
    assert "1/3" in p.stderr and "minimum" in p.stderr
    tsv(d / "empty.tsv", ["genome_id"], [])
    py("prune_tree_to_taxa.py", "--in-tree", d / "gtdb.tree", "--keep", d / "empty.tsv",
       "--out-tree", d / "e.nwk")
    assert (d / "e.nwk").read_text().strip() == "();"
    tsv(d / "bad.tsv", ["genome_id"], [["has space"]])
    fails("prune_tree_to_taxa.py", "--in-tree", d / "gtdb.tree", "--keep", d / "bad.tsv",
          "--out-tree", d / "b.nwk")


@S.test
def prune_gene_tree_exact_mode():
    d = mk(ROOT / "prune3")
    write(d / "gene.nwk", "(a__1:0.1,(b__2:0.2,c__3:0.3)90/95:0.1,(R1:0.4,P1:0.5)100/100:0.2);\n")
    fasta(d / "codon.fa", [("a__1", "ATG"), ("b__2", "ATG"), ("c__3", "ATG")])
    py("prune_tree_to_taxa.py", "--in-tree", d / "gene.nwk", "--keep-fasta", d / "codon.fa",
       "--match", "exact", "--no-relabel", "--min-match-fraction", 1.0, "--out-tree", d / "o.nwk")
    assert tips_of((d / "o.nwk").read_text()) == ["a__1", "b__2", "c__3"]
    fasta(d / "codon2.fa", [("a__1", "ATG"), ("zz__9", "ATG")])
    fails("prune_tree_to_taxa.py", "--in-tree", d / "gene.nwk", "--keep-fasta", d / "codon2.fa",
          "--match", "exact", "--no-relabel", "--min-match-fraction", 1.0, "--out-tree", d / "o2.nwk")


# --------------------------------------------------------------------------- #
# 5. crossdomain_analysis.py (hand-built trees with known parsimony answers)
# --------------------------------------------------------------------------- #
T_NESTED = ("((O1:0.1,O2:0.1)100/100:0.5,(((((((A1:0.1,A2:0.1)100/100:0.1,B1:0.1)100/100:0.1,"
            "B2:0.1)100/100:0.1,B3:0.1)100/100:0.1,B4:0.1)100/100:0.1,B5:0.1)100/100:0.1,"
            "(B6:0.1,B7:0.1)100/100:0.1)100/100:0.5);")
T_SISTER = ("((O1:0.1,O2:0.1)100/100:0.5,((A1:0.1,A2:0.1)100/100:0.2,(((B1:0.1,B2:0.1)100/100:0.1,"
            "(B3:0.1,B4:0.1)100/100:0.1)100/100:0.1,(B5:0.1,(B6:0.1,B7:0.1)100/100:0.1)100/100:0.1)"
            "100/100:0.2)100/100:0.5);")
T_ARCHROOT = ("((O1:0.1,O2:0.1)100/100:0.5,(((((((B1:0.1,B2:0.1)100/100:0.1,A1:0.1)100/100:0.1,"
              "A2:0.1)100/100:0.1,A3:0.1)100/100:0.1,A4:0.1)100/100:0.1,A5:0.1)100/100:0.1,"
              "(A6:0.1,A7:0.1)100/100:0.1)100/100:0.5);")
T_COLLAPSE = ("((O1:0.1,O2:0.1)100/100:0.5,(((A1:0.1,(B1:0.1,B2:0.1)100/100:0.1)50/40:0.1,"
              "A2:0.1)100/100:0.1,(B3:0.1,B4:0.1)100/100:0.1)100/100:0.5);")
T_NONMONO = ("((O1:0.1,O2:0.1)100/100:0.5,((B1:0.1,(O3:0.1,A1:0.1)90/90:0.1)90/90:0.1,"
             "(B2:0.1,B3:0.1)100/100:0.1)100/100:0.5);")
# bacterial ML root that rests on ONE bacterial sister lineage (B1): thinning the
# bacteria to 5 drops B1 in 2/7 of subsamples and the root becomes A/B
T_FRAGILE = ("((O1:0.1,O2:0.1)100/100:0.5,(((A1:0.1,A2:0.1)100/100:0.1,B1:0.1)100/100:0.2,"
             "(((B2:0.1,B3:0.1)100/100:0.1,(B4:0.1,B5:0.1)100/100:0.1)100/100:0.1,"
             "(B6:0.1,B7:0.1)100/100:0.1)100/100:0.2)100/100:0.5);")


def analysis_inputs(d, n_arch=7, n_bact=7, outgroup=("O1", "O2")):
    fasta(d / "arch.faa", [(f"A{i}", "M") for i in range(1, n_arch + 1)])
    fasta(d / "bact.faa", [(f"B{i}", "M") for i in range(1, n_bact + 1)])
    tsv(d / "labels.tsv", ["id", "group"], [[o, "Paralog_ref"] for o in outgroup])


def run_analysis(d, tree, *extra, check=True, tag="s"):
    write(d / f"{tag}.nwk", tree + "\n")
    return py("crossdomain_analysis.py", "--tree", d / f"{tag}.nwk", "--arch-proteins", d / "arch.faa",
              "--bact-proteins", d / "bact.faa", "--ref-labels", d / "labels.tsv",
              "--out", d / f"{tag}.tsv", "--out-tips", d / f"{tag}.tips.tsv",
              "--out-clades", d / f"{tag}.clades.tsv", "--out-rooted-tree", d / f"{tag}.rooted.nwk",
              "--n-subsample", 200, "--seed", 7, *extra, check=check)


@S.test
def analysis_nested_clade_is_robust_bacterial_origin():
    d = mk(ROOT / "an1")
    analysis_inputs(d, n_arch=2)
    run_analysis(d, T_NESTED)
    s = summary(d / "s.tsv")
    assert s["status"] == "complete" and s["rooting"] == "outgroup"
    assert s["ingroup_root_state"] == "B" and s["collapsed_root_state"] == "B"
    assert s["n_crossdomain_transitions"] == "1" and s["n_archaeal_clades"] == "1"
    assert s["archaeal_monophyletic"] == "True" and s["n_archaeal_clades_nested_supported"] == "1"
    assert s["subsample_root_B_fraction"] == "1.0", s["subsample_root_B_fraction"]
    assert s["inferred_direction"] == "bacteria->archaea", s
    assert s["direction_method"] == "topology-based"
    c = read_tsv(d / "s.clades.tsv")
    assert len(c) == 1 and c.loc[0, "nested_in_bacteria"] == "True" and c.loc[0, "n_tips"] == "2"
    assert c.loc[0, "stem_support"] == "100/100" and c.loc[0, "parent_support"] == "100/100"
    tips = read_tsv(d / "s.tips.tsv").set_index("tip")
    assert tips.loc["A1", "archaeal_clade"] == "AC1" and tips.loc["O1", "domain"] == "Paralog"


@S.test
def analysis_sister_clade_is_ambiguous():
    d = mk(ROOT / "an2")
    analysis_inputs(d, n_arch=2)
    run_analysis(d, T_SISTER)
    s = summary(d / "s.tsv")
    assert s["ingroup_root_state"] == "A/B" and s["inferred_direction"] == "ambiguous"
    assert s["inferred_direction_ml"] == "ambiguous" and s["n_crossdomain_transitions"] == "1"
    assert read_tsv(d / "s.clades.tsv").loc[0, "nested_in_bacteria"] == "False"


@S.test
def analysis_archaeal_root_and_imbalance_control():
    d = mk(ROOT / "an3")
    analysis_inputs(d, n_arch=7, n_bact=2)
    run_analysis(d, T_ARCHROOT)
    s = summary(d / "s.tsv")
    assert s["ingroup_root_state"] == "A" and s["inferred_direction"] == "archaea->bacteria", s
    # the ML-tree root is bacterial only because of one sister lineage: the imbalance
    # control (bacteria thinned to 5 of 7) loses it in ~2/7 of subsamples -> no call
    d2 = mk(ROOT / "an3b")
    analysis_inputs(d2, n_arch=2)
    run_analysis(d2, T_FRAGILE)
    s2 = summary(d2 / "s.tsv")
    assert s2["inferred_direction_ml"] == "bacteria->archaea" and s2["ingroup_root_state"] == "B"
    assert 0.55 < float(s2["subsample_root_B_fraction"]) < 0.85, s2["subsample_root_B_fraction"]
    note = read_tsv(d2 / "s.tsv").set_index("metric").loc["inferred_direction", "note"]
    assert s2["inferred_direction"] == "ambiguous" and "balanced subsamples" in note


@S.test
def analysis_support_collapse_counts():
    d = mk(ROOT / "an4")
    analysis_inputs(d, n_arch=2, n_bact=4)
    run_analysis(d, T_COLLAPSE)
    s = summary(d / "s.tsv")
    assert (s["n_archaeal_clades"], s["n_archaeal_clades_supported"]) == ("2", "1")
    assert (s["n_crossdomain_transitions"], s["n_crossdomain_transitions_supported"]) == ("2", "1")
    assert (s["ingroup_root_state"], s["collapsed_root_state"]) == ("A/B", "B")
    assert s["inferred_direction"] == "ambiguous"


@S.test
def analysis_nonmonophyletic_outgroup_unknown_and_trivial():
    d = mk(ROOT / "an5")
    analysis_inputs(d, n_arch=1, n_bact=3, outgroup=("O1", "O2", "O3"))
    run_analysis(d, T_NONMONO)
    s = summary(d / "s.tsv")
    assert s["outgroup_monophyletic"] == "False" and s["rooting"].startswith("outgroup_nonmono")
    assert s["inferred_direction"] == "ambiguous"
    d2 = mk(ROOT / "an6")
    analysis_inputs(d2, n_arch=2, n_bact=6)                          # B7 unknown to the sets
    p = run_analysis(d2, T_NESTED, check=False)
    assert p.returncode != 0 and "none of the input sets" in p.stderr
    run_analysis(d2, "();", tag="t")
    s = summary(d2 / "t.tsv")
    assert s["status"] == "incomplete" and (d2 / "t.rooted.nwk").read_text().strip() == "();"


@S.test
def analysis_bootstrap_direction_support():
    d = mk(ROOT / "an7")
    analysis_inputs(d, n_arch=2)

    def strip(t):
        return re.sub(r"\)[0-9./]+:", "):", t)
    write(d / "boot9.trees", "\n".join([strip(T_NESTED)] * 9 + [strip(T_SISTER)]) + "\n")
    run_analysis(d, T_NESTED, "--ufboot", d / "boot9.trees", tag="b9")
    s = summary(d / "b9.tsv")
    assert (s["n_bootstrap_trees"], s["bootstrap_root_B_fraction"]) == ("10", "0.9")
    assert s["inferred_direction"] == "ambiguous"
    write(d / "boot10.trees", "\n".join([strip(T_NESTED)] * 10) + "\n")
    run_analysis(d, T_NESTED, "--ufboot", d / "boot10.trees", tag="b10")
    s = summary(d / "b10.tsv")
    assert s["bootstrap_root_B_fraction"] == "1.0" and s["inferred_direction"] == "bacteria->archaea"


@S.test
def analysis_fitch_rule_and_parser():
    import crossdomain_analysis as ca
    A, B, AB = frozenset("A"), frozenset("B"), frozenset("AB")
    assert ca.combine([(A, 0), (B, 0)]) == (AB, 1)
    assert ca.combine([(A, 0), (AB, 0)]) == (A, 0)
    assert ca.combine([(A, 1), (A, 0), (B, 2)]) == (AB, 4)            # soft polytomy: one change
    assert ca.combine([(AB, 0), (AB, 0), (frozenset(), 0)]) == (AB, 0)
    t = ca.UTree("((a:1,'b c':2)'90/95':0.5,(d:1,e:1)80/70:1)root:0;")
    assert sorted(t.by_label) == ["a", "b c", "d", "e"]
    # the two root edges are one bipartition: one edge, the more conservative support
    assert list(t.support.values()) == [(80.0, 70.0)]
    t2 = ca.UTree("(a:1,(b:1,c:1)90/95:0.5,(d:1,e:1)80/70:1);")       # IQ-TREE trifurcation
    assert sorted(t2.support.values()) == [(80.0, 70.0), (90.0, 95.0)]
    assert ca.parse_support("100.0:p__X") == (None, 100.0) and ca.parse_support("ab") is None


@S.test
def analysis_saturation_ng86():
    import crossdomain_analysis as ca
    d = mk(ROOT / "an8")
    rng = random.Random(3)
    a1 = "CTT" * 100
    a2 = "CTT" * 90 + "CTC" * 10          # 10 synonymous differences over 100 syn. sites
    b1, b2 = rand_codons(300, rng), rand_codons(300, rng)
    fasta(d / "codon.fa", [("A1", a1), ("A2", a2), ("B1", b1), ("B2", b2)])
    rows = {r[0]: r for r in ca.saturation(str(d / "codon.fa"),
                                           {"A1": "Archaea", "A2": "Archaea", "B1": "Bacteria",
                                            "B2": "Bacteria"}, 50, 1, 3.0)}
    expect = -0.75 * math.log(1 - 4.0 / 3.0 * 0.1)                  # NG86 + Jukes-Cantor
    assert abs(rows["within_archaea"][4] - expect) < 1e-3, rows["within_archaea"]
    assert rows["within_archaea"][2] == 0 and rows["within_archaea"][5] == 0.0
    assert rows["cross_domain"][3] == 1.0, rows["cross_domain"]       # random codons saturate
    empty = ca.saturation("", {}, 10, 1, 3.0)
    assert empty[0][0] == "all" and "not assessable" in empty[0][-1]


# --------------------------------------------------------------------------- #
# 6. plots (Contract v2 presence rule, spec-derived categories, formats)
# --------------------------------------------------------------------------- #
MASTER_HDR = ["genome_id", "pass_primary", "prometheus_present", "prometheus_tier"]


@S.test
def plot_comparison_presence_rule_regression():
    d = mk(ROOT / "plot1")
    tsv(d / "arch.tsv", MASTER_HDR, [["a1", True, True, "A"], ["a2", True, True, "B"],
                                     ["a3", True, False, ""], ["a4", True, False, ""],
                                     ["a5", False, True, "A"]])
    tsv(d / "bact.tsv", MASTER_HDR, [["b1", True, True, "C"], ["b2", True, False, ""]])
    typ_hdr = ["architecture", "n_loci", "example_locus", "has_DUF4011", "has_DUF3320",
               "has_DNA2", "has_nuclease", "has_REase"]
    tsv(d / "atyp.tsv", typ_hdr, [["DUF4011-DNA2-REase", 3, "x", True, False, True, False, True],
                                  ["DUF4011-DNA2-REase-DUF3320", 1, "y", True, True, True, False, True]])
    tsv(d / "btyp.tsv", typ_hdr[:3], [["DUF4011-DNA2", 2, "z"]])     # no has_* -> tokens
    for tiers in ("A", "A,B"):
        out = d / tiers.replace(",", "")
        py("plot_crossdomain_extra.py", "--arch-master", d / "arch.tsv", "--bact-master", d / "bact.tsv",
           "--arch-typology", d / "atyp.tsv", "--bact-typology", d / "btyp.tsv", "--outdir", out,
           "--domains", json.dumps(SPEC), "--primary-tiers", tiers, "--formats", "png,svg,pdf",
           "--out-table", out / "t.tsv")
        for ext in ("png", "svg", "pdf"):
            assert (out / f"crossdomain_comparison.{ext}").stat().st_size > 1000
        t = read_tsv(out / "t.tsv")
        prev = t[t.panel == "prevalence"].set_index("domain_of_life")
        arch = t[(t.panel == "architecture") & (t.domain_of_life == "Archaea")].set_index("category")
        assert list(arch.index) == ["DUF4011", "DNA2", "REase", "DUF3320"]   # spec order, no 'nuclease'
        assert float(arch.loc["DUF3320", "value"]) == 0.25
        bacty = t[(t.panel == "architecture") & (t.domain_of_life == "Bacteria")].set_index("category")
        assert float(bacty.loc["DNA2", "value"]) == 1.0 and float(bacty.loc["REase", "value"]) == 0.0
        if tiers == "A":   # the B-only genome is EXCLUDED, not absent; a5 fails pass_primary
            assert (prev.loc["Archaea", "n"], prev.loc["Archaea", "k"],
                    prev.loc["Archaea", "n_excluded"]) == ("3", "1", "1")
            assert abs(float(prev.loc["Archaea", "value"]) - 1 / 3) < 1e-6
            assert (prev.loc["Bacteria", "n"], prev.loc["Bacteria", "n_excluded"]) == ("1", "1")
        else:              # ... and PRESENT under A+B
            assert (prev.loc["Archaea", "n"], prev.loc["Archaea", "k"],
                    prev.loc["Archaea", "n_excluded"]) == ("4", "2", "0")


@S.test
def plot_viral_presence_and_undefined_fraction():
    d = mk(ROOT / "plot2")
    rows_m, rows_v = [], []
    for i in range(24):
        tier = ["A", "", "B"][i % 3]
        rows_m.append([f"g{i}", True, tier != "", tier])
        nprov = 0 if i in (0, 1) else 2
        rows_v.append([f"g{i}", nprov, 10.0, 90.0, 1 if i % 2 else 0, (0.5 if i % 2 else 0.0)])
    vh = ["genome_id", "n_proviruses", "total_provirus_kb", "mean_checkv_completeness",
          "n_rloop_proviruses", "frac_rloop_proviruses"]
    tsv(d / "m.tsv", MASTER_HDR, rows_m)
    tsv(d / "v.tsv", vh, rows_v)
    tsv(d / "assoc.tsv", ["predictor", "level", "model", "estimate", "std_error", "p_value",
                          "p_adj", "n", "note"],
        [["frac_rloop_proviruses", "continuous", "phyloglm", 1.5, 0.3, 0.01, 0.02, 12, ""]])
    py("plot_viral_crossdomain.py", "--arch-master", d / "m.tsv", "--arch-viral", d / "v.tsv",
       "--bact-master", d / "m.tsv", "--bact-viral", d / "v.tsv", "--bact-assoc", d / "assoc.tsv",
       "--outdir", d / "out", "--primary-tiers", "A", "--out-table", d / "out" / "t.tsv")
    assert (d / "out" / "viral_rloop_crossdomain.png").stat().st_size > 1000
    assert "p_adj=0.02" in (d / "out" / "viral_rloop_crossdomain.svg").read_text()
    t = read_tsv(d / "out" / "t.tsv")
    # 24 genomes - 8 tier-B (excluded) - g0/g1 without proviruses (0/0 undefined) = 14
    assert sum(int(x) for x in t[t.domain_of_life == "Archaea"]["n"]) == 14, t
    py("plot_viral_crossdomain.py", "--arch-master", d / "m.tsv", "--arch-viral", d / "none.tsv",
       "--bact-master", d / "m.tsv", "--bact-viral", d / "none.tsv", "--outdir", d / "empty")
    assert (d / "empty" / "viral_rloop_crossdomain.svg").exists()


@S.test
def plot_tree_figure_and_placeholder():
    d = mk(ROOT / "plot3")
    analysis_inputs(d, n_arch=2)
    run_analysis(d, T_NESTED)
    py("plot_crossdomain_tree.py", "--tree", d / "s.rooted.nwk", "--tips", d / "s.tips.tsv",
       "--summary", d / "s.tsv", "--out-prefix", d / "fig" / "tree", "--formats", "png,svg")
    assert (d / "fig" / "tree.png").stat().st_size > 2000
    svg = (d / "fig" / "tree.svg").read_text()
    assert "<text" in svg and "bacteria-&gt;archaea" in svg            # editable text + title
    run_analysis(d, "();", tag="t")
    py("plot_crossdomain_tree.py", "--tree", d / "t.rooted.nwk", "--tips", d / "t.tips.tsv",
       "--summary", d / "t.tsv", "--out-prefix", d / "fig" / "empty")
    assert (d / "fig" / "empty.png").exists()
    write(d / "other.nwk", "((X1:1,X2:1):1,(X3:1,X4:1):1);\n")
    fails("plot_crossdomain_tree.py", "--tree", d / "other.nwk", "--tips", d / "s.tips.tsv",
          "--summary", d / "s.tsv", "--out-prefix", d / "fig" / "bad")


# --------------------------------------------------------------------------- #
# 7. end-to-end Snakemake runs of crossdomain.smk
# --------------------------------------------------------------------------- #
ARCH_SEQS = ["tr|A0A1H8VVA9|A0A1H8VVA9_9EURY", "tr|A0A2A2F1Z1|A0A2A2F1Z1_9EURY"]
BACT_SEQS = ["tr|A0A060RCK4|A0A060RCK4_9BACT", "tr|A0A073IP09|A0A073IP09_9BACT",
             "tr|A0A0A6X990|A0A0A6X990_ACTUT", "tr|A0A0A7PKR5|A0A0A7PKR5_9SPHN"]
REF_SEQ = "tr|A0ABR6IN60|A0ABR6IN60_9HYPH"
PARA_SEQS = ["sp|P51530|DNA2_HUMAN", "sp|P38859|DNA2_YEAST", "sp|Q6ZQJ5|DNA2_MOUSE"]
BG = [f"GCF_0000100{i:02d}.1" for i in range(1, 9)]
STUB_SNAKEFILE = r'''
include: "rules/common.smk"

# stubs for the archaeal producers: the test writes their outputs; running one is an error
rule validate_prometheus:
    output:
        catalog=f"{RESULTS}/prometheus/prometheus_catalog.tsv",
        proteins=f"{RESULTS}/prometheus/prometheus_proteins.faa",
        sensitivity=f"{RESULTS}/prometheus/detection_sensitivity.tsv",
    shell: "echo 'stub rule must not run' >&2; exit 1"
rule combined_cds:
    output: ffn=f"{RESULTS}/evolution/all_cds.ffn"
    shell: "exit 1"
rule build_master_table:
    output: master=f"{RESULTS}/scaffolds/master_table.tsv"
    shell: "exit 1"
rule architecture_typology:
    output: typology=f"{RESULTS}/diversity/architecture_typology.tsv"
    shell: "exit 1"
rule prometheus_domain_scan:
    output: domtbl=f"{RESULTS}/prometheus/domains.domtbl"
    shell: "exit 1"
rule viral_summary:
    output: summary=f"{RESULTS}/ecology/viral/viral_summary.tsv"
    shell: "exit 1"

include: "rules/crossdomain.smk"

rule cd_e2e:
    input: crossdomain_targets()
'''
FAKE_GENOMAD = r'''#!{py}
# geNomad 1.8.0 end-to-end CLI (genomad/cli.py:350-418): INPUT OUTPUT DATABASE; flags
# --cleanup --restart --verbose; --threads INT; --splits INT >= 0. Unknown options fail.
import os, sys
a = sys.argv[1:]
VALUE, FLAGS = {"--threads", "--splits"}, {"--cleanup", "--restart", "--verbose"}
pos, opts, i = [], set(), 0
while i < len(a):
    if a[i] in VALUE:
        if int(a[i + 1]) < 0: sys.exit(f"Error: Invalid value for {a[i]}")
        opts.add(a[i]); i += 2; continue
    if a[i] in FLAGS: opts.add(a[i]); i += 1; continue
    if a[i].startswith("-"): sys.exit(f"Error: No such option: {a[i]}")
    pos.append(a[i]); i += 1
assert pos[0] == "end-to-end" and len(pos) == 4, pos
inp, out, db = pos[1:4]
if not os.path.exists(db): sys.exit(f"Error: Path '{db}' does not exist.")
stem = os.path.basename(inp).rsplit(".", 1)[0]
sd = os.path.join(out, stem + "_summary"); os.makedirs(sd, exist_ok=True)
viral = int(stem.split(".")[0][-1]) % 2 == 1
hdr = "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
seq = ("GGGGAGGGGTGGGGCAGGGGTCCGGGGAGGGG" + "ATCGATCGTAGCTAGCTAGG" * 10) * 4
with open(os.path.join(sd, stem + "_virus.fna"), "w") as fh:
    if viral: fh.write(">contig1|provirus_1_" + str(len(seq)) + "\n" + seq + "\n")
with open(os.path.join(sd, stem + "_virus_summary.tsv"), "w") as fh:
    fh.write(hdr)
    if viral: fh.write("contig1|provirus_1_%d\t%d\tProvirus\t1-%d\t3\t11\t0.95\tNA\t1\t2.0\tUnclassified\n" % (len(seq), len(seq), len(seq)))
'''
FAKE_CHECKV = r'''#!{py}
# CheckV 1.0.3 end_to_end CLI (checkv/modules/end_to_end.py): input output; -d DB -t N;
# flags --restart --remove_tmp --quiet. Unknown options fail.
import os, sys
a = sys.argv[1:]
assert a[0] == "end_to_end", a
VALUE, FLAGS = {"-d", "-t"}, {"--restart", "--remove_tmp", "--quiet"}
pos, val, i = [], {}, 1
while i < len(a):
    if a[i] in VALUE: val[a[i]] = a[i + 1]; i += 2; continue
    if a[i] in FLAGS: i += 1; continue
    if a[i].startswith("-"): sys.exit(f"checkv end_to_end: error: unrecognized arguments: {a[i]}")
    pos.append(a[i]); i += 1
inp, out = pos[0], pos[1]
db = val["-d"]
if not os.path.isdir(os.path.join(db, "genome_db")): sys.exit(f"checkv: database {db} not found")
os.makedirs(out, exist_ok=True)
names = [l[1:].split()[0] for l in open(inp) if l.startswith(">")]
hdr = ["contig_id", "contig_length", "provirus", "proviral_length", "gene_count", "viral_genes",
       "host_genes", "checkv_quality", "miuvig_quality", "completeness", "completeness_method",
       "contamination", "kmer_freq", "warnings"]
with open(os.path.join(out, "quality_summary.tsv"), "w") as fh:
    fh.write("\t".join(hdr) + "\n")
    for n in names:
        fh.write("\t".join([n, "560", "Yes", "560", "3", "2", "0", "Medium-quality", "Genome-fragment", "60.0", "AAI-based", "0.0", "1.0", ""]) + "\n")
'''
FAIL_TOOL = "#!/bin/sh\necho 'injected failure' >&2\nexit 3\n"
# CheckV 1.0.3 exits 0 without results on some fatal errors (e.g. duplicate FASTA ids)
CHECKV_SILENT = "#!/bin/sh\necho 'Error: duplicate sequence ids (simulated)'\nexit 0\n"
CHECKV_BADHDR = ("#!/bin/sh\nout=$3\nmkdir -p $out\n"
                 "printf 'contig_id\\tcontig_length\\tcompleteness\\tcheckv_quality\\n' > $out/quality_summary.tsv\n")


def duf_starts():
    out = {}
    for line in (FIXTURES / "proA_homologs.pfam_subset.domtbl").read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if f[0] == "DUF4011":
            out.setdefault(f[3], int(f[17]))
    return out


def build_e2e(work: Path, empty=False):
    """Fixture-derived inputs + stub archaeal outputs + config for one scenario. The
    proteins are 650-aa DUF4011-anchored windows of the real ProA homologues (and the
    C-terminal 650 aa of three DNA2s as the outgroup) so IQ-TREE stays fast."""
    rng = random.Random(11)
    shutil.copytree(REPO / "workflow", work / "workflow",
                    ignore=shutil.ignore_patterns("__pycache__", ".snakemake"))
    shutil.copytree(REPO / "config", work / "config")
    write(work / "workflow" / "Snakefile.cdtest", STUB_SNAKEFILE)
    # Snakemake 8.30 hashes every rule's conda env for provenance even when environments
    # are ignored (it runs `conda info` / `conda --version`); this sandbox has no conda
    # frontend, so the PRIVATE copy drops the `conda:` lines. The shell logic under test
    # is untouched (tools come from PROM_TOOLS_BIN / PROM_PYTHON / PROM_RSCRIPT).
    rf = work / "workflow" / "rules" / "crossdomain.smk"
    lines = rf.read_text().splitlines(keepends=True)
    kept = [ln for ln in lines if not re.match(r'^\s+conda:\s*"\.\./envs/[a-z_]+\.yaml"\s*$', ln)]
    assert len(lines) - len(kept) >= 20, "expected to strip the conda: directives"
    rf.write_text("".join(kept))
    homo = dict(read_fa(FIXTURES / "proA_homologs.faa"))
    para = dict(read_fa(FIXTURES / "paralogs.faa"))
    starts = duf_starts()
    prot = {n: homo[n][max(0, starts[n] - 20):max(0, starts[n] - 20) + 650]
            for n in ARCH_SEQS + BACT_SEQS + [REF_SEQ]}
    prot.update({n: para[n][-650:] for n in PARA_SEQS})
    res = work / "results"
    # ---- archaeal stubs --------------------------------------------------------- #
    agen = [f"ARC{i}" for i in range(1, 7)]
    tsv(work / "arch" / "species_metadata.tsv", ["genome_id"], [[g] for g in agen])
    a_ids = ["ARC1__ARC1_1_5", "ARC2__ARC2_2_7"]
    fasta(res / "prometheus" / "prometheus_proteins.faa",
          list(zip(a_ids, [prot[n] for n in ARCH_SEQS])))
    tsv(res / "prometheus" / "prometheus_catalog.tsv", CAT_HDR, [
        [a_ids[0], "ARC1", a_ids[0], 650, 1, "DUF4011-DNA2", 2, "intact", True, "Prometheus", False, "A", ""],
        [a_ids[1], "ARC2", a_ids[1], 650, 1, "DUF4011-DNA2", 2, "intact", True, "Prometheus", False, "B", ""]])
    tsv(res / "prometheus" / "detection_sensitivity.tsv",
        ["genome_id", "detect_defensefinder", "prometheus_present"], [])
    fasta(res / "evolution" / "all_cds.ffn",
          [(i, back_translate(prot[n], rng)) for i, n in zip(a_ids, ARCH_SEQS)] +
          [("ARC3__ARC3_1_1", back_translate("MKKLAAV", rng))])
    mhdr = ["genome_id", "representative", "tree_tip", "gtdb_domain", "gtdb_phylum", "completeness",
            "contamination", "pass_primary", "genome_size", "methanogen", "prometheus_present",
            "prometheus_tier", "prometheus_n_loci", "n_defense_systems_excl_prometheus"]
    tsv(res / "scaffolds" / "master_table.tsv", mhdr,
        [[g, True, g, "Archaea", "Halobacteriota", 95, 1, True, 3000000, False, g in ("ARC1", "ARC2"),
          {"ARC1": "A", "ARC2": "B"}.get(g, ""), int(g in ("ARC1", "ARC2")), 2] for g in agen])
    tsv(res / "diversity" / "architecture_typology.tsv",
        ["architecture", "n_loci", "example_locus", "has_DUF4011", "has_DUF3320", "has_DNA2",
         "has_nuclease", "has_REase"], [["DUF4011-DNA2", 2, a_ids[0], True, False, True, False, False]])
    write(res / "prometheus" / "domains.domtbl", "# stub\n")
    tsv(res / "ecology" / "viral" / "viral_summary.tsv",
        ["genome_id", "n_proviruses", "total_provirus_kb", "mean_checkv_completeness",
         "n_rloop_proviruses", "frac_rloop_proviruses"],
        [[g, i % 3, 10 * (i % 3), 70, i % 2, 0.5 * (i % 2)] for i, g in enumerate(agen)])
    # ---- bacterial inputs ------------------------------------------------------- #
    bi = work / "bact"
    sys_rows, gene_rows = [], []
    faa, cds = mk(bi / "faa"), mk(bi / "cds")
    prom_hits = {}
    for k, g in enumerate(BG):
        decoys = [(f"{g}_1_{j}", "M" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(60)))
                  for j in (1, 2)]
        recs = list(decoys)
        if k < 4 and not empty:
            hid = f"{g}_1_3" if k < 3 else "WP_000000004.1"
            prom_hits[g] = hid
            recs.append((hid, prot[BACT_SEQS[k]]))
            sys_rows.append([g, f"{g}_Prometheus_1", "Prometheus", "Prometheus", "Defense"])
            gene_rows.append([g, f"{g}_1", hid, "ProA", 3, "defense-finder-models/Prometheus/Prometheus",
                              f"{g}_Prometheus_1", 900.0, "Prometheus", "Prometheus", "Defense"])
        for j, (typ, sub) in enumerate([("RM", "RM_Type_I"), ("Gabija", "Gabija"), ("CBASS", "CBASS_I")]):
            if (k + j) % 2 == 0 and k != 7:
                sid = f"{g}_{typ}_1"
                sys_rows.append([g, sid, typ, sub, "Defense"])
                gene_rows.append([g, f"{g}_1", decoys[0][0], typ, 1, f"m/{typ}/{sub}", sid, 100.0,
                                  typ, sub, "Defense"])
        if k == 2:
            sys_rows.append([g, f"{g}_ADF_1", "Anti_CBASS", "Anti_CBASS", "Antidefense"])
        if k == 1:          # GTDB-style Prodigal proteome; the ProA runs off the contig start
            fasta(faa / f"{g}.faa", [
                (f"{n} # {1 + 300 * j} # {300 * (j + 1)} # 1 # ID=1_{j + 1};"
                 f"partial={'10' if n == prom_hits.get(g) else '00'};start_type=ATG;"
                 "rbs_motif=None;rbs_spacer=None;gc_cont=0.50", sq) for j, (n, sq) in enumerate(recs)])
        else:
            fasta(faa / f"{g}.faa", recs)
        crecs = []
        for n, s in recs:
            hdr = (f"lcl|NZ_CP4.1_cds_{n}_1 [locus_tag=L4] [protein_id={n}]"
                   if n.startswith("WP_") else n)
            crecs.append((hdr, back_translate(s, rng)))
        write(cds / f"{g}.fna", "".join(f">{h}\n{s}\n" for h, s in crecs))
        fasta(bi / "genomes" / f"{g}.fasta", [("contig1", "".join(rng.choice("ACGT") for _ in range(3000)))])
    tsv(bi / "df_systems.tsv", ["genome_id", "sys_id", "type", "subtype", "activity"], sys_rows)
    tsv(bi / "df_genes.tsv", ["genome_id", "replicon", "hit_id", "gene_name", "hit_pos", "model_fqn",
                              "sys_id", "hit_score", "type", "subtype", "activity"], gene_rows)
    tsv(bi / "taxonomy.tsv", ["genome_id", "classification"],
        [[g, f"d__Bacteria;p__{'Pseudomonadota' if k < 4 else 'Bacillota'};c__C{k % 2};o__;f__;"
             f"g__G{k};s__G{k} sp"] for k, g in enumerate(BG)])
    tsv(bi / "qc.tsv", ["genome_id", "completeness", "contamination", "genome_size"],
        [[g, 95 + k % 4, 7.5 if k == 6 else 1.0, 3000000 + 10000 * k] for k, g in enumerate(BG)])
    write(bi / "gtdb_bac120.tree",
          "((((RS_GCF_000010001.1:0.1,RS_GCF_000010002.1:0.1)'100.0:g__G':0.1,RS_GCF_000010003.1:0.2)"
          "'99.0:f__F; g__H':0.1,(RS_GCF_000010004.1:0.2,RS_GCF_000099991.1:0.2)'97.2':0.1)"
          "'100.0:p__Pseudomonadota':0.3,(((RS_GCF_000010005.1:0.1,RS_GCF_000010006.1:0.1)'100.0':0.1,"
          "RS_GCF_000010007.1:0.2)'98':0.1,(GB_GCA_000010008.1:0.3,RS_GCF_000099992.1:0.3)'90':0.1)"
          "'100.0:p__Bacillota':0.3)'100.0:d__Bacteria';\n")
    tsv(bi / "niche.tsv", ["genome_id", "oxygen", "host_associated"],
        [[g, ["anaerobe", "aerobe"][k % 2] if k != 6 else "",
          ["no", "yes"][k % 3 == 0] if k != 6 else ""] for k, g in enumerate(BG)])
    # ---- references, Pfam, DBs, fake viral tools -------------------------------- #
    fasta(work / "refs" / "prometheus_reference.faa", [(REF_SEQ + " ProA reference", prot[REF_SEQ])])
    fasta(work / "refs" / "paralog_reference.faa", [(n, prot[n]) for n in PARA_SEQS])
    pf = mk(work / "db" / "pfam")
    shutil.copy(FIXTURES / "pfam_proA_subset.hmm", pf / "Pfam-A.hmm")
    run(["hmmpress", "-f", str(pf / "Pfam-A.hmm")])
    write(work / "db" / "genomad" / "version.txt", "1.7\n")          # genomad download-database layout
    mk(work / "db" / "checkv" / "genome_db")                          # checkv download_database layout
    mk(work / "db" / "checkv" / "hmm_db")
    mk(work / "db" / "checkv_incomplete" / "genome_db")               # no hmm_db
    fb = mk(work / "fakebin")
    for name, body in (("genomad", FAKE_GENOMAD), ("checkv", FAKE_CHECKV)):
        write(fb / name, body.replace("{py}", PYTHON)).chmod(0o755)
    cfgd = {
        "samples": "arch/species_metadata.tsv", "metadata": "arch/species_metadata.tsv",
        "results_dir": "results", "logs_dir": "logs",
        "db": {"pfam_hmm": "db/pfam/Pfam-A.hmm", "genomad": "db/genomad", "checkv": "db/checkv"},
        "refs": {"prometheus_faa": "refs/prometheus_reference.faa",
                 "paralog_faa": "refs/paralog_reference.faa"},
        "prometheus": {"min_len_aa": 500, "validation_tree_model": "LG+G4"},   # MFP: same rule, slower
        "evolution": {"gene_tree_model": "LG+G4", "selection": {"run": False}},
        "tree": {"ufboot": 1000, "alrt": 1000},
        "threads": {"small": 1, "medium": 2, "large": 2},
        "plotting": {"formats": ["png", "svg"]},
        "bacteria": {
            "enabled": True, "df_systems": "bact/df_systems.tsv", "df_genes": "bact/df_genes.tsv",
            "faa_dir": "bact/faa", "cds_dir": "bact/cds", "taxonomy": "bact/taxonomy.tsv",
            "qc": "bact/qc.tsv", "gtdb_bac_tree": "bact/gtdb_bac120.tree",
            "niche_metadata": "" if empty else "bact/niche.tsv",
            "genomes_dir": "" if empty else "bact/genomes", "genome_extension": "fasta",
            "crossdomain": {"n_subsample": 50, "saturation_pairs": 20}},
    }
    write(work / "config" / "test.yaml", json.dumps(cfgd, indent=1))   # JSON is valid YAML
    cfgd["db"]["genomad"] = "db/missing_genomad"
    write(work / "config" / "test_nodb.yaml", json.dumps(cfgd, indent=1))
    cfgd["db"].update(genomad="db/genomad", checkv="db/checkv_incomplete")
    write(work / "config" / "test_badcheckv.yaml", json.dumps(cfgd, indent=1))
    return prom_hits


def smk(work, targets, *extra, prefix_bin=None, check=True, configfile="config/test.yaml",
        timeout=1200):
    path = [str(work / "fakebin"), str(Path(PYTHON).parent)]
    if prefix_bin:
        path.insert(0, str(prefix_bin))
    if TOOLS_BIN:
        path.append(TOOLS_BIN)
    path += [str(Path(RSCRIPT).parent), "/usr/bin", "/bin", "/usr/sbin", "/sbin"]
    e = dict(os.environ, PATH=os.pathsep.join(path))
    cmd = [find_snakemake(), *targets, "-s", "workflow/Snakefile.cdtest", "--configfile",
           configfile, "--cores", "2", "--nolock", "--rerun-triggers", "mtime",
           "--show-failed-logs", *extra]
    p = subprocess.run(cmd, cwd=work, env=e, capture_output=True, text=True, timeout=timeout)
    if check and p.returncode != 0:
        raise AssertionError(f"snakemake failed ({p.returncode}):\n{p.stdout[-4000:]}\n{p.stderr[-6000:]}")
    return p


E2E = {}


def e2e_work():
    """Build + run the full scenario once; later tests reuse it (or re-raise its error)."""
    if "error" in E2E:
        raise AssertionError(f"end-to-end run failed earlier: {E2E['error'][:1500]}")
    if "work" not in E2E:
        work = mk(ROOT / "e2e")
        try:
            E2E["hits"] = build_e2e(work)
            E2E["proc"] = smk(work, ["cd_e2e"])
        except Exception as exc:  # noqa: BLE001
            E2E["error"] = str(exc)
            raise
        E2E["work"] = work
    return E2E["work"]


@S.test
def e2e_bacterial_branch():
    w = e2e_work()
    r = w / "results"
    samples = read_tsv(r / "bacteria" / "bacterial_samples.tsv")
    assert list(samples.genome_id) == BG and "completeness" in samples.columns
    assert (r / "bacteria" / "bacterial_representatives.txt").read_text().split() == BG
    cand = dict(read_fa(r / "bacteria" / "prometheus" / "candidates.faa"))
    assert set(cand) == {f"{g}__{h}" for g, h in E2E["hits"].items()}
    assert set(dict(read_fa(r / "bacteria" / "prometheus" / "bacterial_prometheus_cds.ffn"))) == set(cand)
    dom = (r / "bacteria" / "prometheus" / "domains.domtbl").read_text()
    assert "DUF4011" in dom and "PF13195" in dom                      # real hmmscan --cut_ga
    for g in BG:
        assert (r / "bacteria" / "defense" / g / f"{g}_defense_finder_systems.tsv").exists()
    sysg3 = read_tsv(r / "bacteria" / "defense" / BG[2] / f"{BG[2]}_defense_finder_systems.tsv")
    assert "Antidefense" not in set(sysg3["activity"])
    vt = r / "bacteria" / "prometheus" / "validation_tree"
    assert "sp_P51530_DNA2_HUMAN\tParalog_ref" in (vt / "reference_labels.tsv").read_text()
    vtree = (vt / "validation.treefile").read_text()
    assert "sp_P51530_DNA2_HUMAN" in vtree                              # ids survive IQ-TREE
    # regression: SH-aLRT/UFBoot labels exist (a -fast tree has none, so every bacterial
    # placement would be 'unsupported' and tier A unreachable in bacteria only)
    assert re.search(r"\)[0-9.]+/[0-9.]+:", vtree), vtree[:300]
    cat = read_tsv(r / "bacteria" / "prometheus" / "bacterial_prometheus_catalog.tsv")
    assert set(cat.locus_id) == set(cand)
    # Prodigal partial flags reach validate_prometheus (same evidence as archaea)
    gt = read_tsv(r / "bacteria" / "genes" / BG[1] / f"{BG[1]}.genes.tsv").set_index("protein_id")
    assert len(gt) == 3 and gt.loc[f"{BG[1]}_1_3", "partial"] == "10"
    assert not (r / "bacteria" / "genes" / BG[0]).exists()               # plain headers
    xr = read_tsv(r / "bacteria" / "prometheus" / "extraction_report.tsv").set_index("genome_id")
    assert xr.loc[BG[1], "genes_table"] == "prodigal" and xr.loc[BG[0], "genes_table"] == "none"
    c = cat.set_index("locus_id")
    assert c.loc[f"{BG[1]}__{BG[1]}_1_3", "fragment_evidence"] == "prodigal_partial=10"
    assert c.loc[f"{BG[1]}__{BG[1]}_1_3", "fragmented"] == "True"
    assert c.loc[f"{BG[0]}__{BG[0]}_1_3", "fragment_evidence"].startswith("length_proxy")
    assert "== domain" in (r / "bacteria" / "prometheus" / "domains.hmmscan.txt").read_text()
    m = read_tsv(r / "bacteria" / "bacterial_master_table.tsv").set_index("genome_id")
    assert list(m.index) == BG and m.loc[BG[6], "pass_primary"] == "False"   # contamination 7.5
    tips = read_tsv(r / "bacteria" / "scaffolds" / "bacterial_species_tree_tips.tsv").set_index("genome_id")
    assert tips.loc[BG[7], "match"] == "paired_accession" and tips.loc[BG[0], "match"] == "prefix_stripped"
    sp = (r / "bacteria" / "scaffolds" / "bacterial_species_tree.treefile").read_text()
    assert tips_of(sp) == sorted(BG)
    for f in ("distribution/bacterial_prevalence_by_taxon.tsv",
              "distribution/bacterial_defense_matrix.tsv", "distribution/bacterial_cooccurrence.tsv",
              "diversity/bacterial_architecture_typology.tsv"):
        assert (r / "bacteria" / f).exists(), f
    dm = read_tsv(r / "bacteria" / "distribution" / "bacterial_defense_matrix.tsv")
    assert "Anti_CBASS" not in dm.columns and "Prometheus" in dm.columns


@S.test
def e2e_crossdomain_tree_and_analysis():
    w = e2e_work()
    r = w / "results" / "crossdomain"
    tree = (r / "tree" / "crossdomain.treefile").read_text().strip()
    assert tree != "();" and tree.count(",") == 9, tree[:200]       # regression: real 10-tip tree
    assert len((r / "tree" / "crossdomain.ufboot").read_text().splitlines()) == 1000
    s = summary(r / "crossdomain_summary.tsv")
    assert s["status"] == "complete" and s["n_archaeal_tips"] == "2"
    assert s["n_bacterial_tips"] == "5" and s["n_paralog_tips"] == "3"    # 4 validated + 1 ref
    assert s["n_bootstrap_trees"] == "1000" and s["direction_method"] == "topology-based"
    assert s["inferred_direction"] in ("bacteria->archaea", "archaea->bacteria", "ambiguous")
    tips = read_tsv(r / "crossdomain_tip_domains.tsv")
    assert "Unknown" not in set(tips.domain) and len(tips) == 10
    import crossdomain_analysis as ca                                 # our reader == dendropy
    import dendropy
    t = ca.UTree(tree)
    dt = dendropy.Tree.get(data=tree, schema="newick", preserve_underscores=True)
    assert sorted(t.by_label) == sorted(lf.taxon.label for lf in dt.leaf_node_iter())
    n_lab = sum(1 for nd in dt.internal_nodes() if nd.label)
    assert len(t.support) == n_lab and all(a is not None and u is not None for a, u in t.support.values())
    sat = read_tsv(r / "selection" / "crossdomain_saturation.tsv")
    assert set(sat.pair_class) == {"within_archaea", "within_bacteria", "cross_domain"}
    for f in ("crossdomain_tree.png", "crossdomain_tree.svg", "crossdomain_comparison.png",
              "crossdomain_comparison.svg", "viral_rloop_crossdomain.png"):
        assert (r / "figures" / f).stat().st_size > 1000, f


@S.test
def e2e_codon_alignment_keeps_both_domains():
    w = e2e_work()
    r = w / "results" / "crossdomain" / "selection"
    codon = dict(read_fa(r / "crossdomain.codon.fasta"))
    assert len([n for n in codon if n.startswith("ARC")]) == 2, sorted(codon)   # archaea kept
    assert len([n for n in codon if n.startswith("GCF_")]) == 4, sorted(codon)
    lens = {len(s) for s in codon.values()}
    assert len(lens) == 1 and lens.pop() % 3 == 0
    assert tips_of((r / "crossdomain_selection.nwk").read_text()) == sorted(codon)
    cdsids = set(dict(read_fa(r / "crossdomain_all_cds.ffn")))
    assert cdsids == set(codon) and "ARC3__ARC3_1_1" not in cdsids     # restricted to validated ids
    summ = read_tsv(r / "crossdomain_selection_summary.tsv")
    assert list(summ.columns)[:3] == ["method", "statistic", "value"]
    rep = read_tsv(r / "crossdomain_codon_report.tsv")
    assert set(rep[rep.in_proteins == "1"].status) == {"kept"}           # every expected locus
    assert set(rep[rep.in_proteins == "0"].status) == {"no_cds"}         # refs/paralogs only


@S.test
def e2e_no_bacterial_cds_selection_not_applicable():
    """bacteria.cds_dir unset (config default): no cross-domain codon alignment, and the
    selection summary says why -- instead of codon_align.py aborting on 'no_cds' loci or
    an archaea-only test being reported as cross-domain."""
    w = e2e_work()
    work = ROOT / "e2e_nocds"
    shutil.copytree(w, work, symlinks=True)
    cfgd = json.loads((work / "config" / "test.yaml").read_text())
    cfgd["bacteria"]["cds_dir"] = ""
    cfgd["evolution"]["selection"]["run"] = True
    write(work / "config" / "test_nocds.yaml", json.dumps(cfgd, indent=1))
    target = "results/crossdomain/selection/crossdomain_selection_summary.tsv"
    smk(work, [target], "--forcerun", "crossdomain_codon_align", configfile="config/test_nocds.yaml")
    r = work / "results" / "crossdomain" / "selection"
    assert (r / "crossdomain.codon.fasta").stat().st_size == 0
    assert (r / "crossdomain_selection.nwk").read_text().strip() == "();"
    rep = read_tsv(r / "crossdomain_codon_report.tsv")
    assert len(rep) == 0 and list(rep.columns) == ["id", "status", "in_proteins", "n_aa",
                                                   "n_codons", "note"]
    summ = read_tsv(r / "crossdomain_selection_summary.tsv")
    assert len(summ) >= 4 and set(summ.note) == {
        "not run: bacteria.cds_dir not set (cross-domain dN/dS needs bacterial CDS)"}
    assert "bacteria.cds_dir is not set" in (work / "logs" / "crossdomain" / "codon_align.log").read_text()


@S.test
def e2e_hyphy_accepts_selection_inputs():
    w = e2e_work()
    if not have("hyphy"):
        raise AssertionError("hyphy not on PATH (PROM_TOOLS_BIN)")
    r = w / "results" / "crossdomain" / "selection"
    out = ROOT / "fel.json"
    run(["hyphy", "fel", "--alignment", str(r / "crossdomain.codon.fasta"), "--tree",
         str(r / "crossdomain_selection.nwk"), "--output", str(out), "CPU=2"], timeout=600)
    js = json.loads(out.read_text())
    assert "MLE" in js and js["input"]["number of sequences"] == 6


@S.test
def e2e_niche_and_viral_modules():
    w = e2e_work()
    r = w / "results" / "bacteria" / "ecology"
    nm = read_tsv(r / "bacterial_niche_master.tsv")
    assert len(nm) == 7 and BG[6] not in set(nm.genome_id)            # unlabelled genome dropped
    assert {"niche_oxygen_aerobe", "niche_oxygen_anaerobe", "niche_host_associated"} <= set(nm.columns)
    assoc = read_tsv(r / "bacterial_niche_associations.tsv")
    assert list(assoc.columns)[:9] == ["predictor", "level", "model", "estimate", "std_error",
                                       "p_value", "p_adj", "n", "note"]   # SCHEMA (+ additive)
    assert {"niche_oxygen_aerobe", "niche_oxygen_anaerobe"} <= set(assoc.predictor)
    for f in ("bacterial_niche_associations_sensitivity.tsv",      # A+B re-fit, as archaea
              "bacterial_viral_associations_sensitivity.tsv"):
        assert list(read_tsv(r / f).columns)[:9] == list(assoc.columns)[:9], f
    nlog = (w / "logs" / "bacteria" / "niche.log").read_text()
    assert "== primary tiers A" in nlog and "== sensitivity tiers A,B" in nlog
    assert (r / "viral" / BG[0] / "rloop_rlfs_calls.tsv").exists()
    vs = read_tsv(r / "viral" / "bacterial_viral_summary.tsv").set_index("genome_id")
    assert list(vs.index) == BG
    assert vs.loc[BG[0], "n_proviruses"] == "1" and vs.loc[BG[1], "n_proviruses"] == "0"
    q = read_tsv(r / "viral" / BG[1] / "checkv" / "quality_summary.tsv")
    assert len(q) == 0 and len(q.columns) == 14                        # CheckV header, empty input
    assert (r / "bacterial_viral_associations.tsv").exists()
    logs = w / "logs"
    for rel in ("bacteria/extract.log", "bacteria/domains.log", "crossdomain/tree.log",
                "crossdomain/analysis.log", f"bacteria/genomad/{BG[0]}.log", "bacteria/niche.log",
                "bacteria/samples.log"):
        assert (logs / rel).exists(), rel


@S.test
def e2e_all_empty_scenario():
    work = mk(ROOT / "e2e_empty")
    build_e2e(work, empty=True)
    smk(work, ["cd_e2e"])
    r = work / "results"
    assert (r / "bacteria" / "prometheus" / "candidates.faa").read_text() == ""
    assert (r / "bacteria" / "prometheus" / "domains.domtbl").read_text().startswith("#")
    vt = r / "bacteria" / "prometheus" / "validation_tree" / "validation.treefile"
    assert vt.read_text().strip() == "();"
    assert len(read_tsv(r / "bacteria" / "prometheus" / "bacterial_prometheus_catalog.tsv")) == 0
    assert (r / "crossdomain" / "tree" / "crossdomain.treefile").read_text().strip() == "();"
    assert summary(r / "crossdomain" / "crossdomain_summary.tsv")["status"] == "incomplete"
    assert (r / "crossdomain" / "selection" / "crossdomain_selection.nwk").read_text().strip() == "();"
    assert (r / "crossdomain" / "figures" / "crossdomain_tree.png").exists()
    assert set(read_tsv(r / "bacteria" / "bacterial_master_table.tsv").prometheus_present) == {"False"}


@S.test
def e2e_crashed_tools_fail_their_rules():
    w = e2e_work()
    fb = mk(ROOT / "failbin")
    for tool in ("iqtree2", "hmmscan"):
        write(fb / tool, FAIL_TOOL).chmod(0o755)
    tree = "results/crossdomain/tree/crossdomain.treefile"
    p = smk(w, [tree], "--forcerun", "crossdomain_tree", prefix_bin=fb, check=False)
    assert p.returncode != 0 and "crossdomain_tree" in (p.stdout + p.stderr)
    assert not (w / tree).exists()                                     # no '();' placeholder
    dom = "results/bacteria/prometheus/domains.domtbl"
    p = smk(w, [dom], "--forcerun", "bact_domain_scan", prefix_bin=fb, check=False)
    assert p.returncode != 0 and not (w / dom).exists()
    gm = f"results/bacteria/ecology/viral/{BG[0]}/{BG[0]}_summary/{BG[0]}_virus.fna"
    p = smk(w, [gm], "--forcerun", "bact_genomad", configfile="config/test_nodb.yaml", check=False)
    assert p.returncode != 0
    assert "geNomad database not found" in (w / "logs" / "bacteria" / "genomad" / f"{BG[0]}.log").read_text()
    # CheckV: exit 0 without a table, a non-1.0.3 header, an incomplete database
    q = f"results/bacteria/ecology/viral/{BG[0]}/checkv/quality_summary.tsv"
    qlog = w / "logs" / "bacteria" / "checkv" / f"{BG[0]}.log"
    for name, body in (("silent", CHECKV_SILENT), ("badhdr", CHECKV_BADHDR)):
        fbc = mk(ROOT / f"checkv_{name}")
        write(fbc / "checkv", body).chmod(0o755)
        p = smk(w, [q], "--forcerun", "bact_checkv", prefix_bin=fbc, check=False)
        assert p.returncode != 0 and not (w / q).exists(), name
        assert "CheckV exited 0 but wrote no quality_summary.tsv" in qlog.read_text(), name
    p = smk(w, [q], "--forcerun", "bact_checkv", configfile="config/test_badcheckv.yaml", check=False)
    assert p.returncode != 0 and "CheckV database not found at 'db/checkv_incomplete'" in qlog.read_text()


if __name__ == "__main__":
    try:
        S.main("CROSSDOMAIN TESTS PASSED")
    finally:
        if not os.environ.get("CROSSDOMAIN_KEEP"):
            shutil.rmtree(ROOT, ignore_errors=True)
