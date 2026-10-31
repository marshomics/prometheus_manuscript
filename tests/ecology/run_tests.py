"""Stage 6 (ecology) synthetic test suite.   bash tests/run_leaf.sh ecology

Covers rloop_predict.py (call-for-call against the reference QmRLFS-finder v2.0 output on
real virus genomes + hand-derived constructs), aggregate_viral.py (provirus vs viral-contig
accounting, NA semantics, strict failure), ecology_stats.R (real phylolm::phyloglm runs
with known answers, the Contract v2 item 6 presence rule, bacterial call shapes,
failures), both plotting scripts, and the ecology.smk rules themselves (a Snakemake
mini-workflow with stub geNomad/CheckV binaries that mimic the verified 1.8.0 / 1.0.3
output formats and failure modes). Prints ECOLOGY TESTS PASSED only when all pass.
"""
import atexit
import math
import os
import py_compile
import random
import shutil
import stat
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, rscript, run, read_tsv, REPO,  # noqa: E402
                  SCRIPTS_PY, SCRIPTS_R, PYTHON, RSCRIPT, env)

FIX = HERE / "fixtures"
suite = Suite("ecology")
OWNED_PY = ["rloop_predict.py", "aggregate_viral.py", "plot_ecology.py", "plot_ecology_extra.py"]
SMK = os.environ.get("PROM_SNAKEMAKE") or shutil.which("snakemake") or \
    "/Users/jmarsh/.claude-science/conda/envs/prom-smk/bin/snakemake"

GENOMAD_COLS = ["seq_name", "length", "topology", "coordinates", "n_genes", "genetic_code",
                "virus_score", "fdr", "n_hallmarks", "marker_enrichment", "taxonomy"]
CHECKV_COLS = ["contig_id", "contig_length", "provirus", "proviral_length", "gene_count",
               "viral_genes", "host_genes", "checkv_quality", "miuvig_quality", "completeness",
               "completeness_method", "contamination", "kmer_freq", "warnings"]


def _scratch_root():
    """Absolute scratch root. TMPDIR may be relative (this sandbox sets './.tmp', which
    os.path.abspath leaves relative), and the Snakemake/R subprocesses run with other
    working directories, so every test path must be absolute. Removed at exit unless
    ECO_KEEP_TMP=1."""
    for base in (os.environ.get("ECO_TEST_TMP"), "/tmp", tempfile.gettempdir()):
        if base and os.path.isabs(base) and os.path.isdir(base) and os.access(base, os.W_OK):
            root = Path(tempfile.mkdtemp(prefix="ecology_tests_", dir=base))
            if os.environ.get("ECO_KEEP_TMP") != "1":
                atexit.register(shutil.rmtree, root, True)
            return root
    raise RuntimeError("no absolute writable scratch directory (set ECO_TEST_TMP)")


SCRATCH = _scratch_root()


def tmpdir(prefix="t_"):
    return Path(tempfile.mkdtemp(prefix=prefix, dir=SCRATCH))


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def tsv(path, header, rows):
    lines = ["\t".join(header)] + ["\t".join(str(v) for v in r) for r in rows]
    return write(path, "\n".join(lines) + "\n")


def fnum(x):
    return float("nan") if x in ("", "NA", None) else float(x)


def plot_env():
    e = env()
    e["PROM_DPI"] = "60"
    e["MPLCONFIGDIR"] = str(tmpdir("mpl_"))
    return e


def run_py(script, *args, check=True, extra_env=None):
    import subprocess
    e = plot_env()
    if extra_env:
        e.update(extra_env)
    p = subprocess.run([PYTHON, str(SCRIPTS_PY / script), *map(str, args)], cwd=REPO, env=e,
                       text=True, capture_output=True, timeout=600)
    if check and p.returncode != 0:
        raise AssertionError(f"{script} failed ({p.returncode}):\n{p.stdout[-2000:]}\n{p.stderr[-3000:]}")
    return p


# =============================================================================== static
@suite.test
def test_compile_python():
    for s in OWNED_PY:
        py_compile.compile(str(SCRIPTS_PY / s), doraise=True)
    for s in HERE.glob("*.py"):
        py_compile.compile(str(s), doraise=True)


@suite.test
def test_parse_r():
    p = run([RSCRIPT, "-e", f'invisible(parse("{SCRIPTS_R / "ecology_stats.R"}")); cat("PARSE_OK")'])
    assert "PARSE_OK" in p.stdout


@suite.test
def test_rules_hygiene():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"), str(REPO / "workflow/rules/ecology.smk")])
    assert "RULES HYGIENE OK" in p.stdout, p.stdout


# =============================================================================== R-loops
REF_COLS = ["model", "start_RIZ", "end_RIZ", "length_RIZ", "G_RIZ", "3Gs_RIZ", "4Gs_RIZ",
            "perc_G_RIZ", "Linker", "start_REZ", "end_REZ", "length_REZ", "G_REZ", "3Gs_REZ",
            "4Gs_REZ", "perc_G_REZ", "strand"]


def _rloop(fasta, *extra):
    d = tmpdir("rloop_")
    out, calls = d / "rloop.tsv", d / "calls.tsv"
    py("rloop_predict.py", "--fasta", fasta, "--genome", "G", "--calls", calls, "--out", out, *extra)
    return read_tsv(out), read_tsv(calls)


def _key(row, name_col):
    return tuple([row[name_col].split()[0]] + [str(row[c]) for c in REF_COLS] +
                 [row["location"].split(":")[-1]])


@suite.test
def test_rloop_matches_reference_qmrlfs():
    """Call-for-call identity with the reference QmRLFS-finder v2.0 (fixtures/README.md)."""
    for fasta, ref in (("rloop_real_viruses.fna", "qmrlfs_ref_real_viruses.table.tsv"),
                       ("rloop_reference_input_v1.fna", "qmrlfs_ref_synthetic_v1.table.tsv")):
        per_seq, calls = _rloop(FIX / fasta)
        refdf = read_tsv(FIX / ref)
        mine = [_key(r, "provirus_id") for _, r in calls.iterrows()]
        theirs = [_key(r, "#Name") for _, r in refdf.iterrows()]
        assert mine == theirs, (f"{fasta}: {len(mine)} calls vs reference {len(theirs)}; "
                                f"first diff {next((a, b) for a, b in zip(mine, theirs) if a != b)}")
        assert per_seq["n_rlfs"].astype(int).sum() == len(theirs)
    per_seq, _ = _rloop(FIX / "rloop_real_viruses.fna")
    got = dict(zip(per_seq["provirus_id"], per_seq["n_rlfs"].astype(int)))
    # lambda: one RLFS (van den Berg et al. 2024: lambda carries one predicted R-loop)
    assert got == {"NC_001416.1": 1, "NC_001338.1": 0, "NC_007217.1": 25}, got
    prone = dict(zip(per_seq["provirus_id"], per_seq["rloop_prone"]))
    assert prone == {"NC_001416.1": "True", "NC_001338.1": "False", "NC_007217.1": "True"}
    loci = dict(zip(per_seq["provirus_id"], per_seq["n_rlfs_loci"].astype(int)))
    assert loci["NC_007217.1"] == 13 and loci["NC_001416.1"] == 1


@suite.test
def test_rloop_constructs_known_answers():
    per_seq, calls = _rloop(FIX / "rloop_synthetic.fna")
    s = per_seq.set_index("provirus_id")
    want = {  # n_rlfs, plus, minus, m1, m2
        "syn_plus_m1": (1, 1, 0, 1, 0), "syn_minus_m2": (1, 0, 1, 0, 1),
        "syn_riz_no_rez": (0, 0, 0, 0, 0), "syn_rez_at_end": (0, 0, 0, 0, 0),
        "syn_rez_at_end_plus1": (1, 1, 0, 1, 0), "syn_low_g_riz": (0, 0, 0, 0, 0),
        "syn_n_break": (0, 0, 0, 0, 0), "syn_lowercase": (1, 1, 0, 1, 0),
        "syn_random_gc50_1": (0, 0, 0, 0, 0), "syn_random_gc50_2": (0, 0, 0, 0, 0),
        "syn_random_gc50_3": (0, 0, 0, 0, 0)}
    for name, (n, p, m, m1, m2) in want.items():
        got = tuple(int(s.loc[name, c]) for c in ("n_rlfs", "n_rlfs_plus", "n_rlfs_minus",
                                                   "n_rlfs_m1", "n_rlfs_m2"))
        assert got == (n, p, m, m1, m2), (name, got)
        assert s.loc[name, "rloop_prone"] == str(n > 0)
    c = calls.set_index("provirus_id")
    # RIZ GGGTGGGTGGG at [300,311) (9/11 G = 81.82 %), linker 0, REZ = the LONGEST >=40 %-G
    # window: 'T' + (GA)x100 + 49 A/T = 250 nt with exactly 100 G (40.0 %), ending at 561
    r = c.loc["syn_plus_m1"]
    assert (r["start_RIZ"], r["end_RIZ"], r["perc_G_RIZ"], r["Linker"], r["start_REZ"],
            r["end_REZ"], r["length_REZ"], r["G_REZ"], r["perc_G_REZ"]) == \
        ("300", "311", "81.82", "0", "311", "561", "250", "100", "40.0"), dict(r)
    # minus strand: construct = reverse complement (length 810): RIZ [501,510), REZ [251,501)
    r = c.loc["syn_minus_m2"]
    assert (r["strand"], r["start_RIZ"], r["end_RIZ"], r["start_REZ"], r["end_REZ"],
            r["location"]) == ("-", "501", "510", "251", "501", "syn_minus_m2:252-510"), dict(r)
    # REZ must end strictly before the sequence end: 100 nt left -> none, 101 -> [311,411)
    r = c.loc["syn_rez_at_end_plus1"]
    assert (r["start_REZ"], r["end_REZ"], r["length_REZ"]) == ("311", "411", "100")
    assert dict(c.loc["syn_lowercase"])["end_REZ"] == "561"
    # GC of the 812-nt construct: G = 9 + 100, C = 0 -> 109/812
    assert abs(float(s.loc["syn_plus_m1", "gc"]) - 109 / 812) < 1e-5
    # RIZ-only (permissive, non-QmRLFS) mode counts the REZ-less RIZ
    per_seq2, _ = _rloop(FIX / "rloop_synthetic.fna", "--require-rez", "False")
    s2 = per_seq2.set_index("provirus_id")
    assert int(s2.loc["syn_riz_no_rez", "n_rlfs"]) == 1
    assert int(s2.loc["syn_low_g_riz", "n_rlfs"]) == 0       # %G filter still applies
    # min-rlfs threshold acts on merged loci (m1/m2 duplicates of one locus count once)
    per_seq3, _ = _rloop(FIX / "rloop_synthetic.fna", "--min-rlfs", "2")
    s3 = per_seq3.set_index("provirus_id")
    assert s3.loc["syn_random_gskew", "n_rlfs_loci"] == "1"
    assert s3.loc["syn_random_gskew", "rloop_prone"] == "False"


@suite.test
def test_rloop_crlf_and_wrapping():
    """CRLF line ends and different line wrapping must not change any call."""
    d = tmpdir("rloop_crlf_")
    text = (FIX / "rloop_real_viruses.fna").read_text()
    seqs, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line; seqs[name] = []
        else:
            seqs[name].append(line.strip())
    crlf = "".join(f"{h}\r\n" + "".join(f"{''.join(v)[i:i + 37]}\r\n" for i in range(0, len(''.join(v)), 37))
                   for h, v in seqs.items())
    write(d / "crlf.fna", crlf)
    a, ca = _rloop(FIX / "rloop_real_viruses.fna")
    b, cb = _rloop(d / "crlf.fna")
    assert a.drop(columns=[]).equals(b), (a, b)
    assert ca.equals(cb)


@suite.test
def test_rloop_empty_and_failures():
    d = tmpdir("rloop_e_")
    write(d / "empty.fna", "")
    py("rloop_predict.py", "--fasta", d / "empty.fna", "--genome", "G", "--out", d / "o.tsv",
       "--calls", d / "c.tsv")
    assert (d / "o.tsv").read_text().splitlines() == [
        "genome_id\tprovirus_id\tlength\tn_rlfs\trloop_prone\tn_rlfs_loci\tn_rlfs_plus\t"
        "n_rlfs_minus\tn_rlfs_m1\tn_rlfs_m2\tgc\trlfs_loci_per_kb"]
    assert len((d / "c.tsv").read_text().splitlines()) == 1
    write(d / "dup.fna", ">a\nGGGTGGGTGGG\n>a\nACGT\n")
    write(d / "junk.fna", "this is not fasta\n")
    for args in (["--fasta", d / "dup.fna"], ["--fasta", d / "missing.fna"],
                 ["--fasta", d / "junk.fna"],
                 ["--fasta", d / "empty.fna", "--require-rez", "maybe"],
                 ["--fasta", d / "empty.fna", "--models", "m3"]):
        p = py("rloop_predict.py", *args, "--genome", "G", "--out", d / "x.tsv", check=False)
        assert p.returncode != 0, args


# =============================================================================== aggregate_viral
def _viral_genome(root, gid, vs_rows, cv_rows, rl_rows):
    tsv(root / gid / f"{gid}_summary" / f"{gid}_virus_summary.tsv", GENOMAD_COLS, vs_rows)
    tsv(root / gid / "checkv" / "quality_summary.tsv", CHECKV_COLS, cv_rows)
    tsv(root / gid / "rloop.tsv", ["genome_id", "provirus_id", "length", "n_rlfs", "rloop_prone", "gc"],
        rl_rows)


def _gnm(name, length, topo):
    coord = name.split("provirus_")[1].replace("_", "-") if topo == "Provirus" else "NA"
    return [name, length, topo, coord, 30, 11, 0.95, "NA", 2, 3.1, "Viruses;Duplodnaviria"]


def _cv(name, length, quality, comp):
    return [name, length, "No", "NA", 30, 10, 1, quality, "Genome-fragment", comp, "AAI-based", 0.0, 1.0, ""]


def _viral_fixture():
    root = tmpdir("viral_")
    g1 = "GCA_1"
    p1, p2 = "c1|provirus_101_30100", "c2|provirus_5_20004"
    _viral_genome(root, g1,
                  [_gnm(p1, 30000, "Provirus"), _gnm(p2, 20000, "Provirus"),
                   _gnm("c7", 40000, "DTR"), _gnm("c9", 5000, "No terminal repeats")],
                  [_cv(p1, 30000, "Medium-quality", 80.5), _cv(p2, 20000, "Not-determined", "NA"),
                   _cv("c7", 40000, "Complete", 100.0), _cv("c9", 5000, "Low-quality", 12.0)],
                  [[g1, p1, 30000, 2, "True", 0.42], [g1, p2, 20000, 0, "False", 0.38],
                   [g1, "c7", 40000, 1, "True", 0.5], [g1, "c9", 5000, 0, "False", 0.3]])
    _viral_genome(root, "GCA_2", [], [], [])                 # geNomad ran, found nothing
    _viral_genome(root, "000123", [_gnm("x|provirus_1_9000", 9000, "Provirus")],
                  [_cv("x|provirus_1_9000", 9000, "High-quality", 95.0)],
                  [["000123", "x|provirus_1_9000", 9000, 0, "False", 0.61]])
    samples = tsv(root / "samples.tsv", ["genome_id", "completeness"],
                  [[g1, 95], ["GCA_2", 99], ["000123", 91]])
    return root, samples


@suite.test
def test_aggregate_viral_known_answers():
    root, samples = _viral_fixture()
    out = root / "viral_summary.tsv"
    py("aggregate_viral.py", "--viral-dir", root, "--samples", samples, "--out", out)
    df = read_tsv(out).set_index("genome_id")
    assert list(df.reset_index().columns[:6]) == ["genome_id", "n_proviruses", "total_provirus_kb",
                                                  "mean_checkv_completeness", "n_rloop_proviruses",
                                                  "frac_rloop_proviruses"]
    assert list(df.index) == ["GCA_1", "GCA_2", "000123"]       # order + no int coercion
    g = df.loc["GCA_1"]
    assert g["n_proviruses"] == "2" and fnum(g["total_provirus_kb"]) == 50.0
    assert fnum(g["mean_checkv_completeness"]) == 80.5            # NA completeness ignored
    assert g["n_rloop_proviruses"] == "1" and fnum(g["frac_rloop_proviruses"]) == 0.5
    assert g["n_viral_contigs"] == "2" and fnum(g["total_viral_contig_kb"]) == 45.0
    assert g["n_viral_sequences"] == "4" and g["n_rloop_viral_sequences"] == "2"
    assert fnum(g["frac_rloop_viral_sequences"]) == 0.5
    assert g["n_proviruses_checkv_mq"] == "1"                    # Medium-quality only
    assert fnum(g["mean_provirus_kb"]) == 25.0
    assert abs(fnum(g["mean_provirus_gc"]) - 0.40) < 1e-9
    e = df.loc["GCA_2"]                                          # zeros, fractions undefined
    assert (e["n_proviruses"], e["n_viral_contigs"], e["n_rloop_proviruses"]) == ("0", "0", "0")
    assert e["frac_rloop_proviruses"] == "" and e["mean_checkv_completeness"] == ""
    assert e["frac_rloop_viral_sequences"] == ""
    assert df.loc["000123", "n_proviruses"] == "1"
    assert fnum(df.loc["000123", "frac_rloop_proviruses"]) == 0.0
    # bacterial samples table keyed by user_genome (common.smk fallback)
    s2 = tsv(root / "qc.tsv", ["user_genome", "completeness"], [["GCA_2", 99]])
    py("aggregate_viral.py", "--viral-dir", root, "--samples", s2, "--out", root / "b.tsv")
    assert read_tsv(root / "b.tsv")["genome_id"].tolist() == ["GCA_2"]
    # empty samples -> header-only
    s3 = tsv(root / "none.tsv", ["genome_id"], [])
    py("aggregate_viral.py", "--viral-dir", root, "--samples", s3, "--out", root / "e.tsv")
    assert (root / "e.tsv").read_text().splitlines()[0].startswith("genome_id\tn_proviruses")
    assert len((root / "e.tsv").read_text().splitlines()) == 1


@suite.test
def test_aggregate_viral_failures():
    def fails(mutate, needle):
        root, samples = _viral_fixture()
        mutate(root)
        p = py("aggregate_viral.py", "--viral-dir", root, "--samples", samples, "--out",
               root / "o.tsv", check=False)
        assert p.returncode != 0, needle
        assert needle in p.stderr, (needle, p.stderr[-500:])
        assert not (root / "o.tsv").exists()
    fails(lambda r: (r / "GCA_1" / "rloop.tsv").unlink(), "R-loop table missing")
    fails(lambda r: write(r / "GCA_2" / "checkv" / "quality_summary.tsv", ""), "empty file")
    fails(lambda r: tsv(r / "GCA_1" / "rloop.tsv", ["genome_id", "provirus_id", "length", "n_rlfs", "rloop_prone"],
                        [["GCA_1", "c1|provirus_101_30100", 30000, 1, "True"]]), "disagree")
    fails(lambda r: tsv(r / "000123" / "rloop.tsv", ["genome_id", "provirus_id", "length", "n_rlfs", "rloop_prone"],
                        [["000123", "x|provirus_1_9000", 8000, 0, "False"]]), "lengths differ")
    fails(lambda r: tsv(r / "000123" / "000123_summary" / "000123_virus_summary.tsv", GENOMAD_COLS,
                        [["x|provirus_1_9000", 9000, "Prophage", "1-9000", 1, 11, 0.9, "NA", 0, 0, "NA"]]),
          "unexpected geNomad topology")
    fails(lambda r: tsv(r / "GCA_2" / "GCA_2_summary" / "GCA_2_virus_summary.tsv",
                        ["seq_name", "length", "n_genes", "virus_score"], []), "lacks column")
    fails(lambda r: tsv(r / "000123" / "rloop.tsv", ["genome_id", "provirus_id", "length", "n_rlfs", "rloop_prone"],
                        [["000123", "x|provirus_1_9000", 9000, 0, "False"]] * 2), "duplicated provirus_id")


# =============================================================================== ecology_stats.R
def _tree_newick(tips, rng):
    items = [(t, 0.0) for t in tips]
    h = 0.0
    while len(items) > 1:
        h += rng.expovariate(len(items) * (len(items) - 1) / 2.0)
        i, j = sorted(rng.sample(range(len(items)), 2))
        a, b = items[i], items[j]
        new = (f"({a[0]}:{h - a[1]:.6f},{b[0]}:{h - b[1]:.6f})", h)
        items = [x for k, x in enumerate(items) if k not in (i, j)] + [new]
    return items[0][0] + ";\n"


def _stats_fixture(n=80, n_b_only=4, size_scale=1.0, niche_cols=True, seed=7):
    """Known-answer data: niche_A drives presence (90 % concordant, scattered over the
    tree); niche_B is EXACTLY balanced against presence (log-odds 0 before adjustment);
    group_rare has 2 ones, niche_const none, niche_sep = presence (separation). The last
    n_b_only genomes are B-only (validated locus, best tier B)."""
    rng = random.Random(seed)
    d = tmpdir("estats_")
    tips = [f"G{i:03d}" for i in range(n)]
    write(d / "tree.nwk", _tree_newick(tips, rng))
    rows = []
    for k, t in enumerate(tips):
        a = int(rng.random() < 0.5)
        pres = a if rng.random() < 0.9 else 1 - a
        rows.append(dict(genome_id=t, tree_tip=t, pass_primary="True",
                         prometheus_present=str(bool(pres)), prometheus_tier="A" if pres else "",
                         completeness=f"{rng.uniform(90, 100):.2f}",
                         genome_size=f"{round(rng.uniform(1.5e6, 4e6)) * size_scale:.12g}",
                         methanogen=str(rng.random() < 0.4), gtdb_phylum=rng.choice(["p__A", "p__B"]),
                         niche_A=str(a), niche_B="0", niche_const="0", group_rare=str(int(k < 2)),
                         niche_sep=str(pres)))
    for pres in (0, 1):                      # niche_B: half of each presence class
        idx = [i for i, r in enumerate(rows[:n - n_b_only]) if (r["prometheus_tier"] == "A") == bool(pres)]
        rng.shuffle(idx)
        for i in idx[: len(idx) // 2]:
            rows[i]["niche_B"] = "1"
    for r in rows[n - n_b_only:]:
        r["prometheus_present"], r["prometheus_tier"] = "True", "B"
    cols = list(rows[0].keys())
    if not niche_cols:
        cols = [c for c in cols if not c.startswith(("niche_", "group_"))]
    tsv(d / "master.tsv", cols, [[r[c] for c in cols] for r in rows])
    vrows = []
    for r in rows:
        pres = r["prometheus_tier"] in ("A", "B")
        nprov = rng.choice([0, 1, 2, 3, 4, 6] if pres else [0, 0, 1, 1, 2])
        nr = sum(rng.random() < (0.7 if pres else 0.3) for _ in range(nprov))
        vrows.append([r["genome_id"], nprov, nprov * 35.0, "", nr, f"{nr / nprov:.4f}" if nprov else "",
                      nprov + rng.choice([0, 0, 1]), f"{rng.uniform(20, 50):.1f}" if nprov else "",
                      f"{rng.uniform(0.35, 0.6):.4f}" if nprov else ""])
    tsv(d / "viral.tsv", ["genome_id", "n_proviruses", "total_provirus_kb", "mean_checkv_completeness",
                          "n_rloop_proviruses", "frac_rloop_proviruses", "n_viral_sequences",
                          "mean_provirus_kb", "mean_provirus_gc"], vrows)
    n_prov_pos = sum(1 for v, r in zip(vrows, rows) if v[1] > 0 and r["prometheus_tier"] != "B")
    return d, n_prov_pos


def _stats(d, *args, check=True, out="assoc.tsv"):
    p = rscript("ecology_stats.R", "--master", d / "master.tsv", "--tree", d / "tree.nwk",
                *args, "--out", d / out, check=check)
    return p, (read_tsv(d / out) if (d / out).exists() and p.returncode == 0 else None)


def _bh(p):
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    adj = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * m / rank)
        adj[i] = prev
    return adj


_KNOWN = {}
NOFIT = ("--min-n", "1000")     # counts, labels and predictor resolution only (no fit needed)


def _known_run():
    """The full known-answer run (all predictors + viral), computed once per suite."""
    if not _KNOWN:
        d, n_prov_pos = _stats_fixture()
        _, a = _stats(d, "--viral", d / "viral.tsv")
        _KNOWN.update(d=d, n_prov_pos=n_prov_pos, assoc=a)
    return _KNOWN["d"], _KNOWN["n_prov_pos"], _KNOWN["assoc"].copy()


@suite.test
def test_ecology_stats_known_answers():
    d, n_prov_pos, a = _known_run()
    a = a.set_index("predictor")
    assert list(a.index) == ["group_rare", "niche_A", "niche_B", "niche_const", "niche_sep",
                             "n_proviruses", "frac_rloop_proviruses"], list(a.index)
    A = a.loc["niche_A"]
    assert fnum(A["estimate"]) > 2 and fnum(A["p_value"]) < 1e-4 and fnum(A["ci_low"]) > 0, dict(A)
    assert A["level"] == "1" and A["n"] == "76" and A["n_excluded_tier"] == "4"
    assert "z(completeness) + z(genome_size)" in A["model"]
    B = a.loc["niche_B"]                      # exactly balanced null (glm log-OR = 0)
    assert fnum(B["p_value"]) > 0.5 and abs(fnum(B["estimate"])) < 0.5, dict(B)
    assert B["converged"] == "TRUE"
    assert a.loc["group_rare", "estimate"] == "NA" and "predictor 1=2" in a.loc["group_rare", "note"]
    assert a.loc["niche_const", "estimate"] == "NA" and "predictor 1=0" in a.loc["niche_const", "note"]
    S = a.loc["niche_sep"]
    assert "quasi-separation" in S["note"] and math.isfinite(fnum(S["estimate"])) and fnum(S["estimate"]) > 3
    F = a.loc["frac_rloop_proviruses"]
    assert F["level"] == "continuous" and int(F["n"]) == n_prov_pos, (F["n"], n_prov_pos)
    assert fnum(a.loc["n_proviruses", "estimate"]) > 0 and a.loc["n_proviruses", "n"] == "76"
    assert "per SD of log1p(n_proviruses)" in a.loc["n_proviruses", "note"]
    assert "methanogen" in A["sens_model"] and math.isfinite(fnum(A["sens_estimate"]))
    assert "viral contigs" in a.loc["n_proviruses", "sens_model"]
    assert "provirus GC" in F["sens_model"]
    pv = [(k, fnum(v)) for k, v in a["p_value"].items() if math.isfinite(fnum(v))]
    adj = _bh([v for _, v in pv])
    for (k, _), q in zip(pv, adj):
        assert abs(fnum(a.loc[k, "p_adj"]) - q) < 1e-9, k
    assert a.loc["group_rare", "p_adj"] == "NA"


@suite.test
def test_ecology_stats_nonconvergence_never_reported():
    """phyloglm MPLE stops with L-BFGS-B code 52 at a wrong beta on this exactly balanced
    null (it reported log-OR 1.57, p=0.0016). Only converged fits may be reported: here
    every MPLE start fails, so the flagged Ives-Garland (IG10) estimate (~0) is used."""
    d, _ = _stats_fixture()
    _, a = _stats(d, "--predictors", "niche_B", "--covariates", "")
    r = a.iloc[0]
    assert abs(fnum(r["estimate"])) < 0.25 and fnum(r["p_value"]) > 0.5, dict(r)
    assert r["model"] == "phyloglm(present ~ niche_B, logistic_IG10 fallback)", r["model"]
    assert "MPLE did not converge from 5 starts" in r["note"] and r["converged"] == "TRUE"
    # a strong effect converges from every start and keeps MPLE
    _, b = _stats(d, "--predictors", "niche_A")
    assert b.loc[0, "model"].endswith("logistic_MPLE)") and "MPLE converged from 5/5 starts" in b.loc[0, "note"]
    # a single start reproduces the old behaviour only when it converges; never a code-52 fit
    _, c = _stats(d, "--predictors", "niche_B", "--covariates", "", "--mple-starts", "1")
    assert abs(fnum(c.loc[0, "estimate"])) < 0.25 and "IG10" in c.loc[0, "model"]


@suite.test
def test_ecology_stats_na_niche_values():
    """Blank group_/niche_ cells are unknown: dropped from that predictor and reported."""
    d, _ = _stats_fixture()
    m = read_tsv(d / "master.tsv")
    m.loc[m.index[:5], "niche_A"] = ""
    m.to_csv(d / "master.tsv", sep="\t", index=False)
    _, a = _stats(d, "--predictors", "niche_A,niche_B", *NOFIT)
    a = a.set_index("predictor")
    n_drop = sum(1 for i in range(5) if m.loc[m.index[i], "prometheus_tier"] != "B")
    assert int(a.loc["niche_A", "n"]) == 76 - n_drop and a.loc["niche_B", "n"] == "76"
    assert f"{5} genome(s) with unknown (NA) niche_A dropped" in a.loc["niche_A", "note"], a.loc["niche_A", "note"]


@suite.test
def test_ecology_stats_presence_rule_b_only():
    """Contract v2 item 6: a B-only genome is EXCLUDED under T={A}, PRESENT under T={A,B}."""
    d, _ = _stats_fixture()
    _, a1 = _stats(d, "--predictors", "niche_B", *NOFIT, out="a.tsv")
    _, a2 = _stats(d, "--predictors", "niche_B", "--primary-tiers", "A,B", *NOFIT, out="ab.tsv")
    r1, r2 = a1.iloc[0], a2.iloc[0]
    assert (r1["n"], r1["n_excluded_tier"]) == ("76", "4"), dict(r1)
    assert (r2["n"], r2["n_excluded_tier"]) == ("80", "0"), dict(r2)
    assert int(r2["n_present"]) - int(r1["n_present"]) == 4       # the 4 B-only genomes
    assert r1["n_absent"] == r2["n_absent"]                       # never coded absent
    assert "tiers=A;" in r1["note"] and "4 excluded by the tier rule" in r1["note"]
    assert "tiers=A+B" in r2["note"]


@suite.test
def test_ecology_stats_call_shapes():
    # bacterial niche call (no --viral): niche rows only
    d, _ = _stats_fixture()
    _, a = _stats(d, *NOFIT)
    assert a["predictor"].tolist() == ["group_rare", "niche_A", "niche_B", "niche_const", "niche_sep"]
    # bacterial viral call (no niche columns, --viral): exactly the two viral rows
    d2, _ = _stats_fixture(niche_cols=False)
    _, b = _stats(d2, "--viral", d2 / "viral.tsv", *NOFIT)
    assert b["predictor"].tolist() == ["n_proviruses", "frac_rloop_proviruses"]
    # no predictor at all -> header-only, exit 0
    p, c = _stats(d2)
    assert c is not None and c.empty and list(c.columns)[:9] == [
        "predictor", "level", "model", "estimate", "std_error", "p_value", "p_adj", "n", "note"]
    # explicit predictor spec: exact names + globs; a missing name gets an explicit row
    _, e = _stats(d, "--predictors", "niche_A,biome_envo", *NOFIT)
    assert e["predictor"].tolist() == ["biome_envo", "niche_A"]
    assert "not found" in e.set_index("predictor").loc["biome_envo", "note"]
    _, f = _stats(d, "--predictors", "niche_*", *NOFIT)
    assert all(x.startswith("niche_") for x in f["predictor"]) and len(f) == 4
    # empty species tree placeholder "();" -> header-only (genuinely empty input)
    write(d / "tree.nwk", "();\n")
    _, g = _stats(d)
    assert g is not None and g.empty


@suite.test
def test_ecology_stats_covariate_options():
    """A niche column used as a covariate (temperature-aware tests) is never its own
    covariate; an empty covariate list gives unadjusted models; unknown covariates are
    reported in the notes."""
    d, _ = _stats_fixture()
    _, a = _stats(d, "--predictors", "niche_A,niche_B", "--covariates", "completeness,niche_B,nope", *NOFIT)
    a = a.set_index("predictor")
    assert "z(niche_B)" in a.loc["niche_A", "model"], a.loc["niche_A", "model"]
    assert "z(niche_B)" not in a.loc["niche_B", "model"], a.loc["niche_B", "model"]
    assert "covariate nope not in master" in a.loc["niche_A", "note"]
    _, b = _stats(d, "--predictors", "niche_A", "--covariates", "", *NOFIT)
    assert b.loc[0, "model"] == "phyloglm(present ~ niche_A, logistic_MPLE)", b.loc[0, "model"]


@suite.test
def test_ecology_stats_scale_invariance():
    """Covariates are z-scored: genome_size in bp or Mb gives the same niche_A estimate."""
    d1, _ = _stats_fixture(size_scale=1.0)
    d2, _ = _stats_fixture(size_scale=1e-6)
    _, a1 = _stats(d1, "--predictors", "niche_A")
    _, a2 = _stats(d2, "--predictors", "niche_A")
    e1, e2 = fnum(a1.loc[0, "estimate"]), fnum(a2.loc[0, "estimate"])
    assert abs(e1 - e2) < 1e-3 * max(1.0, abs(e1)), (e1, e2)


@suite.test
def test_ecology_stats_failures():
    d, _ = _stats_fixture()
    unrooted = run([RSCRIPT, "-e", f'cat(ape::write.tree(ape::unroot(ape::read.tree("{d / "tree.nwk"}"))))']).stdout
    assert unrooted.count(",") == 79 and unrooted.endswith(";")
    bad = [["--master", d / "nope.tsv"], ["--tree", write(d / "bad.nwk", "((a,b);\n")],
           ["--tree", write(d / "unrooted.nwk", unrooted + "\n")],
           ["--primary-tiers", "Z"], ["--viral", write(d / "v.tsv", "genome_id\tfoo\nG001\t1\n")],
           ["--tree", write(d / "nan.nwk", "((G000:NaN,G001:1):1,G002:2);\n")],
           # zero branch lengths: parseable, but every phyloglm fit fails -> error, no output
           ["--tree", write(d / "zero.nwk", __import__("re").sub(r":[0-9.]+", ":0",
                                                                  (d / "tree.nwk").read_text()))]]
    import subprocess
    from concurrent.futures import ThreadPoolExecutor

    def one(job):
        i, extra = job
        args = {"--master": d / "master.tsv", "--tree": d / "tree.nwk", "--predictors": "niche_A"}
        it = iter(extra)
        for k in it:
            args[k] = next(it)
        flat = [str(x) for kv in args.items() for x in kv]
        out = d / f"x{i}.tsv"
        p = subprocess.run([RSCRIPT, str(SCRIPTS_R / "ecology_stats.R"), *flat, "--out", str(out)],
                           env=env(), capture_output=True, text=True, timeout=300)
        return extra, p.returncode, out.exists()
    with ThreadPoolExecutor(max_workers=2) as pool:            # <= 2 concurrent R processes
        for extra, rc, made in pool.map(one, list(enumerate(bad))):
            assert rc != 0 and not made, (extra, rc, made)
    tsv(d / "noid.tsv", ["gid", "prometheus_tier"], [["G001", "A"]])
    p = rscript("ecology_stats.R", "--master", d / "noid.tsv", "--tree", d / "tree.nwk",
                "--out", d / "y.tsv", check=False)
    assert p.returncode != 0 and "genome_id" in p.stderr
    # a missing R package is an error, not a skipped test: hide the extra R library
    probe = run([RSCRIPT, "-e", 'cat(requireNamespace("phylolm", quietly=TRUE))'], check=False)
    e = dict(env()); e["R_LIBS"] = ""; e["R_LIBS_USER"] = "/nonexistent"
    import subprocess
    q = subprocess.run([RSCRIPT, "-e", 'cat(requireNamespace("phylolm", quietly=TRUE))'], env=e,
                       capture_output=True, text=True)
    if probe.stdout.strip() == "TRUE" and q.stdout.strip() == "FALSE":
        r = subprocess.run([RSCRIPT, str(SCRIPTS_R / "ecology_stats.R"), "--master", d / "master.tsv",
                            "--tree", d / "tree.nwk", "--out", d / "z.tsv"], env=e,
                           capture_output=True, text=True)
        assert r.returncode != 0 and not (d / "z.tsv").exists(), r.stderr[-400:]


# =============================================================================== plots
@suite.test
def test_plots_known_inputs():
    d, _, _ = _known_run()
    out = tmpdir("plots_")
    p = run_py("plot_ecology.py", "--associations", d / "assoc.tsv", "--viral", d / "viral.tsv",
               "--master", d / "master.tsv", "--outdir", out, "--formats", "png")
    assert "ecology_enrichment: 3 predictor(s) plotted" in p.stderr, p.stderr
    assert "42 present, 34 absent, 4 excluded" in p.stderr or "4 excluded" in p.stderr, p.stderr
    for f in ("ecology_enrichment.png", "viral_rloop_association.png"):
        assert (out / f).stat().st_size > 1000, f
    assert not (out / "ecology_enrichment.svg").exists()               # --formats honoured
    p = run_py("plot_ecology.py", "--associations", d / "assoc.tsv", "--viral", d / "viral.tsv",
               "--master", d / "master.tsv", "--outdir", out / "ab", "--primary-tiers", "A,B")
    assert "0 excluded" in p.stderr and (out / "ab" / "viral_rloop_association.svg").exists()
    p = run_py("plot_ecology_extra.py", "--master", d / "master.tsv", "--viral", d / "viral.tsv",
               "--associations", d / "assoc.tsv", "--outdir", out / "x")
    assert "'excluded': 4" in p.stderr and "niche_effect_forest: 3 predictor(s)" in p.stderr, p.stderr
    for f in ("viral_pressure_prometheus", "niche_effect_forest", "prevalence_niche_phylum"):
        for ext in ("png", "svg"):
            assert (out / "x" / f"{f}.{ext}").stat().st_size > 1000, f


@suite.test
def test_plots_empty_and_missing_inputs():
    d = tmpdir("plot_e_")
    tsv(d / "assoc.tsv", ["predictor", "level", "model", "estimate", "std_error", "p_value",
                          "p_adj", "n", "note"], [])
    tsv(d / "viral.tsv", ["genome_id", "n_proviruses", "frac_rloop_proviruses"], [])
    tsv(d / "master.tsv", ["genome_id", "prometheus_present", "prometheus_tier"], [])
    run_py("plot_ecology.py", "--associations", d / "assoc.tsv", "--viral", d / "viral.tsv",
           "--master", d / "master.tsv", "--outdir", d / "p")
    run_py("plot_ecology_extra.py", "--associations", d / "assoc.tsv", "--viral", d / "viral.tsv",
           "--master", d / "master.tsv", "--outdir", d / "p")
    for f in ("ecology_enrichment", "viral_rloop_association", "viral_pressure_prometheus",
              "niche_effect_forest", "prevalence_niche_phylum"):
        assert (d / "p" / f"{f}.png").exists() and (d / "p" / f"{f}.svg").exists(), f
    p = run_py("plot_ecology.py", "--associations", d / "missing.tsv", "--viral", d / "viral.tsv",
               "--master", d / "master.tsv", "--outdir", d / "q", check=False)
    assert p.returncode != 0
    p = run_py("plot_ecology_extra.py", "--associations", d / "assoc.tsv", "--viral", d / "nope.tsv",
               "--master", d / "master.tsv", "--outdir", d / "q", check=False)
    assert p.returncode != 0


# =============================================================================== rules (Snakemake)
STUB_GENOMAD = r'''#!/usr/bin/env python3
"""Stub of `genomad end-to-end` writing geNomad 1.8.0-format outputs (verified against
genomad/modules/summary.py): {out}/{stem}_summary/{stem}_virus_summary.tsv + _virus.fna."""
import os, sys
args = [a for a in sys.argv[1:]]
assert args[0] == "end-to-end" and {"--cleanup", "--restart", "--threads", "--splits"} <= set(args), args
pos = [a for a in args[1:] if not a.startswith("--") and not a.isdigit()]
inp, out, db = pos[0], pos[1], pos[2]
stem = os.path.basename(inp).rsplit(".", 1)[0]
if stem in os.environ.get("STUB_GENOMAD_FAIL", "").split(","):
    print("Traceback: simulated geNomad crash", file=sys.stderr); sys.exit(1)
if not os.path.isfile(os.path.join(db, "version.txt")):
    sys.exit("database missing")
seqs, name = {}, None
for line in open(inp):
    line = line.strip()
    if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = []
    elif line: seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
sd = os.path.join(out, f"{stem}_summary"); os.makedirs(sd, exist_ok=True)
hdr = "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
with open(os.path.join(sd, f"{stem}_virus_summary.tsv"), "w") as s, open(os.path.join(sd, f"{stem}_virus.fna"), "w") as f:
    s.write(hdr)
    if "novirus" not in stem:
        c1, c2 = list(seqs)[:2]
        prov = seqs[c1][100:2100]
        pn = f"{c1}|provirus_101_2100"
        s.write(f"{pn}\t{len(prov)}\tProvirus\t101-2100\t3\t11\t0.9500\tNA\t1\t2.0\tViruses\n")
        f.write(f">{pn}\n{prov}\n")
        s.write(f"{c2}\t{len(seqs[c2])}\tNo terminal repeats\tNA\t2\t11\t0.9100\tNA\t1\t1.5\tViruses\n")
        f.write(f">{c2} some description\n{seqs[c2]}\n")
'''

STUB_CHECKV = r'''#!/usr/bin/env python3
"""Stub of `checkv end_to_end` writing the CheckV 1.0.3 quality_summary.tsv columns. With
STUB_CHECKV_EXIT0=1 it mimics CheckV's bare sys.exit() paths: exit 0, no output."""
import os, sys
a = sys.argv[1:]
assert a[0] == "end_to_end" and {"--restart", "--remove_tmp", "-d", "-t"} <= set(a), a
inp, out = a[1], a[2]
db = a[a.index("-d") + 1]
for f in ("genome_db", "hmm_db"):
    if not os.path.isdir(os.path.join(db, f)): sys.exit(f"Error: database dir not found '{db}'")
if os.environ.get("STUB_CHECKV_EXIT0") == "1":
    sys.stderr.write("You input FASTA file is empty or not properly formatted."); sys.exit()
os.makedirs(out, exist_ok=True)
ids = [l[1:].split()[0] for l in open(inp) if l.startswith(">")]
cols = ["contig_id","contig_length","provirus","proviral_length","gene_count","viral_genes","host_genes",
        "checkv_quality","miuvig_quality","completeness","completeness_method","contamination","kmer_freq","warnings"]
with open(os.path.join(out, "quality_summary.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for i, c in enumerate(ids):
        fh.write("\t".join([c, "2000", "No", "NA", "3", "2", "0", "Medium-quality" if i == 0 else "Not-determined",
                            "Genome-fragment", "62.5" if i == 0 else "NA", "AAI-based", "0.0", "1.0", ""]) + "\n")
'''

MINI_SNAKEFILE = r'''
import os
configfile: "config.yaml"
RESULTS, LOGS = "results", "logs"
GENOMES = [l.split("\t")[0] for l in open(config["samples"]).read().splitlines()[1:] if l.strip()]
PLOT_FORMATS = ["png"]
T_SMALL, T_MED = 1, 2
def cfg(path, default=None):
    node = config
    for k in path.split("/"):
        if isinstance(node, dict) and k in node: node = node[k]
        else: return default
    return node
def genome_fasta(wc): return f"genomes/{wc.genome}.fna"
def mem_mb(wc, attempt, base=4000, per_attempt=4000): return base + (attempt - 1) * per_attempt
rule all:
    input:
        f"{RESULTS}/ecology/ecology_associations.tsv",
        f"{RESULTS}/ecology/plots/ecology_enrichment.png",
rule build_master_table:
    input: "master_src.tsv"
    output: master=f"{RESULTS}/scaffolds/master_table.tsv"
    shell: "cp {input} {output.master}"
rule species_tree:
    input: "tree_src.nwk"
    output: tree=f"{RESULTS}/scaffolds/species_tree.treefile"
    shell: "cp {input} {output.tree}"
include: "@ECOLOGY_SMK@"
'''


def _mini_workflow(n=12, novirus=("G03",)):
    w = tmpdir("eco_smk_")
    (w / "workflow").mkdir()
    os.symlink(REPO / "workflow" / "scripts", w / "workflow" / "scripts")
    bindir = w / "bin"; bindir.mkdir()
    for name, body in (("genomad", STUB_GENOMAD), ("checkv", STUB_CHECKV)):
        f = bindir / name
        f.write_text(f"#!{PYTHON}\n" + body.split("\n", 1)[1])
        f.chmod(f.stat().st_mode | stat.S_IEXEC)
    (w / "db" / "genomad_db").mkdir(parents=True); write(w / "db" / "genomad_db" / "version.txt", "1.7\n")
    for sub in ("genome_db", "hmm_db"):
        (w / "db" / "checkv-db" / sub).mkdir(parents=True)
    rng = random.Random(3)
    ids = [f"G{i:02d}" for i in range(1, n + 1)]
    rows = []
    for k, g in enumerate(ids):
        stem = f"{g}novirus" if g in novirus else g
        c1 = "".join(rng.choice("ACGGT") for _ in range(3000))
        c2 = "".join(rng.choice("ACGT") for _ in range(1500))
        write(w / "genomes" / f"{stem}.fna", f">{stem}_c1\n{c1}\n>{stem}_c2\n{c2}\n")
        tier = "A" if k % 2 == 0 else ""
        rows.append([stem, stem, "True", str(bool(tier)), tier, f"{90 + k % 7}", str(2_000_000 + 37_000 * k),
                     str(int(k % 3 == 0)), str(int(k % 4 == 0))])
    tsv(w / "master_src.tsv", ["genome_id", "tree_tip", "pass_primary", "prometheus_present",
                               "prometheus_tier", "completeness", "genome_size", "niche_x", "group_y"], rows)
    write(w / "tree_src.nwk", _tree_newick([r[0] for r in rows], rng))
    tsv(w / "samples.tsv", ["genome_id"], [[r[0]] for r in rows])
    write(w / "config.yaml",
          f"samples: samples.tsv\nseed: 1\ndb:\n  genomad: {w}/db/genomad_db\n  checkv: {w}/db/checkv-db\n"
          "stats:\n  confidence_tiers_primary: [A]\n")
    # The real ecology.smk with only its `conda:` directives removed: Snakemake 8.30 touches
    # the conda wrapper of a FAILED job even without --use-conda and then hangs when conda
    # is not on PATH. Shell commands, params, logs and outputs are byte-identical.
    import re
    src = (REPO / "workflow/rules/ecology.smk").read_text()
    stripped, n_conda = re.subn(r'\n[ ]+conda:[ ]*\n[ ]+"\.\./envs/[a-z_]+\.yaml"[ ]*(?=\n)', "", src)
    # 7 = genomad, checkv, rloop_predict, viral_summary, ecology_stats, methanogenesis_stats,
    # plot_ecology
    assert n_conda == 7, f"expected 7 conda directives in ecology.smk, found {n_conda}"
    write(w / "workflow" / "rules" / "ecology.smk", stripped)
    write(w / "Snakefile", MINI_SNAKEFILE.replace("@ECOLOGY_SMK@", str(w / "workflow/rules/ecology.smk")))
    return w, bindir, [r[0] for r in rows]


def _smk(w, bindir, target="all", extra_env=None, config=None):
    import subprocess
    e = dict(env())
    e["PATH"] = os.pathsep.join([str(bindir), str(Path(PYTHON).parent), str(Path(RSCRIPT).parent), e["PATH"]])
    e["HOME"] = str(w / "home"); e["XDG_CACHE_HOME"] = str(w / "home" / ".cache")
    e["MPLCONFIGDIR"] = str(w / "home" / "mpl"); e["PROM_DPI"] = "50"
    if extra_env:
        e.update(extra_env)
    # --drop-metadata: Snakemake 8 otherwise hashes every conda: env file after each job
    # (persistence._software_stack_hash -> `conda info`), which needs conda on PATH even
    # without --use-conda; the stubs replace the conda envs here.
    cmd = [SMK, target, "-s", str(w / "Snakefile"), "-d", str(w), "--cores", "2", "--nolock",
           "--keep-going", "--drop-metadata"]
    if config:
        cmd += ["--config"] + config
    return subprocess.run(cmd, cwd=w, env=e, text=True, capture_output=True, timeout=300)


def _have_snakemake():
    return os.path.exists(SMK) or shutil.which("snakemake") is not None


@suite.test
def test_rules_run_end_to_end_with_stub_tools():
    """ecology.smk executed for real on 12 genomes with geNomad/CheckV stubs."""
    assert _have_snakemake(), f"snakemake not found ({SMK}); set PROM_SNAKEMAKE"
    w, bindir, ids = _mini_workflow()
    p = _smk(w, bindir)
    assert p.returncode == 0, p.stderr[-4000:]
    v = read_tsv(w / "results/ecology/viral/viral_summary.tsv").set_index("genome_id")
    assert list(v.index) == ids
    normal = v.loc["G01"]
    assert (normal["n_proviruses"], normal["n_viral_contigs"], normal["n_viral_sequences"]) == ("1", "1", "2")
    assert fnum(normal["total_provirus_kb"]) == 2.0 and fnum(normal["mean_checkv_completeness"]) == 62.5
    empty = v.loc["G03novirus"]
    assert (empty["n_proviruses"], empty["frac_rloop_proviruses"]) == ("0", "")
    qs = (w / "results/ecology/viral/G03novirus/checkv/quality_summary.tsv").read_text().splitlines()
    assert qs == ["\t".join(CHECKV_COLS)], qs                        # real CheckV header
    assert (w / "results/ecology/viral/G01/rloop_rlfs_calls.tsv").exists()
    for f in ("ecology_associations.tsv", "ecology_associations_sensitivity.tsv"):
        a = read_tsv(w / "results/ecology" / f)
        assert {"group_y", "niche_x", "n_proviruses", "frac_rloop_proviruses"} <= set(a["predictor"]), f
    for log in ("genomad/G01.log", "checkv/G01.log", "rloop/G01.log", "ecology/viral_summary.log",
                "ecology/ecology_stats.log", "ecology/plot_ecology.log"):
        assert (w / "logs" / log).exists(), log
    assert "sensitivity tiers A,B" in (w / "logs/ecology/ecology_stats.log").read_text()


@suite.test
def test_rules_fail_loudly():
    assert _have_snakemake(), f"snakemake not found ({SMK}); set PROM_SNAKEMAKE"
    # a crashed geNomad fails its job and nothing downstream is produced
    w, bindir, _ = _mini_workflow(n=4)
    p = _smk(w, bindir, extra_env={"STUB_GENOMAD_FAIL": "G02"})
    assert p.returncode != 0
    assert not (w / "results/ecology/viral/viral_summary.tsv").exists()
    assert not (w / "results/ecology/viral/G02/checkv/quality_summary.tsv").exists()
    assert "simulated geNomad crash" in (w / "logs/genomad/G02.log").read_text()
    # missing geNomad database
    w, bindir, _ = _mini_workflow(n=4)
    shutil.rmtree(w / "db" / "genomad_db")
    p = _smk(w, bindir)
    assert p.returncode != 0 and "geNomad database not found" in (w / "logs/genomad/G01.log").read_text()
    # CheckV exiting 0 without writing results (its bare sys.exit() paths)
    w, bindir, _ = _mini_workflow(n=4)
    p = _smk(w, bindir, extra_env={"STUB_CHECKV_EXIT0": "1"})
    assert p.returncode != 0
    assert "CheckV exited 0 but wrote no quality_summary.tsv" in (w / "logs/checkv/G01.log").read_text()
    # missing CheckV database
    w, bindir, _ = _mini_workflow(n=4)
    shutil.rmtree(w / "db" / "checkv-db" / "hmm_db")
    p = _smk(w, bindir)
    assert p.returncode != 0 and "CheckV database not found" in (w / "logs/checkv/G01.log").read_text()


if __name__ == "__main__":
    suite.main("ECOLOGY TESTS PASSED")
