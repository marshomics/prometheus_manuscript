"""Synthetic test suite for stage preprocess_qc (leaf-1.1.1): Phase-1 preprocessing
(batch_genomes, gtdbtk_batch, checkm2_batch, bbmap_batch, merge_annotations,
filter_collapse, plot_preprocess) and Stage-0 gene calls / QC (prodigal_genes +
rename_prodigal, genome_metadata, plot_qc, plot_qc_extra).

    bash tests/run_leaf.sh preprocess_qc        # prints PREPROCESS_QC TESTS PASSED

Every script is run on synthetic inputs with KNOWN answers (values are asserted, not just
file existence), on empty input, and on failure input (a crashed tool must fail).
Real tools: Prodigal 2.6.3 (PROM_TOOLS_BIN) is run on a planted-gene genome and at the
20-kb single-mode training boundary; Snakemake (PROM_SNAKEMAKE) runs the real rules of
preprocess.smk + qc.smk end to end, with format-faithful stand-ins for GTDB-Tk / CheckM2 /
statswrapper.sh (their databases are 10-100 GB). A real BBMap 39.06 is used when
PROM_BBMAP_STATSWRAPPER points at statswrapper.sh or statswrapper.sh is on the toolchain
PATH (java on PATH; the JVM is capped at 2 processors); otherwise the parser is still
tested on captured real BBMap 39.06 output.
"""
from __future__ import annotations
import gzip
import math
import os
import py_compile
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from _env import Suite, py, run, tmpdir, have, read_tsv, REPO, SCRIPTS_PY, PYTHON, env as tool_env  # noqa: E402
import synth  # noqa: E402

sys.path.insert(0, str(SCRIPTS_PY))
suite = Suite("preprocess_qc")
OWNED_PY = ["batch_genomes.py", "merge_annotations.py", "filter_collapse.py", "run_bbmap_stats.py",
            "plot_preprocess.py", "rename_prodigal.py", "plot_qc.py", "plot_qc_extra.py"]
RULES = {"preprocess": REPO / "workflow" / "rules" / "preprocess.smk",
         "qc": REPO / "workflow" / "rules" / "qc.smk"}
_TMP = []


def tdir(prefix):
    d = tmpdir(f"pqc_{prefix}_")
    _TMP.append(d)
    return d


def write_tsv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


def run_env(cmd, extra=None, path_prepend=(), cwd=None, check=True, timeout=900):
    """Like _env.run but with extra environment variables / PATH entries."""
    e = tool_env()
    e.update({k: str(v) for k, v in (extra or {}).items()})
    if path_prepend:
        e["PATH"] = os.pathsep.join([str(p) for p in path_prepend] + [e["PATH"]])
    p = subprocess.run(cmd, cwd=cwd or REPO, env=e, text=True, capture_output=True, timeout=timeout,
                       shell=isinstance(cmd, str), executable="/bin/bash" if isinstance(cmd, str) else None)
    if check and p.returncode != 0:
        raise AssertionError(f"command failed ({p.returncode}): {cmd}\n{p.stdout[-2500:]}\n{p.stderr[-2500:]}")
    return p


def expect_fail(proc, pattern):
    out = proc.stdout + proc.stderr
    assert proc.returncode != 0, f"expected failure, got rc=0:\n{out[-1500:]}"
    assert re.search(pattern, out, flags=re.I), f"failure message lacks /{pattern}/:\n{out[-1500:]}"


def rand_dna(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def rule_block(name):
    text = RULES["preprocess"].read_text() + "\n" + RULES["qc"].read_text()
    m = re.search(rf"^rule {name}:\n(.*?)(?=^rule |\Z)", text, flags=re.S | re.M)
    assert m, f"rule {name} not found"
    return m.group(1)


def shell_awk(rule, start_marker):
    """Extract an awk program (single-quoted) from a rule's shell block; undo Snakemake's
    brace doubling."""
    blk = rule_block(rule)
    i = blk.index(start_marker) + len(start_marker)
    j = blk.index("'", i)
    return blk[i:j].replace("{{", "{").replace("}}", "}")


# =========================================================================== #
# static checks
# =========================================================================== #
@suite.test
def test_py_compile_owned_scripts():
    d = tdir("pyc")
    for name in OWNED_PY:
        py_compile.compile(str(SCRIPTS_PY / name), cfile=str(d / (name + "c")), doraise=True)
    for f in ("run_tests.py", "synth.py"):
        py_compile.compile(str(HERE / f), cfile=str(d / (f + "c")), doraise=True)
    # aggregate_qc.py was dead code (no rule called it) and has been removed.
    assert not (SCRIPTS_PY / "aggregate_qc.py").exists()


@suite.test
def test_rule_hygiene_and_static_contract():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"), str(RULES["preprocess"]), str(RULES["qc"])])
    assert "RULES HYGIENE OK" in p.stdout
    pre = RULES["preprocess"].read_text()
    qc = RULES["qc"].read_text()
    for t in (pre, qc):
        assert "|| true" not in t and "|| echo 0" not in t
    assert "find " not in rule_block("gtdbtk_batch")
    assert "classify/gtdbtk.ar53.summary.tsv" in rule_block("gtdbtk_batch")
    assert "--batchfile {input.batchfile}" in rule_block("gtdbtk_batch")
    assert '"../envs/checkm2.yaml"' in rule_block("checkm2_batch")
    assert '"../envs/qc.yaml"' in rule_block("prodigal_genes")
    assert 'genes=f"{RESULTS}/genes/{{genome}}/{{genome}}.genes.tsv"' in rule_block("prodigal_genes")
    assert "-g {params.table}" in rule_block("prodigal_genes")
    assert 'cfg("qc/prodigal_translation_table", 11)' in rule_block("prodigal_genes")
    assert 'cfg("preprocess/bbmap/length_col", "contig_bp")' in rule_block("merge_annotations")
    # statswrapper.sh = one single-threaded JVM per chunk, chunks run sequentially (PQ-35)
    assert re.search(r"^\s*threads: T_SMALL\s*$", rule_block("bbmap_batch"), re.M)
    envs = REPO / "workflow" / "envs"
    qc_env, cm_env, bb_env = ((envs / f).read_text() for f in ("qc.yaml", "checkm2.yaml", "bbmap.yaml"))
    assert re.search(r"^\s*- prodigal=2\.6\.3\s*$", qc_env, re.M) and "checkm2" not in re.sub(r"#.*", "", qc_env)
    assert re.search(r"^\s*- checkm2=1\.0\.2\s*$", cm_env, re.M)
    assert re.search(r"^\s*- bbmap=39\.06\b", bb_env, re.M)


# =========================================================================== #
# batch_genomes.py
# =========================================================================== #
def _count_like_common_smk(path, id_col="genome_id"):
    """Execute common.smk's own _count_input_genomes on `path` (drift guard)."""
    import pandas as pd
    src = (REPO / "workflow" / "rules" / "common.smk").read_text()
    # every top-level import, helper def and private constant (Snakemake directives skipped)
    blocks = re.findall(r"^(?:def \w+\(.*?\):\n(?:[ \t]+.*\n|\n)+|_[A-Z_]+ = .*\n"
                        r"|import [\w.]+(?: as \w+)?\n)", src, flags=re.M)   # not `from snakemake...`
    ns = {"os": os, "pd": pd, "math": math, "ID_COL": id_col}
    exec("".join(blocks), ns)
    assert "_count_input_genomes" in ns, "common.smk no longer defines _count_input_genomes"
    return ns["_count_input_genomes"](str(path))


def _mk_inputs(root, ids, fname=lambda g: f"{g}.fasta", seq_len=500):
    gdir = Path(root) / "in"
    gdir.mkdir(parents=True, exist_ok=True)
    for i, g in enumerate(ids):
        synth.write_fasta(gdir / fname(g), [f"{g}_c1"], [rand_dna(seq_len, i)])
    return gdir


@suite.test
def test_batch_genomes_layout_matches_common_smk():
    d = tdir("batch")
    ids = ["G03", "G01", "G05", "G02", "G04", "G07", "G06"]
    gdir = _mk_inputs(d, ids, fname=lambda g: (f"{g}.fasta" if g == "G05" else f"x_{g}.fa"))
    meta = d / "meta.tsv"
    with open(meta, "w") as fh:
        fh.write("genome_id\tfilename\tsource\tgroup_gut\n# a comment line\n\n")
        for g in ids:
            fh.write(f"{g}\t{'' if g == 'G05' else 'x_' + g + '.fa'}\tNCBI\t1\n")
        fh.write("  \tx_none.fa\tNCBI\t0\n")          # empty id: not a genome
    n_common = _count_like_common_smk(meta)
    assert n_common == 7
    out = d / "batches"
    py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--input-ext", "fasta",
       "--batch-size", 3, "--expected-n", n_common, "--outdir", out)
    nb = math.ceil(n_common / 3)
    dirs = sorted(p.name for p in out.iterdir() if p.name.startswith("batch_"))
    assert dirs == [f"batch_{i:04d}" for i in range(nb)] == ["batch_0000", "batch_0001", "batch_0002"]
    want = {"batch_0000": ["G01", "G02", "G03"], "batch_0001": ["G04", "G05", "G06"], "batch_0002": ["G07"]}
    for b, gids in want.items():
        assert (out / b / ".ready").exists()
        links = sorted(os.listdir(out / b / "genomes"))
        assert links == [f"{g}.fna" for g in gids], links
        for g in gids:
            lk = out / b / "genomes" / f"{g}.fna"
            src = gdir / (f"{g}.fasta" if g == "G05" else f"x_{g}.fa")
            assert lk.is_symlink() and os.readlink(lk) == str(src.resolve())
        bf = (out / b / "gtdbtk_batchfile.tsv").read_text().splitlines()
        assert bf == [f"{os.path.abspath(out / b / 'genomes' / (g + '.fna'))}\tuser_{g}" for g in gids]
    man = read_tsv(out / "manifest.tsv")
    assert list(man.columns) == ["genome_id", "batch", "source_path", "staged_path", "gtdbtk_id"]
    assert man["genome_id"].tolist() == sorted(ids)
    assert man.set_index("genome_id").loc["G07", "batch"] == "0002"
    assert man.set_index("genome_id").loc["G05", "gtdbtk_id"] == "user_G05"


@suite.test
def test_batch_genomes_fails_on_bad_inputs():
    base_ids = ["A1", "A2"]

    def case(name, mutate, pattern, ids=base_ids, extra=()):
        d = tdir(f"bad_{name}")
        gdir = _mk_inputs(d, ids)
        meta = write_tsv(d / "meta.tsv", ["genome_id", "filename"], [[g, f"{g}.fasta"] for g in ids])
        mutate(d, gdir, meta)
        p = py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--batch-size", 5,
               "--outdir", d / "b", *extra, check=False)
        expect_fail(p, pattern)

    case("missing", lambda d, g, m: (g / "A2.fasta").unlink(), r"A2: .*missing")
    case("gzip", lambda d, g, m: (g / "A1.fasta").write_bytes(gzip.compress(b">x\nACGT\n")), r"A1: .*gzip")
    case("empty", lambda d, g, m: (g / "A1.fasta").write_text(""), r"A1: .*empty file")
    case("notfasta", lambda d, g, m: (g / "A1.fasta").write_text("ACGTACGT\n"), r"not FASTA")
    for bad in ("bad id", "G|1", "A__B", "_lead"):
        case(f"id{len(bad)}", lambda d, g, m, b=bad: (
            shutil.copy(g / "A1.fasta", g / f"{b}.fasta"),
            m.write_text(m.read_text() + f"{b}\t{b}.fasta\n")), r"invalid genome id")
    case("dup", lambda d, g, m: m.write_text(m.read_text() + "A1\tA1.fasta\n"), r"duplicated genome id")
    case("hash", lambda d, g, m: m.write_text(m.read_text().replace("A2\tA2.fasta", "A2\tA2.fasta#x")),
         r"inside a data row")
    case("count", lambda d, g, m: None, r"expected 3", extra=("--expected-n", 3))
    # GTDB-Tk 2.4.0 main.py _verify_file_path rejects batchfile paths containing a space:
    # batch_genomes must refuse such an outdir before staging anything.
    d = tdir("bad_space")
    gdir = _mk_inputs(d, base_ids)
    meta = write_tsv(d / "meta.tsv", ["genome_id", "filename"], [[g, f"{g}.fasta"] for g in base_ids])
    out = d / "with space" / "b"
    p = py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--batch-size", 5,
           "--outdir", out, check=False)
    expect_fail(p, r"contains a space; GTDB-Tk 2\.4\.0 rejects")
    assert not out.exists(), "nothing may be staged when the batch directory contains a space"


@suite.test
def test_batch_genomes_rerun_refreshes_only_its_links():
    d = tdir("rerun")
    ids = ["B1", "B2", "B3", "B4"]
    gdir = _mk_inputs(d, ids)
    meta = write_tsv(d / "meta.tsv", ["genome_id"], [[g] for g in ids])
    out = d / "batches"
    py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--batch-size", 10, "--outdir", out)
    assert sorted(os.listdir(out / "batch_0000" / "genomes")) == [f"{g}.fna" for g in ids]
    write_tsv(meta, ["genome_id"], [[g] for g in ids[:3]])
    py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--batch-size", 10, "--outdir", out)
    assert sorted(os.listdir(out / "batch_0000" / "genomes")) == [f"{g}.fna" for g in ids[:3]]
    (out / "batch_0000" / "genomes" / "notes.txt").write_text("user file")
    p = py("batch_genomes.py", "--metadata", meta, "--genomes-dir", gdir, "--batch-size", 10,
           "--outdir", out, check=False)
    expect_fail(p, r"non-symlink")


@suite.test
def test_batch_genomes_empty_metadata():
    d = tdir("bempty")
    meta = write_tsv(d / "meta.tsv", ["genome_id", "filename"], [])
    assert _count_like_common_smk(meta) == 0
    py("batch_genomes.py", "--metadata", meta, "--batch-size", 5, "--expected-n", 0, "--outdir", d / "b")
    assert read_tsv(d / "b" / "manifest.tsv").empty
    assert not [p for p in (d / "b").iterdir() if p.name.startswith("batch_")]


# =========================================================================== #
# gtdbtk_batch: the summary merge (awk program of the rule, run directly)
# =========================================================================== #
def _gtdbtk_summary(path, rows):
    hdr = synth.GTDBTK_HEADER
    write_tsv(path, hdr, [[g, c] + ["N/A"] * (len(hdr) - 2) for g, c in rows])


@suite.test
def test_gtdbtk_summary_merge_awk():
    prog = shell_awk("gtdbtk_batch", "LC_ALL=C awk -v pre=\"{params.prefix}\" '")
    d = tdir("awk")
    ar = d / "classify" / "gtdbtk.ar53.summary.tsv"
    bac = d / "classify" / "gtdbtk.bac120.summary.tsv"
    _gtdbtk_summary(ar, [("user_GCA_000001.1", "d__Archaea;p__Methanobacteriota;s__X y"),
                         ("user_G2", "Unclassified Archaea")])
    _gtdbtk_summary(bac, [("user_G3", "Unclassified")])
    # GTDB-Tk's own root-level symlink and markers summary must NOT be read by the rule.
    os.symlink("classify/gtdbtk.ar53.summary.tsv", d / "gtdbtk.ar53.summary.tsv")
    p = run(["awk", "-v", "pre=user_", prog, str(ar), str(bac)], cwd=d)
    lines = p.stdout.splitlines()
    assert lines[0].split("\t") == synth.GTDBTK_HEADER + ["marker_set"]
    rows = [ln.split("\t") for ln in lines[1:]]
    assert [(r[0], r[1], r[-1]) for r in rows] == [
        ("GCA_000001.1", "d__Archaea;p__Methanobacteriota;s__X y", "ar53"),
        ("G2", "Unclassified Archaea", "ar53"), ("G3", "Unclassified", "bac120")]
    assert all(len(r) == len(synth.GTDBTK_HEADER) + 1 for r in rows)
    # negative controls: header mismatch, an id without the batch prefix
    bad = d / "classify" / "bad.bac120.summary.tsv"
    write_tsv(bad, ["user_genome", "classification", "other"], [["user_G9", "x", "y"]])
    assert run(["awk", "-v", "pre=user_", prog, str(ar), str(bad)], cwd=d, check=False).returncode != 0
    _gtdbtk_summary(bac, [("G3", "Unclassified")])
    assert run(["awk", "-v", "pre=user_", prog, str(ar), str(bac)], cwd=d, check=False).returncode != 0


# =========================================================================== #
# prodigal_genes: training-length rule (awk of the rule) vs real Prodigal
# =========================================================================== #
@suite.test
def test_prodigal_training_length_matches_prodigal():
    prog = shell_awk("prodigal_genes", "bp=$(LC_ALL=C awk '")
    d = tdir("trainlen")
    cases = {  # file -> Prodigal's slen (letters 'A'..'z' + 12 bp per record when > 1 record)
        "one.fna": (">a desc\nACGTNacgtn\r\nAC-GT*\n\n", 14),
        "two.fna": (">a\nACGT\n>b\nACGTACGT\n", 12 + 24),
        "hdr.fna": (">a\n", 0),
    }
    for name, (text, want) in cases.items():
        (d / name).write_text(text)
        got = int(run(["awk", prog, str(d / name)]).stdout.strip())
        assert got == want, (name, got, want)
    if not have("prodigal"):
        print("    (real Prodigal not on PATH: boundary run skipped)")
        return
    # two records; slen = letters + 24 -> 20000 is accepted, 19999 exits 10 (main.c).
    for letters, want_rc in ((19976, 0), (19975, 10)):
        f = d / f"b{letters}.fna"
        s = rand_dna(letters, letters)
        f.write_text(f">c1\n{s[:9000]}\n>c2\n{s[9000:]}\n")
        bp = int(run(["awk", prog, str(f)]).stdout.strip())
        p = run(["prodigal", "-i", str(f), "-o", str(d / "o.gff"), "-f", "gff", "-p", "single",
                 "-g", "11", "-q"], check=False)
        assert p.returncode == want_rc, (letters, p.returncode, p.stderr[-300:])
        assert (bp >= 20000) == (want_rc == 0), (letters, bp)


# =========================================================================== #
# run_bbmap_stats.py
# =========================================================================== #
REAL_BBMAP_39_06 = (  # verbatim stdout of statswrapper.sh in=gA.fna,gB.fna format=6 n_=t (BBMap 39.06)
    "#n_scaffolds\tn_contigs\tscaf_bp\tcontig_bp\tgap_pct\tscaf_N50\tscaf_L50\tctg_N50\tctg_L50\t"
    "scaf_N90\tscaf_L90\tctg_N90\tctg_L90\tscaf_max\tctg_max\tscaf_n_gt50K\tscaf_pct_gt50K\tgc_avg\t"
    "gc_std\tfilename\n"
    "4\t5\t115020\t115000\t0.017\t2\t30000\t2\t30000\t4\t15020\t4\t10000\t50000\t50000\t0\t0.000\t"
    "0.50070\t0.00216\t./bbtest/gA.fna\n"
    "3\t3\t1900\t1900\t0.000\t2\t600\t2\t600\t3\t500\t3\t500\t800\t800\t0\t0.000\t0.30211\t0.01555\t"
    "./bbtest/gB.fna\n")


@suite.test
def test_parse_statswrapper_real_bbmap_output():
    import run_bbmap_stats as rbs
    cols, rows = rbs.parse_statswrapper(REAL_BBMAP_39_06, "fna")
    assert cols == synth.BBMAP_HEADER[:-1]
    rec = dict(rows)
    # gA = contigs 50/30/20 kb + a 10 kb + N20 + 5 kb scaffold: contig N50 LENGTH is
    # ctg_L50 = 30000 bp, ctg_N50 = 2 is the count; contig_bp excludes the 20 gap Ns.
    assert rec["gA"]["ctg_L50"] == "30000" and rec["gA"]["ctg_N50"] == "2"
    assert rec["gA"]["contig_bp"] == "115000" and rec["gA"]["scaf_bp"] == "115020"
    assert rec["gA"]["n_contigs"] == "5" and rec["gA"]["n_scaffolds"] == "4"
    assert float(rec["gA"]["gc_avg"]) == 0.5007          # a fraction, not a percentage
    try:
        rbs.parse_statswrapper(REAL_BBMAP_39_06.replace("\t./bbtest/gB.fna", ""), "fna")
        raise AssertionError("short line accepted")
    except rbs.StatsError:
        pass


def _bbmap_genomes(d):
    gd = Path(d) / "genomes"
    gd.mkdir(parents=True, exist_ok=True)
    a = rand_dna(12000, 1)
    synth.write_fasta(gd / "G1.fna", ["s1", "s2"], [a[:7000], a[7000:]])
    synth.write_fasta(gd / "G2.fna", ["s1"], [rand_dna(9000, 2) + "N" * 25 + rand_dna(3000, 3)])
    for i in (3, 4, 5):
        synth.write_fasta(gd / f"G{i}.fna", ["s1"], [rand_dna(4000 + 1000 * i, i)])
    return gd


@suite.test
def test_run_bbmap_stats_chunks_and_failures():
    d = tdir("bbmap")
    gd = _bbmap_genomes(d)
    fakes = synth.write_fake_tools(d / "bin", PYTHON)
    out = d / "stats.tsv"
    run_env([PYTHON, str(SCRIPTS_PY / "run_bbmap_stats.py"), "--genomes-dir", str(gd), "--chunk", "2",
             "--out", str(out)], path_prepend=[fakes])
    t = read_tsv(out).set_index("genome_id")
    assert t.index.tolist() == ["G1", "G2", "G3", "G4", "G5"]
    assert t.loc["G2", "n_scaffolds"] == "1" and t.loc["G2", "n_contigs"] == "2"     # split at N25
    assert t.loc["G2", "contig_bp"] == "12000" and t.loc["G2", "scaf_bp"] == "12025"
    assert t.loc["G2", "ctg_L50"] == "9000" and t.loc["G1", "ctg_L50"] == "7000"
    base = [PYTHON, str(SCRIPTS_PY / "run_bbmap_stats.py"), "--genomes-dir", str(gd), "--chunk", "2",
            "--out", str(d / "x.tsv")]
    expect_fail(run_env(base, {"FAKE_BBMAP_FAIL": 1}, [fakes], check=False), r"exited 1")
    expect_fail(run_env(base, {"FAKE_BBMAP_OMIT": "G4"}, [fakes], check=False), r"without stats.*G4")
    expect_fail(run_env(base + ["--statswrapper", "no_such_statswrapper.sh"], check=False), r"not found on PATH")
    empty = d / "empty"
    empty.mkdir()
    run_env([PYTHON, str(SCRIPTS_PY / "run_bbmap_stats.py"), "--genomes-dir", str(empty),
             "--out", str(d / "e.tsv")])
    assert (d / "e.tsv").read_text() == "genome_id\n"


def _real_statswrapper():
    """statswrapper.sh of a real BBMap: $PROM_BBMAP_STATSWRAPPER, else one on the toolchain PATH."""
    exe = (os.environ.get("PROM_BBMAP_STATSWRAPPER", "")
           or shutil.which("statswrapper.sh", path=tool_env()["PATH"]) or "")
    return exe if exe and os.path.exists(exe) else ""


@suite.test
def test_run_bbmap_stats_real_bbmap_if_available():
    exe = _real_statswrapper()
    if not exe or not shutil.which("java", path=tool_env()["PATH"]):
        print("    (real BBMap not available: put statswrapper.sh and java on PATH, "
              "or set PROM_BBMAP_STATSWRAPPER)")
        return
    d = tdir("bbreal")
    gd = _bbmap_genomes(d)
    fakes = synth.write_fake_tools(d / "bin", PYTHON)
    jvm = {"JAVA_TOOL_OPTIONS": "-XX:ActiveProcessorCount=2"}  # at most 2 threads per tool (shared cores)
    run_env([PYTHON, str(SCRIPTS_PY / "run_bbmap_stats.py"), "--genomes-dir", str(gd), "--chunk", "2",
             "--statswrapper", exe, "--out", str(d / "real.tsv")], extra=jvm)
    run_env([PYTHON, str(SCRIPTS_PY / "run_bbmap_stats.py"), "--genomes-dir", str(gd), "--chunk", "2",
             "--out", str(d / "fake.tsv")], path_prepend=[fakes])
    real, fake = read_tsv(d / "real.tsv"), read_tsv(d / "fake.tsv")
    assert list(real.columns) == list(fake.columns) == ["genome_id"] + synth.BBMAP_HEADER[:-1]
    for c in ("genome_id", "n_scaffolds", "n_contigs", "scaf_bp", "contig_bp", "ctg_N50", "ctg_L50",
              "scaf_L50", "scaf_max", "ctg_max"):
        assert real[c].tolist() == fake[c].tolist(), (c, real[c].tolist(), fake[c].tolist())
    assert all(abs(float(a) - float(b)) < 1e-4 for a, b in zip(real["gc_avg"], fake["gc_avg"]))
    print("    (real BBMap 39.06 statswrapper.sh agrees with the stand-in)")


# =========================================================================== #
# merge_annotations.py
# =========================================================================== #
SMITHII = ("d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;"
           "f__Methanobacteriaceae;g__Methanobrevibacter;s__Methanobrevibacter smithii")
NOVEL = ("d__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;"
         "f__Methanobacteriaceae;g__Methanobrevibacter;s__")
ECOLI = ("d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;"
         "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli")
BB_HDR = ["genome_id"] + synth.BBMAP_HEADER[:-1]


def _bb_row(gid, n_contigs, contig_bp, ctg_l50, gc, scaf_bp=None, ctg_n50=2):
    scaf_bp = scaf_bp or contig_bp
    r = dict(zip(BB_HDR, [gid, 1, n_contigs, scaf_bp, contig_bp, "0.000", 1, scaf_bp, ctg_n50, ctg_l50,
                          1, scaf_bp, 1, ctg_l50, scaf_bp, ctg_l50, 0, "0.000", gc, "0.01"]))
    return [r[c] for c in BB_HDR]


def _cm_row(gid, comp, cont, n50=1000):
    return [gid, comp, cont, "Neural Network (Specific Model)", 11, 0.9, n50, 300, 100000, 0.5, 90, 3,
            5000, "None"]


def _merge_fixture(d):
    """7 genomes in 2 batches with every GTDB-Tk/CheckM2/BBMap situation."""
    ids = ["M1", "M2", "M3", "M4", "M5", "M6", "M7"]
    meta = d / "meta.tsv"
    with open(meta, "w") as fh:
        fh.write("genome_id\tfilename\tsource\tgroup_gut\n# comment\n\n")
        for g in ids:
            fh.write(f"{g}\t{g}.fasta\tNCBI\t1\n")
    man = write_tsv(d / "manifest.tsv", ["genome_id", "batch", "source_path", "staged_path", "gtdbtk_id"],
                    [[g, "0000" if g <= "M3" else "0001", f"/data/{g}.fasta", f"/b/{g}.fna", f"user_{g}"]
                     for g in ids])
    gh = synth.GTDBTK_HEADER + ["marker_set"]
    na = ["N/A"] * (len(gh) - 3)
    g0 = write_tsv(d / "b0" / "gtdbtk_summary.tsv", gh, [
        ["M1", SMITHII] + na + ["ar53"], ["M2", "Unclassified Archaea"] + na + ["ar53"],
        ["M3", NOVEL] + na + ["ar53"]])
    g1 = write_tsv(d / "b1" / "gtdbtk_summary.tsv", gh, [
        ["M4", ECOLI] + na + ["bac120"], ["M5", SMITHII.replace("smithii", "oralis")] + na + ["ar53"],
        ["M7", "Unclassified"] + na + ["bac120"]])                       # M6 absent -> missing
    c0 = write_tsv(d / "b0" / "checkm2_quality.tsv", synth.CHECKM2_HEADER,
                   [_cm_row("M1", "95.12", "1.03", 30000), _cm_row("M2", "30.00", "0.00"),
                    _cm_row("M3", "64.15", "2.83")])
    c1 = write_tsv(d / "b1" / "checkm2_quality.tsv", synth.CHECKM2_HEADER,
                   [_cm_row("M4", "99.00", "0.00"), _cm_row("M5", "75.00", "5.00"),
                    _cm_row("M7", "10.00", "0.00")])                     # M6 skipped by CheckM2
    b0 = write_tsv(d / "b0" / "bbmap_stats.tsv", BB_HDR,
                   [_bb_row("M1", 5, 115000, 30000, "0.50070", scaf_bp=115020), _bb_row("M2", 3, 9000, 4000, "0.4"),
                    _bb_row("M3", 10, 2000000, 250000, "0.31")])
    b1 = write_tsv(d / "b1" / "bbmap_stats.tsv", BB_HDR,
                   [_bb_row(g, 2, 1500000, 900000, "0.3") for g in ("M4", "M5", "M6", "M7")])
    return meta, man, (g0, g1), (c0, c1), (b0, b1)


def _merge_cmd(meta, man, g, c, b, out, extra=()):
    return ["merge_annotations.py", "--metadata", meta, "--manifest", man, "--gtdbtk", *g, "--checkm2", *c,
            "--bbmap", *b, "--out", out, *extra]


@suite.test
def test_merge_annotations_values():
    d = tdir("merge")
    meta, man, g, c, b = _merge_fixture(d)
    out = d / "annotated.tsv"
    py(*_merge_cmd(meta, man, g, c, b, out))
    t = read_tsv(out).set_index("genome_id")
    assert t.index.tolist() == ["M1", "M2", "M3", "M4", "M5", "M6", "M7"]
    m1 = t.loc["M1"]
    assert (m1["gtdbtk_status"], m1["gtdb_domain"], m1["gtdb_genus"], m1["gtdb_species"]) == (
        "classified", "Archaea", "Methanobrevibacter", "Methanobrevibacter smithii")
    assert m1["classification"] == SMITHII and m1["gtdbtk_marker_set"] == "ar53"
    assert float(m1["completeness"]) == 95.12 and float(m1["contamination"]) == 1.03
    assert m1["quality_score"] == "89.97"                         # 95.12 - 5 * 1.03, exactly
    assert float(m1["n50"]) == 30000 and float(m1["genome_size"]) == 115000   # ctg_L50, contig_bp
    assert float(m1["n_contigs"]) == 5 and float(m1["gc"]) == 0.5007
    assert (m1["batch"], m1["source_path"]) == ("0000", "/data/M1.fasta")
    assert t.loc["M2", "gtdbtk_status"] == "unclassified" and t.loc["M2", "gtdb_domain"] == ""
    assert t.loc["M3", "gtdbtk_status"] == "classified" and t.loc["M3", "gtdb_species"] == ""
    assert t.loc["M3", "quality_score"] == "50.0"                 # 64.15 - 5*2.83: float noise removed
    assert t.loc["M5", "quality_score"] == "50.0"
    assert t.loc["M4", "gtdb_domain"] == "Bacteria" and t.loc["M4", "gtdbtk_marker_set"] == "bac120"
    assert t.loc["M6", "gtdbtk_status"] == "missing" and t.loc["M6", "checkm2_status"] == "missing"
    assert t.loc["M6", "completeness"] == "" and t.loc["M6", "bbmap_status"] == "ok"
    assert t.loc["M7", "gtdbtk_status"] == "unclassified" and t.loc["M7", "checkm2_status"] == "ok"
    assert t.loc["M1", "checkm2_contig_n50"] == "30000" and t.loc["M1", "checkm2_model"].startswith("Neural")
    for col in ("group_gut", "filename", "source"):                # input metadata kept verbatim
        assert col in t.columns


@suite.test
def test_merge_annotations_hard_errors():
    d = tdir("mergeerr")
    meta, man, g, c, b = _merge_fixture(d)
    out = d / "x.tsv"
    dup = write_tsv(d / "dup" / "gtdbtk_summary.tsv", synth.GTDBTK_HEADER + ["marker_set"],
                    [["M1", SMITHII] + ["N/A"] * 18 + ["ar53"]])
    expect_fail(py(*_merge_cmd(meta, man, (*g, dup), c, b, out), check=False), r"more than once")
    stem = write_tsv(d / "stem" / "checkm2_quality.tsv", synth.CHECKM2_HEADER, [_cm_row("M1.fna", "95", "1")])
    expect_fail(py(*_merge_cmd(meta, man, g, (stem,), b, out), check=False), r"not genomes of the input")
    nocol = write_tsv(d / "nocol" / "bbmap_stats.tsv", ["genome_id", "n_contigs", "contig_bp", "gc_avg"],
                      [["M1", 5, 100, 0.5]])
    expect_fail(py(*_merge_cmd(meta, man, g, c, (nocol,), out), check=False), r"BBMap column 'ctg_L50'")
    kb = write_tsv(d / "kb" / "bbmap_stats.tsv", BB_HDR, [_bb_row("M1", 5, 115000, "30 KB", "0.5")])
    expect_fail(py(*_merge_cmd(meta, man, g, c, (kb,), out), check=False), r"non-numeric")
    man2 = write_tsv(d / "man2.tsv", ["genome_id", "batch", "source_path", "staged_path", "gtdbtk_id"],
                     [["M1", "0000", "/x", "/y", "user_M1"]])
    expect_fail(py(*_merge_cmd(meta, man2, g, c, b, out), check=False), r"does not match")
    swapped = ("--n50-col", "scaf_bp", "--length-col", "ctg_L50")
    expect_fail(py(*_merge_cmd(meta, man, g, c, b, out, swapped), check=False), r"n50 > genome_size")


@suite.test
def test_merge_annotations_empty_input():
    d = tdir("mergeempty")
    meta = write_tsv(d / "meta.tsv", ["genome_id", "filename", "group_gut"], [])
    man = write_tsv(d / "man.tsv", ["genome_id", "batch", "source_path", "staged_path", "gtdbtk_id"], [])
    py("merge_annotations.py", "--metadata", meta, "--manifest", man, "--gtdbtk", "--checkm2", "--bbmap",
       "--out", d / "a.tsv")
    t = read_tsv(d / "a.tsv")
    assert t.empty
    for col in ("genome_id", "classification", "gtdb_species", "completeness", "quality_score", "n50",
                "gtdbtk_status", "checkm2_status", "bbmap_status"):
        assert col in t.columns, col


# =========================================================================== #
# filter_collapse.py
# =========================================================================== #
ANN_COLS = ["genome_id", "filename", "source", "group_gut", "niche_rumen", "batch", "source_path",
            "gtdbtk_status", "classification", "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
            "gtdb_family", "gtdb_genus", "gtdb_species", "completeness", "contamination",
            "quality_score", "gc", "genome_size", "n_contigs", "n50"]


def _ann(gid, species, comp, cont, n_contigs=10, n50=50000, gut=0, rumen=0, domain="Archaea",
         status="classified", src_dir=None, size=2000000):
    cls = f"d__{domain};p__P;c__C;o__O;f__F;g__G;s__{species}" if status == "classified" else "Unclassified"
    src = str(Path(src_dir) / f"{gid}.fasta") if src_dir else ""
    qs = "" if comp == "" or cont == "" else round(float(comp) - 5 * float(cont), 6)
    return [gid, f"{gid}.fasta", "NCBI", gut, rumen, "0000", src, status, cls,
            domain if status == "classified" else "", "P", "C", "O", "F", "G",
            species if status == "classified" else "", comp, cont, qs, 0.45, size, n_contigs, n50]


def _run_fc(d, rows, *extra, check=True, stage=None, thresholds=None):
    ann = write_tsv(d / "annotated.tsv", ANN_COLS, rows)
    th = {"--min-completeness": 90, "--max-contamination": 5, "--max-contigs": 1000, "--min-n50": 5000,
          "--min-quality-score": 50}
    th.update(thresholds or {})
    args = ["filter_collapse.py", "--annotated", ann, "--stage-dir", stage or (d / "stage"),
            "--out-species", d / "species.tsv", "--out-reps-list", d / "reps.txt",
            "--out-genome-meta", d / "gmeta.tsv", "--out-taxonomy", d / "tax.tsv",
            "--out-staged", d / "stage" / ".staged", "--report", d / "report.tsv"]
    for k, v in th.items():
        args += [k, v]
    return py(*args, *extra, check=check)


def _fasta_dir(d, ids):
    src = d / "src"
    src.mkdir(parents=True, exist_ok=True)
    for i, g in enumerate(ids):
        synth.write_fasta(src / f"{g}.fasta", ["c1"], [rand_dna(300, i)])
    return src


@suite.test
def test_filter_collapse_thresholds_are_exact():
    d = tdir("fcth")
    rows_spec = [  # gid, comp, cont, n_contigs, n50, expected_pass, expected reason
        ("P_comp_eq", "90.00", "0.00", 10, 50000, True, ""),
        ("F_comp", "89.99", "0.00", 10, 50000, False, "completeness<90"),
        ("P_cont_eq", "100.00", "5.00", 10, 50000, True, ""),
        ("F_cont", "100.00", "5.01", 10, 50000, False, "contamination>5"),
        ("P_contigs", "95.00", "1.00", 999, 50000, True, ""),
        ("F_contigs", "95.00", "1.00", 1000, 50000, False, "n_contigs>=1000"),
        ("P_n50", "95.00", "1.00", 10, 5001, True, ""),
        ("F_n50", "95.00", "1.00", 10, 5000, False, "n50<=5000"),
        ("F_comp_na", "", "1.00", 10, 50000, False, "completeness_missing;quality_score_missing"),
        ("F_n50_na", "95.00", "1.00", 10, "", False, "n50_missing"),
    ]
    ids = [r[0] for r in rows_spec]
    src = _fasta_dir(d, ids)
    rows = [_ann(g, f"S {g}", c, k, nc, n5, src_dir=src) for g, c, k, nc, n5, _, _ in rows_spec]
    _run_fc(d, rows)
    rep = read_tsv(d / "report.tsv").set_index("genome_id")
    for g, *_x, want, reason in rows_spec:
        assert rep.loc[g, "passing"] == str(want), (g, rep.loc[g].to_dict())
        assert rep.loc[g, "reasons"] == reason, (g, rep.loc[g, "reasons"])
        assert rep.loc[g, "status"] == ("representative" if want else "failed_quality")
    # quality_score: exactly 50 fails `> 50` (64.15 - 5*2.83 is 50.00000000000001 in raw floats)
    d2 = tdir("fcqs")
    src2 = _fasta_dir(d2, ["Q_eq", "Q_float", "Q_pass"])
    rows2 = [_ann("Q_eq", "S1", "75.00", "5.00", src_dir=src2), _ann("Q_float", "S2", "64.15", "2.83", src_dir=src2),
             _ann("Q_pass", "S3", "75.01", "5.00", src_dir=src2)]
    rows2[1][ANN_COLS.index("quality_score")] = 64.15 - 5 * 2.83          # a raw-float table value
    _run_fc(d2, rows2, thresholds={"--min-completeness": 60})
    rep2 = read_tsv(d2 / "report.tsv").set_index("genome_id")
    assert rep2.loc["Q_eq", "passing"] == "False" and rep2.loc["Q_eq", "reasons"] == "quality_score<=50"
    assert rep2.loc["Q_float", "passing"] == "False" and rep2.loc["Q_float", "quality_score"] == "50.0"
    assert rep2.loc["Q_pass", "passing"] == "True"


@suite.test
def test_filter_collapse_representative_tie_breaks():
    d = tdir("fcties")
    ids = ["s1_a", "s1_b", "s1_c1", "s1_c0", "s2_x"]
    src = _fasta_dir(d, ids)
    rows = [  # same quality_score/completeness/contamination within S1
        _ann("s1_a", "S1", "95.00", "1.00", n_contigs=50, n50=10000, src_dir=src),
        _ann("s1_b", "S1", "95.00", "1.00", n_contigs=60, n50=20000, src_dir=src),
        _ann("s1_c1", "S1", "95.00", "1.00", n_contigs=40, n50=20000, src_dir=src),
        _ann("s1_c0", "S1", "95.00", "1.00", n_contigs=40, n50=20000, src_dir=src),
        _ann("s2_x", "S2", "99.00", "0.50", src_dir=src),
    ]
    for order in (rows, rows[::-1], [rows[2], rows[0], rows[4], rows[3], rows[1]]):
        _run_fc(d, order)
        sp = read_tsv(d / "species.tsv").set_index("gtdb_species")
        assert sp.loc["S1", "genome_id"] == "s1_c0", sp.loc["S1", "genome_id"]   # n50, contigs, then id
        assert sp.loc["S1", "n_genomes_in_species"] == "4" and sp.loc["S1", "n_passing_in_species"] == "4"
        assert (d / "reps.txt").read_text() == "s1_c0\ns2_x\n"


@suite.test
def test_filter_collapse_niche_union_modes():
    d = tdir("fcniche")
    ids = ["X", "Y", "Z", "W", "P1", "P2", "Q1", "Q2", "R1", "R2"]
    src = _fasta_dir(d, ids)
    rows = [_ann("X", "S1", "95.00", "1.00", gut=1, rumen=0, src_dir=src),
            _ann("Y", "S1", "80.00", "1.00", gut=0, rumen=1, src_dir=src),        # fails completeness
            _ann("Z", "S2", "95.00", "1.00", gut=0, rumen="", src_dir=src),       # blank -> unknown
            _ann("W", "", "95.00", "1.00", gut=0, rumen=1, src_dir=src),          # novel species
            _ann("P1", "S3", "96.00", "1.00", gut="", rumen="", src_dir=src),     # blank + 0 -> 0
            _ann("P2", "S3", "95.00", "1.00", gut=0, rumen="NA", src_dir=src),    # blank + blank -> unknown
            _ann("Q1", "S4", "96.00", "1.00", gut="", rumen="", src_dir=src),     # blank + 1 -> 1
            _ann("Q2", "S4", "95.00", "1.00", gut=1, rumen="", src_dir=src),
            _ann("R1", "S5", "96.00", "1.00", gut="", rumen=0, src_dir=src),      # passing: blank
            _ann("R2", "S5", "70.00", "1.00", gut=1, rumen=0, src_dir=src)]       # failing: 1
    _run_fc(d, rows, "--niche-aggregation", "all")
    sp = read_tsv(d / "species.tsv").set_index("genome_id")
    assert sp.index.tolist() == ["X", "Z", "P1", "Q1", "R1"], sp.index.tolist()
    assert (sp.loc["X", "group_gut"], sp.loc["X", "niche_rumen"]) == ("1", "1")
    assert (sp.loc["X", "n_genomes_in_species"], sp.loc["X", "n_passing_in_species"]) == ("2", "1")
    # unknown is not absence (SCHEMA: species-level {1, 0, NA}; blank = NA, never 0)
    assert (sp.loc["Z", "group_gut"], sp.loc["Z", "niche_rumen"]) == ("0", "")
    assert (sp.loc["P1", "group_gut"], sp.loc["P1", "niche_rumen"]) == ("0", "")
    assert (sp.loc["Q1", "group_gut"], sp.loc["Q1", "niche_rumen"]) == ("1", "")
    assert (sp.loc["R1", "group_gut"], sp.loc["R1", "niche_rumen"]) == ("1", "0")
    _run_fc(d, rows, "--niche-aggregation", "passing")
    sp = read_tsv(d / "species.tsv").set_index("genome_id")
    assert (sp.loc["X", "group_gut"], sp.loc["X", "niche_rumen"]) == ("1", "0")
    assert (sp.loc["X", "n_genomes_in_species"], sp.loc["X", "n_passing_in_species"]) == ("2", "1")
    assert (sp.loc["R1", "group_gut"], sp.loc["R1", "niche_rumen"]) == ("", "0")   # only R1 passes
    bad = [r[:] for r in rows]
    bad[0][ANN_COLS.index("group_gut")] = 2
    expect_fail(_run_fc(d, bad, check=False), r"group_gut.*0/1")


@suite.test
def test_filter_collapse_domain_unclassified_tiers_and_files():
    d = tdir("fcdom")
    ids = ["A1", "B1", "N1", "U1", "M1", "L1", "C1"]
    src = _fasta_dir(d, ids)
    rows = [_ann("A1", "Sa", "95.00", "1.00", src_dir=src),
            _ann("B1", "Sb", "99.00", "0.00", domain="Bacteria", src_dir=src),
            _ann("N1", "", "96.00", "0.00", src_dir=src),                    # novel species (s__)
            _ann("U1", "", "", "", status="unclassified", src_dir=src),
            _ann("M1", "", "", "", status="missing", src_dir=src),
            _ann("L1", "Sl", "60.00", "2.00", src_dir=src),                  # tier B
            _ann("C1", "Sc", "45.00", "1.00", src_dir=src)]                  # tier C, not pass_primary
    relaxed = {"--min-completeness": 40, "--max-contamination": 10, "--min-quality-score": 0}
    _run_fc(d, rows, "--keep-domain", "d__Archaea", thresholds=relaxed)
    rep = read_tsv(d / "report.tsv").set_index("genome_id")
    assert rep.loc["B1", "status"] == "excluded_domain" and rep.loc["B1", "reasons"] == "domain=Bacteria"
    assert rep.loc["U1", "reasons"].startswith("gtdbtk_unclassified")
    assert rep.loc["M1", "reasons"].startswith("gtdbtk_missing")
    assert rep.loc["N1", "status"] == "no_species" and rep.loc["N1", "reasons"] == "no_gtdb_species"
    assert (d / "reps.txt").read_text() == "A1\nC1\nL1\n"                 # sorted by species key
    gm = read_tsv(d / "gmeta.tsv")
    assert "\t".join(gm.columns) == ("genome_id\tsource\tcompleteness\tcontamination\tgenome_size\t"
                                     "n_contigs\tn50\tgc\tquality_tier\tpass_primary")
    gm = gm.set_index("genome_id")
    assert [gm.loc[g, "quality_tier"] for g in ("A1", "L1", "C1")] == ["A", "B", "C"]
    assert [gm.loc[g, "pass_primary"] for g in ("A1", "L1", "C1")] == ["True", "True", "False"]
    assert gm.loc["A1", "genome_size"] == "2000000" and gm.loc["A1", "n50"] == "50000"   # no '.0'
    tax = read_tsv(d / "tax.tsv")
    assert list(tax.columns) == ["user_genome", "classification"] and tax["user_genome"].tolist() == ["A1", "C1", "L1"]
    sp = read_tsv(d / "species.tsv")
    for c in ("genome_id", "filename", "source", "classification", "gtdb_domain", "gtdb_species",
              "assembly_path", "completeness", "contamination", "gc", "genome_size", "n_contigs", "n50",
              "quality_score", "n_genomes_in_species", "n_passing_in_species", "group_gut", "niche_rumen"):
        assert c in sp.columns, c
    sp = sp.set_index("genome_id")
    assert sp.loc["A1", "assembly_path"] == str(d / "stage" / "A1.fna")
    assert os.readlink(d / "stage" / "A1.fna") == str((src / "A1.fasta").resolve())
    assert (d / "stage" / ".staged").read_text() == "genome_id\nA1\nC1\nL1\n"
    # same result for 'Archaea'; singleton mode keeps the novel-species genome
    _run_fc(d, rows, "--keep-domain", "Archaea", thresholds=relaxed)
    assert (d / "reps.txt").read_text() == "A1\nC1\nL1\n"
    _run_fc(d, rows, "--unclassified", "singleton", thresholds=relaxed)
    assert (d / "reps.txt").read_text() == "A1\nC1\nL1\nN1\n"


@suite.test
def test_filter_collapse_staging_is_safe():
    d = tdir("fcstage")
    src = _fasta_dir(d, ["R1", "R2"])
    stage = d / "stage"
    stage.mkdir()
    os.symlink(str(src / "R2.fasta"), stage / "other.fna")                # a link the script did not make
    _run_fc(d, [_ann("R1", "S1", "95.00", "1.00", src_dir=src)], stage=stage)
    assert (stage / "R1.fna").is_symlink()
    _run_fc(d, [_ann("R2", "S2", "95.00", "1.00", src_dir=src)], stage=stage)
    assert not (stage / "R1.fna").exists() and (stage / "R2.fna").is_symlink()   # stale rep link removed
    assert (stage / "other.fna").is_symlink()                                  # untouched
    (stage / "R2.fna").unlink()
    (stage / "R2.fna").write_text(">real file\nACGT\n")
    p = _run_fc(d, [_ann("R2", "S2", "95.00", "1.00", src_dir=src)], stage=stage, check=False)
    expect_fail(p, r"refusing to overwrite")
    assert (stage / "R2.fna").read_text() == ">real file\nACGT\n"
    p = _run_fc(d, [_ann("R3", "S3", "95.00", "1.00", src_dir=src)], stage=stage, check=False)
    expect_fail(p, r"not found")                                                # missing source FASTA


@suite.test
def test_filter_collapse_empty_input():
    d = tdir("fcempty")
    _run_fc(d, [])
    assert (d / "reps.txt").read_text() == ""
    assert (d / "stage" / ".staged").read_text() == "genome_id\n"
    assert read_tsv(d / "tax.tsv").columns.tolist() == ["user_genome", "classification"]
    sp = read_tsv(d / "species.tsv")
    assert sp.empty and "assembly_path" in sp.columns and "group_gut" in sp.columns
    assert read_tsv(d / "report.tsv").empty and read_tsv(d / "gmeta.tsv").empty


# =========================================================================== #
# rename_prodigal.py (+ real Prodigal)
# =========================================================================== #
PLANT_NAMES = ["ctgA some description", "ctgEmpty", "ctgB", "ctgC len=short"]


def _expected_gene_rows(genome, recs, genes):
    names = [n.split()[0] for n in PLANT_NAMES]
    per = {}
    rows = []
    for ci, s, e, st, part, ncod in genes:
        per[ci] = per.get(ci, 0) + 1
        rows.append([f"{genome}_{ci}_{per[ci]}", names[ci - 1], str(ci), str(per[ci]), str(s), str(e), st, part,
                     str(len(recs[ci - 1])), str(ncod), "11"])
    return rows


def _read_fa(path):
    import utils
    return list(utils.read_fasta(str(path)))


@suite.test
def test_prodigal_real_run_contract_v2_ids_and_genes_tsv():
    if not have("prodigal"):
        raise AssertionError("prodigal not on PATH (PROM_TOOLS_BIN) -- the real-tool run is required")
    import utils
    d = tdir("prod")
    recs, genes = synth.planted_genome()
    synth.write_fasta(d / "G1.fna", PLANT_NAMES, recs)
    run(["prodigal", "-i", str(d / "G1.fna"), "-a", str(d / "raw.faa"), "-d", str(d / "raw.fna"), "-f", "gff",
         "-o", str(d / "raw.gff"), "-p", "single", "-g", "11", "-q"])
    py("rename_prodigal.py", "--genome", "G1", "--in-faa", d / "raw.faa", "--in-fna", d / "raw.fna",
       "--in-gff", d / "raw.gff", "--faa", d / "G1.faa", "--fna", d / "G1.fna.cds", "--gff", d / "G1.gff",
       "--genes-tsv", d / "G1.genes.tsv")
    t = read_tsv(d / "G1.genes.tsv")
    assert list(t.columns) == ["protein_id", "contig", "contig_index", "gene_index", "start", "end", "strand",
                               "partial", "contig_length", "length_aa", "transl_table"]
    assert t.values.tolist() == _expected_gene_rows("G1", recs, genes)
    assert not t["protein_id"].str.contains("_2_").any()             # gene-less contig 2 keeps its index
    assert t["partial"].tolist().count("01") == 1                       # the right-edge gene
    faa, fna = _read_fa(d / "G1.faa"), _read_fa(d / "G1.fna.cds")
    assert [i for i, _ in faa] == [i for i, _ in fna] == t["protein_id"].tolist()
    for (pid, prot), (_, cds) in zip(faa, fna):
        assert synth.translate(cds).lstrip("M") == prot.lstrip("M"), pid
        assert utils.split_gene_id(pid)[0] == "G1"
        assert utils.make_gene_id("G1", *utils.split_gene_id(pid)[1:]) == pid
    gff_ids = re.findall(r"\tID=([^;]+);", (d / "G1.gff").read_text())
    assert gff_ids == t["protein_id"].tolist()
    assert (d / "G1.gff").read_text().count("# Sequence Data:") == 4


def _tiny_prodigal(d, n_genes=2):
    faa = ["".join(f">c1_{i} # {1 + 300 * (i - 1)} # {300 * i} # 1 # ID=1_{i};partial=00;start_type=ATG\n"
                   f"{'M' + 'A' * 98 + '*'}\n") for i in range(1, n_genes + 1)]
    fna = ["".join(f">c1_{i} # {1 + 300 * (i - 1)} # {300 * i} # 1 # ID=1_{i};partial=00;start_type=ATG\n"
                   f"{'ATG' + 'GCT' * 98 + 'TAA'}\n") for i in range(1, n_genes + 1)]
    gff = ['##gff-version  3\n# Sequence Data: seqnum=1;seqlen=900;seqhdr="c1 x;y=z"\n'
           '# Model Data: version=Prodigal.v2.6.3;run_type=Single;transl_table=11;uses_sd=1\n']
    gff += [f"c1\tProdigal_v2.6.3\tCDS\t{1 + 300 * (i - 1)}\t{300 * i}\t9.9\t+\t0\tID=1_{i};partial=00;\n"
            for i in range(1, n_genes + 1)]
    (d / "r.faa").write_text("".join(faa))
    (d / "r.fna").write_text("".join(fna))
    (d / "r.gff").write_text("".join(gff))


def _rename(d, genome="T1"):
    return py("rename_prodigal.py", "--genome", genome, "--in-faa", d / "r.faa", "--in-fna", d / "r.fna",
              "--in-gff", d / "r.gff", "--faa", d / "o.faa", "--fna", d / "o.fna", "--gff", d / "o.gff",
              "--genes-tsv", d / "o.tsv", check=False)


@suite.test
def test_rename_prodigal_failure_modes():
    d = tdir("renfail")
    _tiny_prodigal(d)
    assert _rename(d).returncode == 0
    t = read_tsv(d / "o.tsv")
    assert t["protein_id"].tolist() == ["T1_1_1", "T1_1_2"] and t["contig"].tolist() == ["c1", "c1"]
    assert t["length_aa"].tolist() == ["99", "99"] and t["contig_length"].tolist() == ["900", "900"]
    _tiny_prodigal(d, n_genes=0)                                          # Prodigal found no gene
    expect_fail(_rename(d), r"0 genes")
    _tiny_prodigal(d)
    (d / "r.fna").write_text((d / "r.fna").read_text().split(">c1_2")[0])  # CDS record missing
    expect_fail(_rename(d), r"disagree")
    _tiny_prodigal(d)
    (d / "r.gff").write_text((d / "r.gff").read_text().replace("\t301\t600\t", "\t301\t603\t"))
    expect_fail(_rename(d), r"GFF coordinates")
    _tiny_prodigal(d)
    (d / "r.faa").write_text(">c1_1 not a prodigal header\nMAAA*\n")
    expect_fail(_rename(d), r"not Prodigal format")
    _tiny_prodigal(d)
    expect_fail(_rename(d, genome="bad|id"), r"outside")


# =========================================================================== #
# plots
# =========================================================================== #
def _svg_text(path):
    return Path(path).read_text(errors="replace")


@suite.test
def test_plot_preprocess_values_formats_and_empty():
    d = tdir("pp")
    ann = write_tsv(d / "ann.tsv", ["genome_id", "completeness", "contamination"],
                    [["a", 95, 1], ["b", 80, 2], ["c", 99, 30], ["d", "", ""]])
    sp = write_tsv(d / "sp.tsv", ["genome_id"], [["a"]])
    rep = write_tsv(d / "rep.tsv", ["genome_id", "gtdbtk_status", "passing", "status"],
                    [["a", "classified", "True", "representative"], ["b", "classified", "False", "failed_quality"],
                     ["c", "classified", "False", "failed_quality"], ["d", "missing", "False", "excluded_domain"]])
    py("plot_preprocess.py", "--annotated", ann, "--species", sp, "--report", rep, "--outdir", d / "o")
    svg = _svg_text(d / "o" / "preprocess_quality.svg")
    assert (d / "o" / "preprocess_quality.png").stat().st_size > 1000
    assert "pass (all 5 QC thresholds) (n=1)" in svg and "fail (n=2)" in svg
    assert "Filtering funnel" in svg and "1 genome(s) without CheckM2 estimates" in svg
    py("plot_preprocess.py", "--annotated", ann, "--species", sp, "--report", rep, "--formats", "png",
       "--outdir", d / "o2")
    assert sorted(os.listdir(d / "o2")) == ["preprocess_quality.png"]
    empty = write_tsv(d / "e.tsv", ["genome_id", "completeness"], [])
    py("plot_preprocess.py", "--annotated", empty, "--species", sp, "--outdir", d / "o3")
    assert "no data" in _svg_text(d / "o3" / "preprocess_quality.svg")


@suite.test
def test_plot_qc_units_thresholds_and_empty():
    d = tdir("pqc")
    gm = write_tsv(d / "gm.tsv", ["genome_id", "source", "completeness", "contamination", "genome_size",
                                  "n_contigs", "n50", "gc", "quality_tier", "pass_primary"],
                   [["a", "x", 95, 1, 2500000, 30, 150000, 0.41, "A", True],
                    ["b", "x", 91, 4, 1800000, 120, 30000, 0.62, "A", True],
                    ["c", "x", 60, 2, 900000, 400, 6000, 0.35, "B", True]])
    py("plot_qc.py", "--meta", gm, "--outdir", d / "o", "--completeness-floor", 50, "--contamination-ceiling", 10,
       "--tier-high-completeness", 90, "--tier-high-contamination", 5)
    s1, s2 = (_svg_text(d / "o" / f"{n}.svg") for n in ("qc_completeness_contamination", "qc_assembly_stats"))
    assert "tier A: &gt;= 90% / &lt;= 5%" in s1 and "floor: &gt;= 50% / &lt;= 10%" in s1
    assert "tier A (n=2)" in s1 and "tier B (n=1)" in s1
    assert "GC content (%)" in s2 and "Genome size (Mb)" in s2 and "median 41" in s2   # 0.41 -> 41 %
    assert "median 1.8" in s2                                                            # 1.8 Mb
    hdr = write_tsv(d / "h.tsv", ["genome_id", "completeness", "contamination"], [])
    py("plot_qc.py", "--meta", hdr, "--outdir", d / "e")
    assert "no data" in _svg_text(d / "e" / "qc_assembly_stats.svg")


@suite.test
def test_plot_qc_extra_degenerate_groups():
    d = tdir("pqx")
    m = write_tsv(d / "m.tsv", ["genome_id", "prometheus_present", "genome_size",
                                "n_defense_systems_excl_prometheus", "completeness"],
                  [["a", "True", 2000000, 3, 100.0], ["b", "False", 1500000, 5, 100.0],
                   ["c", "False", 1700000, 5, 100.0], ["e", "False", 2100000, 7, 100.0]])
    py("plot_qc_extra.py", "--master", m, "--outdir", d / "o")      # constant / single-value groups
    svg = _svg_text(d / "o" / "prometheus_covariates.svg")
    assert "present\n(n=1)" in svg or "(n=1)" in svg
    assert "defense systems (excl. Prometheus)" in svg and "genome size (Mb)" in svg
    one = write_tsv(d / "one.tsv", ["genome_id", "prometheus_present", "genome_size"], [["a", "True", 1]])
    py("plot_qc_extra.py", "--master", one, "--outdir", d / "p")
    assert "no data" in _svg_text(d / "p" / "prometheus_covariates.svg")


# =========================================================================== #
# Snakemake: the real rules of preprocess.smk + qc.smk, end to end
# =========================================================================== #
STADTMANAE = SMITHII.replace("g__Methanobrevibacter;s__Methanobrevibacter smithii",
                             "g__Methanosphaera;s__Methanosphaera stadtmanae")
LUMINYENSIS = ("d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;"
               "f__Methanomassiliicoccaceae;g__Methanomassiliicoccus;s__Methanomassiliicoccus luminyensis")
INTEG = {  # gid: (GTDB-Tk classification, CheckM2 completeness, contamination, group_gut, niche_rumen)
    "ARC_A1": (SMITHII, "96.00", "6.00", 1, 0),        # fails contamination; same species as A2
    "ARC_A2": (SMITHII, "97.00", "0.50", 0, 1),        # representative of M. smithii
    "ARC_B1": (STADTMANAE, "92.00", "1.00", 0, 0),     # scaffold with a 15-N gap
    "ARC_C1": (NOVEL, "99.00", "0.00", 1, 1),          # novel species (s__) -> dropped
    "BAC_D1": (ECOLI, "99.00", "0.00", 1, 1),          # Bacteria -> excluded
    "ARC_E1": (LUMINYENSIS, "91.00", "2.00", 1, 0),    # 7.7 kb: Prodigal meta fallback
    "ARC_F1": (SMITHII, "99.00", "0.00", 0, 0),        # absent from GTDB-Tk output
}


def _snakemake_exe():
    exe = os.environ.get("PROM_SNAKEMAKE", "")
    return exe if exe and os.path.exists(exe) else shutil.which("snakemake")


def _mk_project(root):
    """Minimal workflow: common.smk + preprocess.smk + qc.smk with the repo's scripts/config."""
    (root / "workflow" / "rules").mkdir(parents=True)
    for f in ("common.smk", "preprocess.smk", "qc.smk"):
        text = (REPO / "workflow" / "rules" / f).read_text()
        # Snakemake 8 reads the env file and runs `conda info` to record job metadata for
        # rules with a conda: directive, even without --use-conda. The tools come from PATH
        # here (env files are checked by tests/check_env_solve.sh), so drop the directive.
        text = re.sub(r"^[ \t]+conda:[ \t]*\n[ \t]+\"[^\"]+\.yaml\"[ \t]*\n", "", text, flags=re.M)
        text = re.sub(r"^[ \t]+conda:[ \t]*\"[^\"]+\.yaml\"[ \t]*\n", "", text, flags=re.M)
        assert "conda:" not in text, f
        (root / "workflow" / "rules" / f).write_text(text)
    shutil.copytree(REPO / "workflow" / "schemas", root / "workflow" / "schemas")
    shutil.copytree(REPO / "workflow" / "scripts" / "py", root / "workflow" / "scripts" / "py",
                    ignore=shutil.ignore_patterns("__pycache__"))
    (root / "config").mkdir()
    shutil.copy(REPO / "config" / "config.yaml", root / "config" / "config.yaml")
    (root / "workflow" / "Snakefile").write_text(
        'include: "rules/common.smk"\ninclude: "rules/preprocess.smk"\ninclude: "rules/qc.smk"\n\n'
        "rule preprocess:\n    input:\n        preprocess_targets()\n")
    gdir = root / "input" / "genomes"
    gdir.mkdir(parents=True)
    planted = {}
    for i, gid in enumerate(sorted(INTEG)):
        layout = ((7, False),) if gid == "ARC_E1" else ((14, False), (0, False), (6, False), (1, True))
        recs, genes = synth.planted_genome(seed=100 + i, layout=layout)
        if gid == "ARC_B1":
            recs[1] = recs[1][:300] + "N" * 15 + recs[1][300:]
        names = PLANT_NAMES[:len(recs)]
        synth.write_fasta(gdir / f"{gid}.fasta", names, recs)
        planted[gid] = (recs, genes)
    with open(root / "input" / "meta.tsv", "w") as fh:
        fh.write("genome_id\tfilename\tsource\tgroup_gut\tniche_rumen\n# commented row\n\n")
        for gid, (_c, _a, _b, gut, rumen) in INTEG.items():
            fh.write(f"{gid}\t{gid}.fasta\tNCBI\t{gut}\t{rumen}\n")
    (root / "db" / "gtdbtk").mkdir(parents=True)
    (root / "db" / "checkm2.dmnd").write_bytes(b"\0" * 64)
    with open(root / "fake_tax.tsv", "w") as fh:
        fh.writelines(f"{g}\t{v[0]}\n" for g, v in INTEG.items())
    with open(root / "fake_checkm2.tsv", "w") as fh:
        fh.writelines(f"{g} {v[1]} {v[2]}\n" for g, v in INTEG.items())
    synth.write_fake_tools(root / "fakebin", PYTHON)
    synth.write_fake_gtdbtk_module(root / "fakepy")
    _write_cfg(root, "cfg.yaml")
    return planted


def _write_cfg(root, name, results="results", logs="logs", checkm2_db="db/checkm2.dmnd", gtdb="db/gtdbtk"):
    (root / name).write_text(
        f'samples: "{results}/preprocess/species_metadata.tsv"\n'
        f'metadata: "{results}/preprocess/species_metadata.tsv"\n'
        f'genomes_dir: "{results}/preprocess/genomes"\ngenome_extension: "fna"\n'
        f'results_dir: "{results}"\nlogs_dir: "{logs}"\n'
        'input: {metadata: "input/meta.tsv", genomes_dir: "input/genomes", genome_extension: "fasta"}\n'
        "preprocess: {gtdbtk_batch_size: 2}\n"
        f'db: {{checkm2: "{checkm2_db}", gtdbtk_data: "{gtdb}"}}\n'
        "threads: {small: 1, medium: 2, large: 2}\n")


def _smk(root, targets, cfg="cfg.yaml", extra_env=None, check=True):
    home = root / "home"
    home.mkdir(exist_ok=True)
    env_extra = {"HOME": home, "XDG_CACHE_HOME": home / ".cache", "MPLCONFIGDIR": home / ".mpl",
                 "PYTHONPATH": root / "fakepy", "FAKE_GTDBTK_TAX": root / "fake_tax.tsv",
                 "FAKE_CHECKM2_TABLE": root / "fake_checkm2.tsv", "FAKE_GTDBTK_DROP": "ARC_F1"}
    env_extra.update(extra_env or {})
    cmd = [_snakemake_exe(), *targets, "--cores", "2", "--nolock", "-s", "workflow/Snakefile",
           "--configfile", cfg]
    return run_env(cmd, env_extra, [root / "fakebin", Path(PYTHON).parent], cwd=root, check=check, timeout=1200)


@suite.test
def test_snakemake_phase1_phase2_end_to_end():
    if not _snakemake_exe():
        raise AssertionError("snakemake not found (set PROM_SNAKEMAKE) -- the rule-level test is required")
    if not have("prodigal"):
        raise AssertionError("prodigal not on PATH (PROM_TOOLS_BIN)")
    root = tdir("smk")
    planted = _mk_project(root)
    _smk(root, ["preprocess"])
    R = root / "results" / "preprocess"
    batches = sorted(p.name for p in (R / "batches").iterdir() if p.name.startswith("batch_"))
    assert batches == ["batch_0000", "batch_0001", "batch_0002", "batch_0003"]    # ceil(7 / 2)
    seen = []
    for b in batches:
        s = read_tsv(R / "batches" / b / "gtdbtk_summary.tsv")
        assert list(s.columns) == synth.GTDBTK_HEADER + ["marker_set"]
        assert not s["user_genome"].duplicated().any() and not s["user_genome"].str.startswith("user_").any()
        seen += s["user_genome"].tolist()
    assert sorted(seen) == sorted(g for g in INTEG if g != "ARC_F1")              # no duplicated rows
    ann = read_tsv(R / "annotated_metadata.tsv").set_index("genome_id")
    assert ann.loc["ARC_F1", "gtdbtk_status"] == "missing"
    assert ann.loc["BAC_D1", "gtdb_domain"] == "Bacteria" and ann.loc["BAC_D1", "gtdbtk_marker_set"] == "bac120"
    b1 = planted["ARC_B1"][0]
    assert float(ann.loc["ARC_B1", "genome_size"]) == sum(len(s) for s in b1) - 15   # contig_bp: gap Ns out
    assert float(ann.loc["ARC_B1", "n_contigs"]) == len(b1) + 1                       # split at the N run
    assert float(ann.loc["ARC_A2", "n50"]) == max(len(s) for s in planted["ARC_A2"][0])
    assert ann.loc["ARC_A1", "quality_score"] == "66.0" and ann.loc["ARC_A2", "quality_score"] == "94.5"
    sp = read_tsv(R / "species_metadata.tsv").set_index("genome_id")
    assert sp.index.tolist() == ["ARC_A2", "ARC_E1", "ARC_B1"]                       # sorted by species
    assert (sp.loc["ARC_A2", "group_gut"], sp.loc["ARC_A2", "niche_rumen"]) == ("1", "1")   # union, 'all'
    assert (sp.loc["ARC_A2", "n_genomes_in_species"], sp.loc["ARC_A2", "n_passing_in_species"]) == ("2", "1")
    assert (R / "representatives.txt").read_text() == "ARC_A2\nARC_E1\nARC_B1\n"
    for g in sp.index:
        link = root / sp.loc[g, "assembly_path"]
        assert link.is_symlink() and link.resolve() == (root / "input" / "genomes" / f"{g}.fasta").resolve()
    rep = read_tsv(R / "filter_report.tsv").set_index("genome_id")
    assert rep["status"].to_dict() == {
        "ARC_A1": "failed_quality", "ARC_A2": "representative", "ARC_B1": "representative",
        "ARC_C1": "no_species", "ARC_E1": "representative", "ARC_F1": "excluded_domain",
        "BAC_D1": "excluded_domain"}
    assert rep.loc["ARC_A1", "reasons"] == "contamination>5"
    assert rep.loc["ARC_F1", "reasons"].startswith("gtdbtk_missing")
    assert (R / "plots" / "preprocess_quality.png").stat().st_size > 1000
    assert "gtdbtk 2.4.0" in (root / "logs" / "preprocess" / "gtdbtk_batch_0000.log").read_text()

    # ---- Phase 2 (the species table now exists): real Prodigal through the rule ----
    G = root / "results" / "genes"
    _smk(root, [f"results/genes/{g}/{g}.genes.tsv" for g in ("ARC_A2", "ARC_E1")]
         + ["results/qc/genome_metadata.tsv", "results/qc/plots/qc_assembly_stats.png"])
    t = read_tsv(G / "ARC_A2" / "ARC_A2.genes.tsv")
    assert t.values.tolist() == _expected_gene_rows("ARC_A2", *planted["ARC_A2"])
    assert [i for i, _ in _read_fa(G / "ARC_A2" / "ARC_A2.faa")] == t["protein_id"].tolist()
    assert [i for i, _ in _read_fa(G / "ARC_A2" / "ARC_A2.fna")] == t["protein_id"].tolist()
    assert not list(G.glob("*/*.raw"))
    assert "mode=single" in (root / "logs" / "prodigal" / "ARC_A2.log").read_text()
    log_e1 = (root / "logs" / "prodigal" / "ARC_E1.log").read_text()
    assert "mode=meta" in log_e1 and "training_length=" in log_e1
    t_e1 = read_tsv(G / "ARC_E1" / "ARC_E1.genes.tsv")
    assert len(t_e1) == 7 and t_e1["protein_id"].str.match(r"^ARC_E1_1_\d+$").all()
    gm = read_tsv(root / "results" / "qc" / "genome_metadata.tsv")
    assert gm["genome_id"].tolist() == ["ARC_A2", "ARC_E1", "ARC_B1"] and set(gm["quality_tier"]) == {"A"}

    # ---- failures: a crashed tool / missing database / crashed Prodigal fails the rule ----
    gfx = root / "results" / "preprocess" / "genomes"
    (gfx / "BOGUS.fna").write_text(">empty\n")
    p = _smk(root, ["results/genes/BOGUS/BOGUS.faa"], check=False)
    assert p.returncode != 0 and not (G / "BOGUS" / "BOGUS.faa").exists()
    assert "no input sequences" in (root / "logs" / "prodigal" / "BOGUS.log").read_text()
    synth.write_fasta(gfx / "NOGENE.fna", ["c1"], [synth.spacer(25000)])
    p = _smk(root, ["results/genes/NOGENE/NOGENE.faa"], check=False)
    assert p.returncode != 0 and not (G / "NOGENE" / "NOGENE.faa").exists()
    assert "0 genes" in (root / "logs" / "prodigal" / "NOGENE.log").read_text()
    scen = [("gt", {"FAKE_GTDBTK_FAIL": 1}, {}, "gtdbtk_summary.tsv", "gtdbtk_batch_0000", "simulated crash"),
            ("cm", {"FAKE_CHECKM2_FAIL": 1}, {}, "checkm2_quality.tsv", "checkm2_batch_0000", "simulated crash"),
            ("cmdb", {}, {"checkm2_db": "db/absent.dmnd"}, "checkm2_quality.tsv", "checkm2_batch_0000",
             "missing or empty"),
            ("gtdb", {}, {"gtdb": "db/absent_dir"}, "gtdbtk_summary.tsv", "gtdbtk_batch_0000", "not a directory"),
            ("bb", {"FAKE_BBMAP_OMIT": "ARC_A1"}, {}, "bbmap_stats.tsv", "bbmap_batch_0000", "without stats")]
    for tag, env_x, cfg_x, out, logname, msg in scen:
        _write_cfg(root, f"cfg_{tag}.yaml", results=f"res_{tag}", logs=f"logs_{tag}", **cfg_x)
        target = f"res_{tag}/preprocess/batches/batch_0000/{out}"
        p = _smk(root, [target], cfg=f"cfg_{tag}.yaml", extra_env=env_x, check=False)
        assert p.returncode != 0, (tag, p.stdout[-800:])
        assert not (root / target).exists(), tag
        assert msg in (root / f"logs_{tag}" / "preprocess" / f"{logname}.log").read_text(), tag


@suite.test
def test_snakemake_empty_input_set():
    """Zero input genomes -> zero batches: Phase 1 still completes with header-only tables."""
    if not _snakemake_exe():
        raise AssertionError("snakemake not found (set PROM_SNAKEMAKE)")
    root = tdir("smkempty")
    _mk_project(root)
    (root / "input" / "meta.tsv").write_text("genome_id\tfilename\tsource\tgroup_gut\tniche_rumen\n")
    _smk(root, ["preprocess"])
    R = root / "results" / "preprocess"
    assert not [p for p in (R / "batches").iterdir() if p.name.startswith("batch_")]
    assert read_tsv(R / "annotated_metadata.tsv").empty
    sp = read_tsv(R / "species_metadata.tsv")
    assert sp.empty and {"genome_id", "assembly_path", "group_gut", "niche_rumen"} <= set(sp.columns)
    assert (R / "representatives.txt").read_text() == ""
    assert "no data" in (R / "plots" / "preprocess_quality.svg").read_text()


if __name__ == "__main__":
    try:
        suite.main("PREPROCESS_QC TESTS PASSED")
    except SystemExit as exc:
        if exc.code:
            print(f"(temporary directories kept for inspection: {[str(p) for p in _TMP]})")
        raise
    for p in _TMP:
        shutil.rmtree(p, ignore_errors=True)
