"""Synthetic + real-tool test suite for the opt-in intrinsic-disorder stage (leaf 1.1.9).

Run:  bash tests/run_leaf.sh disorder        (prints DISORDER TESTS PASSED on success)

Owned code under test: workflow/rules/disorder.smk, workflow/envs/disorder.yaml,
workflow/scripts/py/{predict_disorder,map_disorder_to_alignment,plot_disorder}.py,
workflow/scripts/R/disorder_stats.R.

Real tools (used when installed; located next to PROM_PYTHON's conda env, or via env vars):
  PROM_DISORDER_PYTHON   python with metapredict==3.1.1   (default <envs>/prom-disorder-metapredict)
  PROM_DISORDER_BLAST    dir holding BLAST+ segmasker      (default <envs>/prom-disorder-blast/bin)
  PROM_RSCRIPT           Rscript with ape + phylolm        (tests/_env.py)
  PROM_TOOLS_BIN         mafft                             (tests/_env.py)
  PROM_SNAKEMAKE         snakemake >= 8                    (default <envs>/prom-smk/bin/snakemake)
IUPred3 is licence-gated and cannot be installed here: its library contract is exercised
through a stub module that mimics the documented iupred(seq, mode[, smoothing]) call.
Threads: every tool is capped at 2 (torch, OMP/MKL, MAFFT).
"""
from __future__ import annotations

import atexit
import json
import math
import os
import py_compile
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True                    # never write __pycache__ outside tests/disorder
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _env import (Suite, py, rscript, run, tmpdir, FIXTURES, have, read_tsv,  # noqa: E402
                  REPO, SCRIPTS_PY, SCRIPTS_R, PYTHON, RSCRIPT, env)

# ---- private scratch (not the sandbox TMPDIR, never the repo root) ----------------------
SCRATCH = Path(__file__).resolve().parent / "_scratch"
shutil.rmtree(SCRATCH, ignore_errors=True)
SCRATCH.mkdir(parents=True)
tempfile.tempdir = str(SCRATCH)
os.environ["TMPDIR"] = str(SCRATCH)
os.environ["MPLCONFIGDIR"] = str(SCRATCH / "mpl")
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
if not os.environ.get("DISORDER_KEEP_SCRATCH"):
    atexit.register(shutil.rmtree, SCRATCH, True)

sys.path.insert(0, str(SCRIPTS_PY))
import predict_disorder as PD  # noqa: E402
import map_disorder_to_alignment as MAP  # noqa: E402
import utils  # noqa: E402

ENVS = Path(PYTHON).resolve().parents[2] if Path(PYTHON).exists() else Path("/nonexistent")
META_PY = os.environ.get("PROM_DISORDER_PYTHON") or str(ENVS / "prom-disorder-metapredict" / "bin" / "python")
BLAST_BIN = os.environ.get("PROM_DISORDER_BLAST") or str(ENVS / "prom-disorder-blast" / "bin")
SNAKEMAKE = os.environ.get("PROM_SNAKEMAKE") or str(ENVS / "prom-smk" / "bin" / "snakemake")
HAVE_META = Path(META_PY).is_file()
HAVE_SEG = (Path(BLAST_BIN) / "segmasker").is_file()
OWNED_PY = ["predict_disorder.py", "map_disorder_to_alignment.py", "plot_disorder.py"]
OWNED_R = ["disorder_stats.R"]
SPEC = {"DUF4011": {"names": ["DUF4011"], "accessions": ["PF13195"], "required": True},
        "DNA2": {"names": ["AAA_11", "AAA_12", "AAA_19", "AAA_30"],
                 "accessions": ["PF13086", "PF13087", "PF13245", "PF13604"], "required": True},
        "REase": {"names": ["MTES_1575"], "accessions": ["PF18741"], "required": True},
        "DUF3320": {"names": ["DUF3320"], "accessions": ["PF11784"], "required": False}}
SPEC_JSON = json.dumps(SPEC)
P53 = ("MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPA"
       "APAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQS"
       "QHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRR"
       "PILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQI"
       "RGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD")
assert len(P53) == 393
REAL_RUNS = []          # real-tool tests that actually ran (printed at the end)

suite = Suite("disorder")


# =========================================================================
# helpers
# =========================================================================
def sh(cmd, extra_env=None, path=None, timeout=900, cwd=None):
    """subprocess with the toolchain env (+ overrides); PATH replaced when `path` given."""
    e = env()
    e.update(extra_env or {})
    if path is not None:
        e["PATH"] = path
    return subprocess.run([str(c) for c in cmd], cwd=str(cwd or REPO), env=e, text=True,
                          capture_output=True, timeout=timeout)


def must_fail(p, needle=None):
    assert p.returncode != 0, (f"expected failure, got exit 0\nstdout={p.stdout[-1500:]}\n"
                               f"stderr={p.stderr[-1500:]}")
    if needle:
        assert re.search(needle, p.stdout + p.stderr, re.I), \
            f"failure message lacks {needle!r}: {(p.stdout + p.stderr)[-1500:]}"


def must_pass(p):
    assert p.returncode == 0, f"exit {p.returncode}\nstdout={p.stdout[-2000:]}\nstderr={p.stderr[-2000:]}"


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def fasta(path, recs):
    return write(path, "".join(f">{k}\n{v}\n" for k, v in recs.items()))


def catalog(path, rows):
    """rows: {locus: (catalytic_status, fragmented, confidence_tier)}"""
    head = "locus_id\tgenome_id\tprotein_id\tlength_aa\tcatalytic_status\tfragmented\tconfidence_tier\n"
    body = "".join(f"{k}\t{k.split('__')[0]}\t{k}\t0\t{c}\t{f}\t{t}\n" for k, (c, f, t) in rows.items())
    return write(path, head + body)


def fake_modules(d):
    """Stub metapredict + torch + iupred3_lib + segmasker mimicking the verified APIs."""
    d = Path(d)
    write(d / "torch" / "__init__.py",
          "import os\n"
          "def set_num_threads(n):\n"
          "    log = os.environ.get('FAKE_TORCH_LOG')\n"
          "    if log:\n"
          "        open(log, 'a').write(str(n) + '\\n')\n")
    write(d / "metapredict" / "__init__.py", r'''
import os
__version__ = "stub"
MODE = os.environ.get("FAKE_META_MODE", "ok")
if MODE == "noimport":
    raise ImportError("simulated: metapredict not installed")
def _s(c):
    return 0.9 if c in "PESKQGDN" else 0.1
def predict_disorder(inputs, version="V3", device=None, **kw):
    assert isinstance(inputs, dict) and version == "V3" and device == "cpu", (version, device)
    if MODE == "crash":
        raise RuntimeError("simulated metapredict crash")
    out = {}
    for k, s in inputs.items():
        sc = [_s(c) for c in s]
        if MODE == "short":
            sc = sc[:-1]
        if MODE == "nan":
            sc[0] = float("nan")
        out[k] = [s, sc]
    return out
def predict_pLDDT(inputs, pLDDT_version="V2", device=None, **kw):
    assert pLDDT_version == "V1", pLDDT_version
    return {k: [s, [30.0 if c in "PESK" else 85.0 for c in s]] for k, s in inputs.items()}
''')
    write(d / "iupred" / "iupred3_lib.py", r'''
import json, os
def iupred(seq, mode="long", smoothing=None):
    log = os.environ.get("FAKE_IUPRED_LOG")
    if log:
        open(log, "a").write(json.dumps({"mode": mode, "smoothing": smoothing}) + "\n")
    m = os.environ.get("FAKE_IUPRED_MODE", "tuple")
    sc = [0.7 if c in "PESKQG" else 0.2 for c in seq]
    if m == "crash":
        raise ValueError("simulated IUPred3 crash")
    if m == "short":
        return (sc[:-1], "")
    if m == "list":
        return sc
    return (sc, "")
''')
    seg = write(d / "bin" / "segmasker", r'''#!/bin/bash
in=""
while [ $# -gt 0 ]; do case "$1" in -in) in="$2"; shift 2;; *) shift;; esac; done
case "${FAKE_SEG_MODE:-ok}" in
  crash) echo "simulated segmasker crash" >&2; exit 3;;
  missing) awk '/^>/{n++; h[n]=$0} END{for(i=1;i<n;i++){print h[i]; print "0 - 2"}}' "$in";;
  beyond) awk '/^>/{print $0; print "0 - 99999"}' "$in";;
  *) awk '/^>/{print $0; print "0 - 2"}' "$in";;
esac
''')
    seg.chmod(0o755)
    return d


def base_path():
    """PATH with the interpreter but without any segmasker."""
    return os.pathsep.join([str(Path(PYTHON).parent), "/usr/bin", "/bin"])


def predict(args, fake_dir=None, extra_env=None, path=None, python=None):
    e = dict(extra_env or {})
    if fake_dir is not None:
        e["PYTHONPATH"] = str(Path(fake_dir))
    return sh([python or PYTHON, SCRIPTS_PY / "predict_disorder.py", *args], extra_env=e,
              path=path if path is not None else base_path())


def bh(p):
    """Benjamini-Hochberg adjusted p-values (independent reference implementation)."""
    n = len(p)
    order = sorted(range(n), key=lambda i: p[i])
    adj, prev = [0.0] * n, 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * n / rank)
        adj[i] = prev
    return adj


def domtbl_line(query, qlen, name, acc, ali_from, ali_to, score=100.0):
    """A 22-field hmmscan --domtblout row (+ description)."""
    f = [name, acc, "100", query, "-", str(qlen), "1e-30", str(score), "0.0", "1", "1",
         "1e-30", "1e-30", str(score), "0.0", "1", "100", str(ali_from), str(ali_to),
         str(ali_from), str(ali_to), "0.95", "desc"]
    return " ".join(f) + "\n"


# =========================================================================
# 1. static checks
# =========================================================================
@suite.test
def test_py_compile_owned():
    for f in OWNED_PY:
        py_compile.compile(str(SCRIPTS_PY / f), cfile=str(SCRATCH / (f + "c")), doraise=True)


@suite.test
def test_r_parse_owned():
    for f in OWNED_R:
        p = sh([RSCRIPT, "-e", f"invisible(parse(file='{SCRIPTS_R / f}')); cat('PARSED')"])
        must_pass(p)
        assert "PARSED" in p.stdout


@suite.test
def test_rule_hygiene():
    p = sh([PYTHON, REPO / "tests" / "check_rules.py", REPO / "workflow" / "rules" / "disorder.smk"])
    must_pass(p)
    assert "RULES HYGIENE OK" in p.stdout
    txt = (REPO / "workflow" / "rules" / "disorder.smk").read_text()
    assert "|| true" not in txt and "nuclease" not in txt
    for rule in ("disorder_predict", "disorder_map", "disorder_predict_bact", "disorder_stats",
                 "plot_disorder"):
        assert re.search(rf"^rule {rule}:", txt, re.M), rule
    assert "rules.family_align.output.aln" in txt, "mapping must use the UNtrimmed alignment"
    assert "family_align.output.trim" not in txt


@suite.test
def test_env_yaml_pins():
    y = (REPO / "workflow" / "envs" / "disorder.yaml").read_text()
    code = "\n".join(l.split("#", 1)[0] for l in y.splitlines())
    assert re.search(r"-\s*metapredict==3\.1\.1\s*$", code, re.M), "metapredict must be pinned to 3.1.1"
    assert "2.6.5" not in code, "metapredict 2.6.5 does not exist on PyPI"
    assert not re.search(r"alphapredict", code, re.I), "an alphaPredict pin conflicts with metapredict"
    assert re.search(r"-\s*blast=2\.17\.0", code) and re.search(r"-\s*pytorch=2\.5\.1=cpu_\*", code)
    assert re.search(r"-\s*c-compiler", code), "metapredict's sdist needs a C compiler on linux-64"
    if HAVE_META:   # the version under test is the pinned one
        p = sh([META_PY, "-c", "import importlib.metadata as m; print(m.version('metapredict'))"])
        must_pass(p)
        assert p.stdout.strip() == "3.1.1", p.stdout


# =========================================================================
# 2. predict_disorder.py -- pure logic with known answers
# =========================================================================
def _rejects(fn, *a):
    try:
        fn(*a)
    except PD.DisorderError:
        return True
    return False


@suite.test
def test_sanitize_sequence():
    assert PD.sanitize_sequence("a", "mkv*") == ("MKV", 0)
    assert PD.sanitize_sequence("a", "MXBZU*") == ("MGNQC", 4)
    for bad in ("MK*V", "MKOV", "MK-V", "MK1", "*"):
        assert _rejects(PD.sanitize_sequence, "a", bad), f"{bad!r} should be rejected"


@suite.test
def test_consensus_rules_and_thresholds():
    meta = [0.6, 0.5, 0.49, 0.9]
    cons, mb, sb = PD.consensus_binary(meta, [0.51, 0.5, 0.9, 0.7], 0.5, 0.5, False, 50)
    assert mb == [1, 1, 0, 1]                 # metapredict: score >= threshold
    assert sb == [1, 0, 1, 1]                 # IUPred: score > threshold
    assert cons == [1, 0, 0, 1]
    cons, mb, sb = PD.consensus_binary(meta, [40, 50, 49.9, 80], 0.5, 0.5, True, 50)
    assert sb == [1, 0, 1, 0] and cons == [1, 0, 0, 0]      # pLDDT: < cutoff
    cons, mb, sb = PD.consensus_binary(meta, None, 0.5, 0.5, False, 50)
    assert cons == mb == [1, 1, 0, 1] and sb == [None] * 4
    for sec, plddt in (([0.5, 0.5, 0.5], False), ([0.5, float("nan"), 0.5, 0.5], False),
                       ([1.5, 0.5, 0.5, 0.5], False), ([120, 1, 1, 1], True)):
        assert _rejects(PD.consensus_binary, meta, sec, 0.5, 0.5, plddt, 50), sec


@suite.test
def test_call_idrs_runs():
    b = [1] * 12 + [0] + [1] * 11 + [0] * 3 + [1] * 15
    assert PD.call_idrs(b, 12) == [(1, 12), (28, 42)]
    assert PD.call_idrs(b, 1) == [(1, 12), (14, 24), (28, 42)]
    assert PD.call_idrs([0, 0], 1) == [] and PD.call_idrs([], 3) == []


@suite.test
def test_summarize_lcr_semantics():
    cons = [1, 1, 0, 0, 1, 1, 1, 1, 0, 0]
    lcr = [1, 0, 0, 0, 0, 0, 0, 0, 1, 1]
    s = PD.summarize(cons, cons, [None] * 10, lcr, 3)
    assert s["frac_disordered"] == 0.6 and s["frac_lcr"] == 0.3
    assert s["frac_disordered_nonlcr"] == round(5 / 7, 4)      # LCR out of numerator AND denominator
    assert s["n_idrs"] == 1 and s["longest_idr"] == 4 and s["frac_in_idrs"] == 0.4
    assert s["frac_disordered_second"] == ""
    s = PD.summarize([1, 1], [1, 1], [1, 0], [1, 1], 1)
    assert s["frac_disordered_nonlcr"] == "" and s["frac_disordered_second"] == 0.5
    assert PD.summarize([], [], [], None, 3) is None


@suite.test
def test_iupred_return_shapes():
    assert PD.first_profile(([0.1, 0.2], "glob text")) == [0.1, 0.2]
    assert PD.first_profile(([0.1, 0.2], None)) == [0.1, 0.2]
    assert PD.first_profile([0.1, 0.2]) == [0.1, 0.2]
    for bad in ("abc", 3.0, {"a": 1}):
        assert _rejects(PD.first_profile, bad), bad


@suite.test
def test_segmasker_interval_parser():
    txt = ">q1 desc\n30 - 45\n136 - 156\n>lcl|q2\n>q3\n0 - 24\n"
    assert PD.parse_segmasker_interval(txt) == {"q1": [(30, 45), (136, 156)], "q2": [],
                                                "q3": [(0, 24)]}
    for bad in ("30 - 45\n", ">q1\n3 to 5\n", ">q1\n9 - 3\n", ">q1\n>q1\n"):
        assert _rejects(PD.parse_segmasker_interval, bad), bad


@suite.test
def test_catalog_and_background_loading():
    d = tmpdir()
    cat = catalog(d / "c.tsv", {"G1__a": ("intact", "False", "A"), "G2__b": ("intact", "", "B"),
                                "G3__c": ("partial", "False", "C"), "G4__d": ("INTACT", "false", "A")})
    c = PD.load_catalog(str(cat))
    assert [k for k in c if PD.is_intact(c[k])] == ["G1__a", "G4__d"]    # unknown fragmented -> out
    write(d / "bad.tsv", "locus_id\tcatalytic_status\nx\tintact\n")
    assert _rejects(PD.load_catalog, str(d / "bad.tsv"))
    assert _rejects(PD.load_catalog, str(d / "none.tsv"))
    write(d / "bg1.tsv", "genome_id\tbackground_frac\nG1\t0.25\n")
    write(d / "bg2.tsv", "genome_id\tbackground_frac_disordered\nG1\t0.1\nG2\t0.2\n")
    assert PD.load_background(str(d / "bg1.tsv")) == {"G1": 0.25}
    assert PD.load_background(str(d / "bg2.tsv")) == {"G1": 0.1, "G2": 0.2}
    assert PD.load_background("") == {}
    for txt in ("genome_id\tfrac\nG1\t0.1\n", "genome_id\tbackground_frac\nG1\t1.5\n",
                "genome_id\tbackground_frac\nG1\t0.1\nG1\t0.2\n", "genome_id\tbackground_frac\nG1\tx\n"):
        write(d / "bgx.tsv", txt)
        assert _rejects(PD.load_background, str(d / "bgx.tsv")), txt


# =========================================================================
# 3. predict_disorder.py -- CLI: empty input, failures, known answers (stub tools)
# =========================================================================
def _cli_inputs(d):
    seqs = {"G1__G1_1_1": "MPESKQAAAAPESKQAAAAPESK*",       # 23 residues once '*' is stripped
            "G2__G2_1_2": "MAAAXAAAAAAAPPPPPPPPPPPPP",        # X -> G
            "G3__G3_1_3": "MPPPPPPPPPPPPPPPPPP"}               # not intact -> skipped
    fasta(d / "p.faa", seqs)
    catalog(d / "c.tsv", {"G1__G1_1_1": ("intact", "False", "A"),
                          "G2__G2_1_2": ("intact", "False", "B"),
                          "G3__G3_1_3": ("degenerate", "False", "C")})
    return seqs


def _outs(d):
    return ["--out-locus", d / "loc.tsv", "--out-residue", d / "res.tsv"]


@suite.test
def test_predict_empty_and_all_filtered():
    d = tmpdir()
    write(d / "empty.faa", "")
    p = predict(["--proteins", d / "empty.faa", "--intact-only", "False", *_outs(d)])
    must_pass(p)
    assert (d / "loc.tsv").read_text() == PD.LOCUS_HEADER + "\n"
    assert (d / "res.tsv").read_text() == PD.RESIDUE_HEADER + "\n"
    fasta(d / "p.faa", {"G3__x": "MPPPP"})
    catalog(d / "c.tsv", {"G3__x": ("partial", "True", "C")})
    p = predict(["--proteins", d / "p.faa", "--catalog", d / "c.tsv", *_outs(d)])
    must_pass(p)       # genuinely empty after the intact filter; no predictor needed
    assert read_tsv(d / "loc.tsv").shape[0] == 0 and "1 skipped" in p.stderr
    must_fail(predict(["--proteins", d / "missing.faa", *_outs(d)]), "not found")


@suite.test
def test_predict_fails_on_missing_or_crashing_tools():
    d = tmpdir()
    fk = fake_modules(d / "fake")
    _cli_inputs(d)
    core = ["--proteins", d / "p.faa", "--catalog", d / "c.tsv"]
    args = [*core, "--mask-lcr", "False", *_outs(d)]
    must_fail(predict(args, fk, {"FAKE_META_MODE": "noimport"}), "metapredict could not be imported")
    must_fail(predict(args, fk, {"FAKE_META_MODE": "crash"}), "predict_disorder failed")
    must_fail(predict(args, fk, {"FAKE_META_MODE": "short"}), "returned 22 scores for a 23-residue")
    must_fail(predict(args, fk, {"FAKE_META_MODE": "nan"}), "outside")
    lcr_args = [*core, "--mask-lcr", "True", *_outs(d)]
    must_fail(predict(lcr_args, fk), "segmasker")          # LCR on, no segmasker -> error
    fpath = os.pathsep.join([str(fk / "bin"), base_path()])
    must_fail(predict(lcr_args, fk, {"FAKE_SEG_MODE": "crash"}, path=fpath), "segmasker failed")
    must_fail(predict(lcr_args, fk, {"FAKE_SEG_MODE": "missing"}, path=fpath), "reported 1 sequences for 2")
    must_fail(predict(lcr_args, fk, {"FAKE_SEG_MODE": "beyond"}, path=fpath), "beyond")
    must_fail(predict(["--proteins", d / "p.faa", *_outs(d)], fk), "needs --catalog")
    catalog(d / "c2.tsv", {"G1__G1_1_1": ("intact", "False", "A")})
    must_fail(predict(["--proteins", d / "p.faa", "--catalog", d / "c2.tsv", *_outs(d)], fk),
              "absent from the catalog")
    write(d / "dup.faa", ">G1__a\nMK\n>G1__a\nMK\n")
    must_fail(predict(["--proteins", d / "dup.faa", "--intact-only", "False", *_outs(d)], fk), "duplicated")
    must_fail(predict([*args, "--second-predictor", "bogus"], fk), "second-predictor")
    must_fail(predict([*args, "--second-predictor", "iupred3"], fk), "needs --iupred3-path")
    must_fail(predict([*args, "--iupred3-path", d / "nope"], fk), "not a directory")
    # fail fast: a bad IUPred3 path / missing segmasker is reported before metapredict runs
    must_fail(predict([*args, "--iupred3-path", d / "nope"], fk, {"FAKE_META_MODE": "crash"}),
              "not a directory")
    must_fail(predict(lcr_args, fk, {"FAKE_META_MODE": "crash"}), "segmasker \\(BLAST\\+\\) is not on PATH")
    must_fail(predict([*core, "--intact-only", "maybe", *_outs(d)], fk), "boolean")
    must_fail(predict([*args, "--meta-threshold", "1.5"], fk), "meta-threshold")


@suite.test
def test_predict_known_answer_alphapredict_fallback():
    """Stub metapredict: P/E/S/K/Q/G/D/N -> 0.9 else 0.1; stub pLDDT: P/E/S/K -> 30 else 85;
    stub segmasker masks residues 1-3 of every sequence."""
    d = tmpdir()
    fk = fake_modules(d / "fake")
    _cli_inputs(d)
    write(d / "bg.tsv", "genome_id\tbackground_frac\nG1\t0.1\nG2\t0.05\n")
    fpath = os.pathsep.join([str(fk / "bin"), base_path()])
    p = predict(["--proteins", d / "p.faa", "--catalog", d / "c.tsv", "--mask-lcr", "True",
                 "--min-idr-len", "3", "--host-background", d / "bg.tsv", "--threads", "2",
                 *_outs(d)], fk, {"FAKE_TORCH_LOG": str(d / "torch.log")}, path=fpath)
    must_pass(p)
    assert (d / "torch.log").read_text().split() == ["2"], "torch threads not capped"
    loc = read_tsv(d / "loc.tsv")
    assert list(loc.columns)[:14] == PD.LOCUS_COLUMNS[:14], "SCHEMA column order changed"
    assert list(loc["locus_id"]) == ["G1__G1_1_1", "G2__G2_1_2"]          # G3 not intact
    r1 = loc.iloc[0]
    # G1 = MPESKQAAAAPESKQAAAAPESK (23): metapredict-disordered = PESKQ x2 + PESK = 14;
    # pLDDT < 50 = PESK x3 = 12; consensus = 12; LCR (stub) = residues 1-3 (M,P,E)
    assert r1["length_aa"] == "23" and r1["predictors"] == "metapredict+alphapredict"
    assert float(r1["frac_disordered"]) == round(12 / 23, 4)
    assert float(r1["frac_disordered_metapredict"]) == round(14 / 23, 4)
    assert float(r1["frac_disordered_second"]) == round(12 / 23, 4)
    assert float(r1["frac_lcr"]) == round(3 / 23, 4)
    assert float(r1["frac_disordered_nonlcr"]) == round(10 / 20, 4)     # P,E of 1-3 masked
    assert r1["n_idrs"] == "3" and r1["longest_idr"] == "4"
    assert float(r1["disorder_minus_background"]) == round(12 / 23 - 0.1, 4)
    assert r1["host_background"] == "0.1" and r1["confidence_tier"] == "A"
    assert r1["intact"] == "True" and r1["n_nonstandard"] == "0"
    assert "metapredict stub V3 (score>=0.5)" in r1["predictor_detail"]
    assert "pLDDT V1 = alphaPredict network (pLDDT<50)" in r1["predictor_detail"]
    assert "min_idr_len=3; LCR=segmasker" in r1["predictor_detail"]
    r2 = loc.iloc[1]
    assert r2["n_nonstandard"] == "1" and r2["confidence_tier"] == "B" and r2["length_aa"] == "25"
    res = read_tsv(d / "res.tsv")
    assert list(res.columns) == PD.RESIDUE_COLUMNS
    g2 = res[res["locus_id"] == "G2__G2_1_2"]
    assert "".join(g2["aa"]) == "MAAAGAAAAAAAPPPPPPPPPPPPP"            # X converted, positions kept
    assert list(g2["position"].astype(int)) == list(range(1, 26))
    g1 = res[res["locus_id"] == "G1__G1_1_1"]
    assert list(g1["lcr"].astype(int))[:4] == [1, 1, 1, 0]
    assert set(g1.loc[g1["aa"].isin(list("PESK")), "second_score"]) == {"30"}
    write(d / "bg_partial.tsv", "genome_id\tbackground_frac\nG1\t0.1\n")
    must_fail(predict(["--proteins", d / "p.faa", "--catalog", d / "c.tsv", "--mask-lcr", "False",
                       "--host-background", d / "bg_partial.tsv", *_outs(d)], fk), "lacks 1 scored")


@suite.test
def test_predict_iupred3_library_contract():
    d = tmpdir()
    fk = fake_modules(d / "fake")
    _cli_inputs(d)
    log = d / "iup.log"
    base = ["--proteins", d / "p.faa", "--catalog", d / "c.tsv", "--mask-lcr", "False",
            "--iupred3-path", fk / "iupred", "--min-idr-len", "3", *_outs(d)]
    must_pass(predict([*base, "--iupred3-smoothing", "medium"], fk, {"FAKE_IUPRED_LOG": str(log)}))
    calls = [json.loads(l) for l in log.read_text().splitlines()]
    assert calls == [{"mode": "long", "smoothing": "medium"}] * 2, calls
    loc = read_tsv(d / "loc.tsv")
    assert set(loc["predictors"]) == {"metapredict+iupred3"}
    assert "IUPred3 long smoothing=medium (score>0.5)" in loc.iloc[0]["predictor_detail"]
    # G1: metapredict-disordered PESKQ x2 + PESK (14); IUPred > 0.5 for PESKQG -> same 14
    assert float(loc.iloc[0]["frac_disordered"]) == round(14 / 23, 4)
    res = read_tsv(d / "res.tsv")
    assert set(res.loc[res["aa"] == "P", "second_score"]) == {"0.7"}
    log.unlink()
    must_pass(predict(base, fk, {"FAKE_IUPRED_LOG": str(log), "FAKE_IUPRED_MODE": "list"}))
    assert json.loads(log.read_text().splitlines()[0])["smoothing"] is None   # library default
    must_fail(predict(base, fk, {"FAKE_IUPRED_MODE": "short"}), "IUPred3 returned 22 scores")
    must_fail(predict(base, fk, {"FAKE_IUPRED_MODE": "crash"}), "iupred3_lib.iupred")
    must_pass(predict([*base, "--second-predictor", "none"], fk))
    assert set(read_tsv(d / "loc.tsv")["predictors"]) == {"metapredict_only"}


# =========================================================================
# 4. real tools: segmasker and metapredict 3.1.1 (+ pLDDT V1)
# =========================================================================
@suite.test
def test_real_segmasker_positions():
    if not HAVE_SEG:
        print("    (segmasker not installed: real segmasker run skipped)")
        return
    seqs = {"G1__polyQ": "MKVLITGAGS" + "Q" * 25 + "GIGLELAKRLLEEGHEVIVTDLNPEK",
            "G2__p53": P53, "G3__none": "MKVLITGAGSGIGLELAKRLLE"}
    d = tmpdir()
    old = os.environ["PATH"]
    os.environ["PATH"] = BLAST_BIN + os.pathsep + old
    try:
        got = PD.low_complexity_positions(seqs, str(d))
    finally:
        os.environ["PATH"] = old
    fasta(d / "in.faa", seqs)
    p = sh([Path(BLAST_BIN) / "segmasker", "-in", d / "in.faa", "-outfmt", "fasta"])
    must_pass(p)
    masked, cur = {}, None
    for line in p.stdout.splitlines():
        if line.startswith(">"):
            cur = line[1:].split()[0]
            masked[cur] = ""
        else:
            masked[cur] += line.strip()
    assert set(masked) == set(seqs)
    for k, s in masked.items():   # lower case = masked; independent of the interval parser
        exp = {i + 1 for i, c in enumerate(s) if c.islower()}
        assert got[k] == exp, (k, sorted(got[k] ^ exp)[:10])
    assert set(range(11, 36)) <= got["G1__polyQ"], "poly-Q stretch not masked"
    assert got["G2__p53"], "p53 has SEG low-complexity segments"
    REAL_RUNS.append("segmasker")


@suite.test
def test_real_metapredict_end_to_end():
    if not HAVE_META:
        print("    (metapredict env not installed: real metapredict run skipped)")
        return
    d = tmpdir()
    fasta(d / "p.faa", {"G1__p53": P53 + "*",
                        "G2__polyQ": "MKVLITGAGS" + "Q" * 25 + "GIGLELAKRLLEEGHEVIVTDLNPEK"})
    catalog(d / "c.tsv", {"G1__p53": ("intact", "False", "A"), "G2__polyQ": ("intact", "False", "A")})
    dirs = ([BLAST_BIN] if HAVE_SEG else []) + [str(Path(META_PY).parent), "/usr/bin", "/bin"]
    path = os.pathsep.join(dirs)
    p = sh([META_PY, SCRIPTS_PY / "predict_disorder.py", "--proteins", d / "p.faa",
            "--catalog", d / "c.tsv", "--mask-lcr", str(HAVE_SEG), "--threads", "2", *_outs(d)],
           path=path, timeout=600)
    must_pass(p)
    assert "metapredict 3.1" in p.stderr and "V3 (score>=0.5)" in p.stderr and "device=cpu" in p.stderr, \
        p.stderr[-800:]
    res = read_tsv(d / "res.tsv")
    g = res[res["locus_id"] == "G1__p53"]
    assert len(g) == 393 and "".join(g["aa"]) == P53
    meta = [float(v) for v in g["meta_score"]]
    sec = [float(v) for v in g["second_score"]]
    cons = [int(v) for v in g["consensus"]]

    def frac(a, b):
        return sum(1 for v in meta[a:b] if v >= 0.5) / (b - a)
    # p53 (P04637): transactivation domain and C-terminal tail disordered, DNA-binding
    # domain folded (measured with metapredict V3: 1.00 / 0.053 / 1.00)
    assert frac(0, 60) > 0.9 and frac(100, 290) < 0.15 and frac(360, 393) > 0.9, \
        (frac(0, 60), frac(100, 290), frac(360, 393))
    assert all(c == int(m >= 0.5 and s < 50) for c, m, s in zip(cons, meta, sec))
    q = sh([META_PY, "-c", "import json,sys,metapredict as m;s=sys.argv[1];"
            "r=m.predict_disorder({'x':s},version='V3',device='cpu')['x'][1];"
            "p=m.predict_pLDDT({'x':s},pLDDT_version='V1',device='cpu')['x'][1];"
            "print(json.dumps([[float(v) for v in r],[float(v) for v in p]]))", P53],
           path=path, timeout=300)
    must_pass(q)
    direct, plddt = json.loads(q.stdout)
    assert max(abs(a - b) for a, b in zip(direct, meta)) < 1e-4
    assert max(abs(a - b) for a, b in zip(plddt, sec)) < 1e-3
    loc = read_tsv(d / "loc.tsv")
    assert set(loc["predictors"]) == {"metapredict+alphapredict"}
    assert int(loc.iloc[0]["length_aa"]) == 393                       # trailing '*' stripped
    assert loc.iloc[0]["predictor_detail"].startswith("metapredict 3.1.01 V3"), loc.iloc[0]["predictor_detail"]
    assert float(loc.iloc[0]["frac_disordered"]) == round(sum(cons) / 393, 4)
    if HAVE_SEG:
        gq = res[res["locus_id"] == "G2__polyQ"]
        assert all(int(v) for v in list(gq["lcr"])[10:35]), "poly-Q not LCR-masked"
        REAL_RUNS.append("segmasker (inside predict_disorder)")
    REAL_RUNS.append("metapredict 3.1.1 V3 + pLDDT V1")


# =========================================================================
# 5. map_disorder_to_alignment.py
# =========================================================================
def _map_case(d, residue_rows=None, aln=None, dom=None):
    """3 aligned loci, 2 scored (G3 aligned but not scored). Known answers in the tests."""
    aln = aln or {"G1__a": "MK-VLDE--", "G2__b": "MKAVL-EGG", "G3__c": "-KAVLDEGG"}
    fasta(d / "aln.faa", aln)
    if residue_rows is None:
        residue_rows = []
        for locus, seq, cons in (("G1__a", "MKVLDE", [1, 1, 0, 0, 1, 1]),
                                 ("G2__b", "MKAVLEGG", [0, 1, 0, 0, 0, 1, 1, 1])):
            for i, (aa, c) in enumerate(zip(seq, cons), 1):
                residue_rows.append(f"{locus}\t{i}\t{0.9 if c else 0.1}\t\t{c}\t0\t{aa}\n")
    write(d / "res.tsv", PD.RESIDUE_HEADER + "\n" + "".join(residue_rows))
    if dom is None:
        dom = (domtbl_line("G1__a", 6, "DUF4011", "PF13195.12", 1, 3)
               + domtbl_line("G1__a", 6, "MTES_1575", "PF18741.8", 4, 6)
               + domtbl_line("G2__b", 8, "DUF4011", "PF13195.12", 1, 4)
               + domtbl_line("G3__c", 7, "AAA_11", "PF13086.13", 1, 7)
               + domtbl_line("G2__b", 8, "Pkinase", "PF00069.32", 5, 8))    # not a spec category
    write(d / "d.domtbl", "# hmmscan\n" + dom)
    return ["--residue", d / "res.tsv", "--alignment", d / "aln.faa", "--domtbl", d / "d.domtbl",
            "--domains", SPEC_JSON, "--out", d / "bycol.tsv"]


def _residue_rows(d):
    return (d / "res.tsv").read_text().splitlines(True)[1:]


@suite.test
def test_map_known_answer():
    d = tmpdir()
    must_pass(py("map_disorder_to_alignment.py", *_map_case(d), check=False))
    t = read_tsv(d / "bycol.tsv")
    assert list(t.columns) == MAP.HEADER_COLUMNS and list(t.columns)[:5] == [
        "alignment_col", "n_loci", "mean_disorder", "frac_loci_disordered", "domain"]
    rows = {int(r["alignment_col"]): r for _, r in t.iterrows()}
    # col: 1 M/M/-  2 K/K/K  3 -/A/A  4 V/V/V  5 L/L/L  6 D/-/D  7 E/E/E  8 -/G/G  9 -/G/G
    assert {c: int(r["n_loci"]) for c, r in rows.items()} == {1: 2, 2: 2, 3: 1, 4: 2, 5: 2, 6: 1,
                                                              7: 2, 8: 1, 9: 1}
    # G1 residues M1 K2 V3 L4 D5 E6 sit in cols 1,2,4,5,6,7 with consensus 1,1,0,0,1,1;
    # G2 residues M1 K2 A3 V4 L5 E6 G7 G8 sit in cols 1,2,3,4,5,7,8,9 with 0,1,0,0,0,1,1,1
    assert {c: float(r["mean_disorder"]) for c, r in rows.items()} == {
        1: 0.5, 2: 1.0, 3: 0.0, 4: 0.0, 5: 0.0, 6: 1.0, 7: 1.0, 8: 1.0, 9: 1.0}
    assert all(r["mean_disorder"] == r["frac_loci_disordered"] for r in rows.values())
    # domain votes over the 3 annotated loci. G1: DUF4011 res 1-3 (cols 1,2,4), REase res
    # 4-6 (cols 5,6,7); G2: DUF4011 res 1-4 (cols 1-4), Pkinase (no category) res 5-8;
    # G3: AAA_11 = DNA2 res 1-7 (cols 2-8), res 8 (col 9) unannotated.
    # col1 DUF4011 2/2; col2 DUF4011 2/3; col3 DUF4011 vs DNA2 1:1 TIE -> NA (support 0.5);
    # col4 DUF4011 2/3; col5 REase 1, DNA2 1 of 3 -> NA; col6 REase vs DNA2 1:1 TIE -> NA;
    # col7 1/3 -> NA; col8 DNA2 1 of 2 (G2 uncategorised) = 0.5 -> DNA2; col9 no vote -> NA
    assert {c: r["domain"] for c, r in rows.items()} == {
        1: "DUF4011", 2: "DUF4011", 3: "NA", 4: "DUF4011", 5: "NA", 6: "NA", 7: "NA",
        8: "DNA2", 9: "NA"}, {c: r["domain"] for c, r in rows.items()}
    assert {c: float(r["domain_support"]) for c, r in rows.items()} == {
        1: 1.0, 2: round(2 / 3, 4), 3: 0.5, 4: round(2 / 3, 4), 5: round(1 / 3, 4), 6: 0.5,
        7: round(1 / 3, 4), 8: 0.5, 9: 0.0}
    assert float(rows[3]["occupancy"]) == 0.5 and float(rows[1]["occupancy"]) == 1.0
    assert float(rows[2]["mean_meta_score"]) == 0.9 and float(rows[3]["mean_meta_score"]) == 0.1


@suite.test
def test_map_validation_failures_and_empty():
    d = tmpdir()
    args = _map_case(d)
    bad = [r.replace("\t3\t0.1\t\t0\t0\tV", "\t3\t0.1\t\t0\t0\tW") for r in _residue_rows(d)]
    must_fail(py("map_disorder_to_alignment.py", *_map_case(d, residue_rows=bad), check=False),
              "alignment residue 3")
    must_fail(py("map_disorder_to_alignment.py",
                 *_map_case(d, aln={"G1__a": "MK-VLD---", "G2__b": "MKAVL-EGG", "G3__c": "-KAVLDEGG"}),
                 check=False), "5 residues in the alignment but 6")
    must_fail(py("map_disorder_to_alignment.py",
                 *_map_case(d, aln={"G2__b": "MKAVL-EGG", "G3__c": "-KAVLDEGG", "G4__d": "MKAVLDEGG"}),
                 check=False), "absent from the alignment")
    must_fail(py("map_disorder_to_alignment.py",
                 *_map_case(d, aln={"G1__a": "MKVLDE", "G2__b": "MKAVLEGG", "G3__c": "KAVLDEGG"}),
                 check=False), "not an alignment")
    # ragged with 2 sequences = family_align's < 3 pass-through -> header-only, exit 0
    must_pass(py("map_disorder_to_alignment.py",
                 *_map_case(d, aln={"G1__a": "MKVLDE", "G2__b": "MKAVLEGG"}), check=False))
    assert (d / "bycol.tsv").read_text() == MAP.HEADER + "\n"
    must_fail(py("map_disorder_to_alignment.py",
                 *_map_case(d, dom=domtbl_line("G2__b", 8, "Pkinase", "PF00069.32", 1, 8)),
                 check=False), "domain track cannot be built")
    write(d / "res0.tsv", PD.RESIDUE_HEADER + "\n")
    args = _map_case(d)
    args[1] = d / "res0.tsv"
    must_pass(py("map_disorder_to_alignment.py", *args, check=False))
    assert (d / "bycol.tsv").read_text() == MAP.HEADER + "\n"
    # an aligner that keeps the trailing '*' is tolerated
    must_pass(py("map_disorder_to_alignment.py",
                 *_map_case(d, aln={"G1__a": "MK-VLDE*-", "G2__b": "MKAVL-EGG", "G3__c": "-KAVLDEGG"}),
                 check=False))
    assert int(read_tsv(d / "bycol.tsv").iloc[-1]["alignment_col"]) == 9


@suite.test
def test_map_real_fixture_with_mafft():
    """Real ProA homologues + their real hmmscan domtbl, aligned with MAFFT (2 threads): the
    consensus domain track must follow the canonical N->C order, and synthetic disorder
    (1 outside each locus' own domains) must be lower in domain columns."""
    if not have("mafft"):
        print("    (mafft not installed: fixture alignment test skipped)")
        return
    d = tmpdir()
    recs = dict(utils.read_fasta(str(FIXTURES / "proA_homologs.faa")))
    run(f"mafft --retree 2 --thread 2 --quiet '{FIXTURES / 'proA_homologs.faa'}' > '{d / 'aln.faa'}'")
    hits = utils.parse_domtbl(str(FIXTURES / "proA_homologs.pfam_subset.domtbl"))
    spec = utils.load_domain_spec(SPEC_JSON)
    lines = []
    for locus, seq in recs.items():
        clean, _ = PD.sanitize_sequence(locus, seq)
        cats = MAP.locus_residue_categories(hits.get(locus, []), spec)
        for i, aa in enumerate(clean, 1):
            c = 0 if i in cats else 1
            lines.append(f"{locus}\t{i}\t{0.8 if c else 0.2}\t\t{c}\t0\t{aa}\n")
    write(d / "res.tsv", PD.RESIDUE_HEADER + "\n" + "".join(lines))
    must_pass(py("map_disorder_to_alignment.py", "--residue", d / "res.tsv", "--alignment",
                 d / "aln.faa", "--domtbl", FIXTURES / "proA_homologs.pfam_subset.domtbl",
                 "--domains", SPEC_JSON, "--out", d / "bycol.tsv", check=False))
    t = read_tsv(d / "bycol.tsv")
    cols = {c: sorted(int(x) for x in t.loc[t["domain"] == c, "alignment_col"]) for c in SPEC}
    assert all(cols[c] for c in SPEC), {c: len(v) for c, v in cols.items()}
    med = [cols[c][len(cols[c]) // 2] for c in SPEC]
    assert med == sorted(med), f"domain track not in canonical N->C order: {dict(zip(SPEC, med))}"
    y = t["mean_disorder"].astype(float)
    dom_mask = t["domain"] != "NA"
    assert y[dom_mask].mean() < y[~dom_mask].mean()
    assert int(t["n_loci"].astype(int).max()) <= len(recs)
    REAL_RUNS.append("mafft + real hmmscan domtbl")


# =========================================================================
# 6. plot_disorder.py
# =========================================================================
def _plot_inputs(d):
    rows = ["\t".join(MAP.HEADER_COLUMNS) + "\n"]
    doms = ["DUF4011"] * 20 + ["NA"] * 10 + ["DNA2"] * 20 + ["REase"] * 10 + ["DUF3320"] * 10 + ["NA"] * 5
    for i, dom in enumerate(doms, 1):
        v = 0.8 if dom == "NA" else 0.1
        occ = 0.2 if i in (21, 22) else 1.0
        rows.append(f"{i}\t5\t{v}\t{v}\t{dom}\t{occ}\t1.0\t{v}\n")
    write(d / "bycol.tsv", "".join(rows))
    write(d / "loc.tsv", "locus_id\tgenome_id\tfrac_disordered\tfrac_disordered_nonlcr\tconfidence_tier\n"
          "GA__1\tGA\t0.2\t0.2\tA\nGB__1\tGB\t0.3\t0.3\tB\nGC__1\tGC\t0.25\t0.25\tA\n"
          "GD__1\tGD\t0.22\t0.22\tA\nGE__1\tGE\t0.21\t0.21\tA\nGA__2\tGA\t0.9\t0.9\tC\n")
    write(d / "master.tsv", "genome_id\tpass_primary\tmethanogen\tis_anme\tprometheus_tier\n"
          "GA\tTrue\tTrue\tFalse\tA\nGB\tTrue\tTrue\tFalse\tB\nGC\tTrue\tFalse\tFalse\tA\n"
          "GD\tTrue\tFalse\tTrue\tA\nGE\tTrue\tFalse\tFalse\tA\nGF\tTrue\tFalse\tFalse\t\n")
    write(d / "bact.tsv", "locus_id\tgenome_id\tfrac_disordered_nonlcr\tconfidence_tier\n"
          "B1__1\tB1\t0.3\tA\nB2__1\tB2\t0.35\tA\n")


def _plot(d, out, *extra):
    return py("plot_disorder.py", "--by-column", d / "bycol.tsv", "--per-locus", d / "loc.tsv",
              "--master", d / "master.tsv", "--outdir", out, *extra, check=False)


@suite.test
def test_plot_outputs_categories_and_tier_rule():
    d = tmpdir()
    _plot_inputs(d)
    out = d / "plots"
    must_pass(_plot(d, out, "--bact-per-locus", d / "bact.tsv", "--domains", SPEC_JSON,
                    "--formats", "png,svg", "--tiers", "A"))
    for fig in ("disorder_along_protein", "disorder_by_group", "disorder_linker_vs_domain"):
        for ext in ("png", "svg"):
            f = out / f"{fig}.{ext}"
            assert f.is_file() and f.stat().st_size > 1000, f
    svg = (out / "disorder_along_protein.svg").read_text()
    for c in SPEC:
        assert c in svg, f"legend lacks spec category {c}"
    assert "nuclease" not in svg
    grp = (out / "disorder_by_group.svg").read_text()
    # tier A: GB (best tier B) EXCLUDED + counted; GA keeps only its tier-A locus;
    # groups: non-methanogen GC,GE (n=2), ANME GD (n=1), methanogen GA (n=1)
    assert "1 genome(s) outside the tier set excluded" in grp
    assert "(n=2)" in grp and "(n=1)" in grp
    must_pass(_plot(d, out, "--domains", SPEC_JSON, "--formats", "svg", "--tiers", "A,B"))
    assert "0 genome(s) outside the tier set excluded" in (out / "disorder_by_group.svg").read_text()
    # categories come from the spec (no hardcoded list): a custom spec adds its own category
    spec2 = dict(SPEC, XYZ={"names": ["XYZ"], "accessions": []})
    write(d / "bycol.tsv", (d / "bycol.tsv").read_text().replace("DUF3320", "XYZ"))
    must_pass(_plot(d, out, "--domains", json.dumps(spec2), "--formats", "svg"))
    assert "XYZ" in (out / "disorder_along_protein.svg").read_text()
    # a domain value that is not a spec category is an error, not a silent grey bar
    must_fail(_plot(d, out, "--domains", SPEC_JSON, "--formats", "svg"), "not categories")


@suite.test
def test_plot_tier_filter_regression():
    """Contract v2 item 6: a B-only genome is EXCLUDED under tier A and PRESENT under A+B."""
    import pandas as pd
    import plot_disorder as PL
    pl = pd.DataFrame({"genome_id": ["GA", "GB", "GA"], "confidence_tier": ["A", "B", "C"],
                       "frac_disordered_nonlcr": ["0.2", "0.3", "0.9"]})
    m = pd.DataFrame({"genome_id": ["GA", "GB", "GZ"], "prometheus_tier": ["A", "B", ""]})
    a, n_a = PL.tier_filter_loci(pl, m, ["A"])
    ab, n_ab = PL.tier_filter_loci(pl, m, ["A", "B"])
    assert list(a["genome_id"]) == ["GA"] and n_a == 1          # GB excluded + counted; GZ absent
    assert sorted(ab["genome_id"]) == ["GA", "GB"] and n_ab == 0
    assert list(a["frac_disordered_nonlcr"]) == ["0.2"]          # GA's tier-C locus not averaged


@suite.test
def test_plot_empty_inputs_placeholders():
    d = tmpdir()
    write(d / "bycol.tsv", MAP.HEADER + "\n")
    write(d / "loc.tsv", PD.LOCUS_HEADER + "\n")
    write(d / "master.tsv", "genome_id\tmethanogen\n")
    must_pass(_plot(d, d / "plots", "--formats", "png"))
    assert sorted(f.name for f in (d / "plots").iterdir()) == [
        "disorder_along_protein.png", "disorder_by_group.png", "disorder_linker_vs_domain.png"]
    must_fail(py("plot_disorder.py", "--by-column", d / "nope.tsv", "--per-locus", d / "loc.tsv",
                 "--master", d / "master.tsv", "--outdir", d / "plots", check=False), "not found")


# =========================================================================
# 7. disorder_stats.R (real ape + phylolm)
# =========================================================================
def _balanced_newick(tips):
    """Deterministic rooted, non-ultrametric bifurcating tree."""
    nodes = [f"{t}:{0.05 + 0.01 * (i % 5):.3f}" for i, t in enumerate(tips)]
    lvl = 0
    while len(nodes) > 1:
        lvl += 1
        nxt = [f"({nodes[i]},{nodes[i + 1]}):{0.04 * lvl:.3f}" for i in range(0, len(nodes) - 1, 2)]
        if len(nodes) % 2:
            nxt.append(nodes[-1])
        nodes = nxt
    return nodes[0].rsplit(":", 1)[0] + ";"


def _stats_case(d, n=24):
    """24 genomes: G11 pass_primary False, G12 no tree tip, G10 tier '' (no Prometheus),
    G07/G08 best tier B, G09 tier C, G03/G05 ANME; methanogen = even index;
    niche_hot = second half; disorder = 0.10 + 0.05*methanogen + 0.03*hot + small wiggle."""
    tips = [f"G{i:02d}" for i in range(1, n + 1)]
    write(d / "tree.nwk", _balanced_newick(tips) + "\n")
    mrows, lrows = [], []
    for i, g in enumerate(tips, 1):
        meth = i % 2 == 0
        anme = i in (3, 5)
        tier = "B" if i in (7, 8) else ("C" if i == 9 else ("" if i == 10 else "A"))
        hot = 1 if i > n // 2 else 0
        marine = 1 if i % 3 == 0 else 0
        mrows.append(f"{g}\t{'' if i == 12 else g}\t{'False' if i == 11 else 'True'}\t{tier}\t"
                     f"{meth}\t{anme}\t{marine}\t{hot}\n")
        if tier:
            y = 0.10 + 0.05 * meth + 0.03 * hot + 0.004 * math.sin(i * 1.7)
            lrows.append(f"{g}__{g}_1_1\t{g}\t{y:.4f}\t{y:.4f}\t{tier}\n")
            if i == 1:   # extra tier-C locus of a tier-A genome: not averaged under tier A
                lrows.append(f"{g}__{g}_1_2\t{g}\t0.9000\t0.9000\tC\n")
    write(d / "master.tsv", "genome_id\ttree_tip\tpass_primary\tprometheus_tier\tmethanogen\t"
          "is_anme\tgroup_marine\tniche_hot\n" + "".join(mrows))
    write(d / "dis.tsv", "locus_id\tgenome_id\tfrac_disordered\tfrac_disordered_nonlcr\tconfidence_tier\n"
          + "".join(lrows))
    return tips


@suite.test
def test_stats_known_answer_and_tier_rule():
    d = tmpdir()
    _stats_case(d)
    p = rscript("disorder_stats.R", "--master", d / "master.tsv", "--disorder", d / "dis.tsv",
                "--tree", d / "tree.nwk", "--temp-cols", "niche_hot", "--tier-sets", "A;A,B",
                "--out", d / "out.tsv", check=False)
    must_pass(p)
    t = read_tsv(d / "out.tsv")
    assert list(t.columns) == ["predictor", "level", "model", "estimate", "std_error", "p_value",
                               "p_adj", "n", "note", "response", "tier_set", "n_level1",
                               "n_excluded_tier", "lambda", "ci_low", "ci_high"]
    assert list(t["predictor"]) == ["phylogenetic_signal", "methanogen", "group_marine"] * 2
    a = t[t["tier_set"] == "A"].set_index("predictor")
    ab = t[t["tier_set"] == "A,B"].set_index("predictor")
    # population 24 - G11 (pass_primary) - G12 (no tip) = 22; G10 has no Prometheus;
    # tier A: G07, G08 (B) and G09 (C) EXCLUDED -> 3, analysed 18
    # tier A,B: only G09 excluded -> 1, analysed 20 (the B-only genomes are PRESENT)
    assert a.loc["phylogenetic_signal", "n"] == "18" and a.loc["methanogen", "n_excluded_tier"] == "3"
    assert ab.loc["phylogenetic_signal", "n"] == "20" and ab.loc["methanogen", "n_excluded_tier"] == "1"
    assert a.loc["methanogen", "n"] == "16" and ab.loc["methanogen", "n"] == "18"   # minus 2 ANME
    assert "2 ANME genome(s) excluded" in a.loc["methanogen", "note"]
    assert "niche_hot" in a.loc["methanogen", "model"] and "niche_hot" not in list(t["predictor"])
    assert a.loc["methanogen", "response"] == "frac_disordered_nonlcr"
    est = float(a.loc["methanogen", "estimate"])
    assert abs(est - 0.05) < 0.01, est                      # simulated effect 0.05
    assert float(a.loc["methanogen", "ci_low"]) < est < float(a.loc["methanogen", "ci_high"])
    # identical to an independent phylolm fit on the expected tier-A data set
    q = sh([RSCRIPT, "-e", f"""
suppressPackageStartupMessages({{library(ape); library(phylolm)}})
tr <- read.tree('{d / "tree.nwk"}'); m <- read.delim('{d / "master.tsv"}', colClasses='character')
x <- read.delim('{d / "dis.tsv"}', colClasses='character'); x <- x[x$confidence_tier == 'A', ]
m <- m[m$pass_primary == 'True' & m$tree_tip != '' & m$prometheus_tier == 'A' & m$is_anme != 'True', ]
m$y <- as.numeric(x$frac_disordered_nonlcr[match(m$genome_id, x$genome_id)])
m$pred <- as.integer(m$methanogen == 'True'); m$niche_hot <- as.integer(m$niche_hot)
t2 <- keep.tip(tr, m$tree_tip); rownames(m) <- m$tree_tip; m <- m[t2$tip.label, ]
f <- phylolm(y ~ pred + niche_hot, data=m, phy=t2, model='lambda'); s <- summary(f)$coefficients
q <- qt(0.975, f$n - f$d)
pchk <- 2 * pt(-abs(s['pred','t.value']), f$n - f$d)
cat(sprintf('%.12g %.12g %.12g %.12g %.12g %.12g', s['pred','Estimate'], s['pred','StdErr'],
            s['pred','p.value'], f$optpar, s['pred','Estimate'] - q * s['pred','StdErr'], pchk))
"""])
    must_pass(q)
    e2, se2, p2, lam2, lo2, pchk = (float(v) for v in q.stdout.split())
    assert math.isclose(p2, pchk, rel_tol=1e-8), "phylolm's t-test df is not n - d"
    assert math.isclose(est, e2, rel_tol=1e-8)
    assert math.isclose(float(a.loc["methanogen", "std_error"]), se2, rel_tol=1e-8)
    assert math.isclose(float(a.loc["methanogen", "p_value"]), p2, rel_tol=1e-6)
    assert math.isclose(float(a.loc["methanogen", "lambda"]), lam2, rel_tol=1e-6)
    assert math.isclose(float(a.loc["methanogen", "ci_low"]), lo2, rel_tol=1e-8)
    # lambda LRT: p = 0.5 * P(chi2_1 > stat); stat recorded in the note
    from scipy.stats import chi2
    stat = float(re.search(r"LRT stat ([0-9.eE+-]+)", a.loc["phylogenetic_signal", "note"]).group(1))
    pexp = 0.5 * chi2.sf(stat, 1) if stat > 0 else 1.0
    assert math.isclose(float(a.loc["phylogenetic_signal", "p_value"]), pexp, rel_tol=1e-3)
    pv = [float(v) for v in t["p_value"]]
    for got, exp in zip(t["p_adj"].astype(float), bh(pv)):
        assert math.isclose(got, exp, rel_tol=1e-9), (got, exp)
    # phylolm warnings are captured into the note (none leak to stderr); a lambda estimate
    # at its bound is flagged
    assert "Warn" not in p.stderr and "warning message" not in p.stderr.lower(), p.stderr[-600:]
    for _, r in t.iterrows():
        if r["lambda"] not in ("", "NA") and (float(r["lambda"]) > 0.999999 or float(r["lambda"]) < 1e-6) \
                and r["predictor"] != "phylogenetic_signal":
            assert "bound" in r["note"], r["note"]
    REAL_RUNS.append("phylolm (R)")


@suite.test
def test_stats_lambda0_loglik_matches_phylolm():
    """The closed-form lambda = 0 log-likelihood used by the LRT equals phylolm's own fit
    with lambda pinned near 0 (same tip-variance-preserving transform)."""
    d = tmpdir()
    _stats_case(d)
    src = (SCRIPTS_R / "disorder_stats.R").read_text()
    start = src.index("loglik_lambda0 <- function")
    fn = src[start:src.index("\n}\n", start) + 3]          # just that function
    q = sh([RSCRIPT, "-e", f"""
suppressPackageStartupMessages({{library(ape); library(phylolm)}})
{fn}
tr <- read.tree('{d / "tree.nwk"}'); set.seed(3); y <- rnorm(Ntip(tr), 0.2, 0.05)
d <- data.frame(y = y, row.names = tr$tip.label)
f <- phylolm(y ~ 1, data = d, phy = tr, model = 'lambda', lower.bound = 1e-7, upper.bound = 1e-6,
             starting.value = 5e-7)
cat(sprintf('%.10g %.10g', loglik_lambda0(d$y, tr), f$logLik))
"""])
    must_pass(q)
    a, b = (float(v) for v in q.stdout.split())
    assert math.isclose(a, b, rel_tol=1e-5, abs_tol=1e-4), (a, b)


@suite.test
def test_stats_skips_low_prevalence_and_handles_empty():
    d = tmpdir()
    _stats_case(d)
    must_pass(rscript("disorder_stats.R", "--master", d / "master.tsv", "--disorder", d / "dis.tsv",
                      "--tree", d / "tree.nwk", "--tier-sets", "A", "--min-group", "9",
                      "--out", d / "out.tsv", check=False))
    t = read_tsv(d / "out.tsv").set_index("predictor")
    assert t.loc["group_marine", "p_value"] == "NA" and "skipped" in t.loc["group_marine", "note"]
    assert "niche_hot" in list(t.index)         # no --temp-cols -> niche_hot is a predictor
    write(d / "dis0.tsv", "locus_id\tgenome_id\tfrac_disordered_nonlcr\tconfidence_tier\n")
    for dis, tree in ((d / "dis0.tsv", d / "tree.nwk"), (d / "dis.tsv", write(d / "t0.nwk", "();\n"))):
        must_pass(rscript("disorder_stats.R", "--master", d / "master.tsv", "--disorder", dis,
                          "--tree", tree, "--out", d / "o0.tsv", check=False))
        lines = (d / "o0.tsv").read_text().splitlines()
        assert len(lines) == 1 and lines[0].split("\t")[:9] == [
            "predictor", "level", "model", "estimate", "std_error", "p_value", "p_adj", "n", "note"], lines


@suite.test
def test_stats_failures_are_errors():
    d = tmpdir()
    _stats_case(d)
    base = ["--master", d / "master.tsv", "--disorder", d / "dis.tsv", "--tree", d / "tree.nwk",
            "--out", d / "out.tsv"]
    must_fail(rscript("disorder_stats.R", *base, "--temp-cols", "niche_volcano", check=False),
              "temperature column")
    must_fail(rscript("disorder_stats.R", *base, "--fdr", "bogus", check=False), "unknown --fdr")
    m = (d / "master.tsv").read_text()
    write(d / "m_notier.tsv", m.replace("prometheus_tier", "tier"))
    must_fail(rscript("disorder_stats.R", "--master", d / "m_notier.tsv", *base[2:], check=False),
              "lacks 'prometheus_tier'")
    write(d / "m_nonbin.tsv", m.replace("\tTrue\tFalse\t0\t0\n", "\tTrue\tFalse\t2\t0\n", 1))
    assert (d / "m_nonbin.tsv").read_text() != m
    must_fail(rscript("disorder_stats.R", "--master", d / "m_nonbin.tsv", *base[2:], check=False),
              "not binary")
    lines = (d / "dis.tsv").read_text().splitlines(True)
    f = lines[1].rstrip("\n").split("\t")
    f[3] = "abc"                                   # frac_disordered_nonlcr = the response
    write(d / "dis_bad.tsv", "".join([lines[0], "\t".join(f) + "\n"] + lines[2:]))
    must_fail(rscript("disorder_stats.R", "--master", d / "master.tsv", "--disorder", d / "dis_bad.tsv",
                      *base[4:], check=False), "non-numeric")
    write(d / "t_nobl.nwk", "((G01,G02),(G03,G04));\n")
    must_fail(rscript("disorder_stats.R", *base[:4], "--tree", d / "t_nobl.nwk", "--out", d / "o.tsv",
                      check=False), "branch lengths")
    # a missing required package is an error (not a header-only output)
    empty_lib = tmpdir()
    e = {"R_LIBS": str(empty_lib), "R_LIBS_USER": str(empty_lib)}
    probe = sh([RSCRIPT, "-e", "cat(requireNamespace('phylolm', quietly=TRUE))"], extra_env=e)
    if probe.stdout.strip() == "FALSE":
        must_fail(sh([RSCRIPT, SCRIPTS_R / "disorder_stats.R", *base], extra_env=e),
                  "required R package 'phylolm'")
    else:
        print("    (phylolm is in the base R library here: missing-package check not simulable)")


# =========================================================================
# 8. Snakemake DAG (private copy of the repo; nothing is written to the repo)
# =========================================================================
@suite.test
def test_snakemake_dry_run_disorder_target():
    if not Path(SNAKEMAKE).is_file():
        print("    (snakemake not installed: dry run skipped)")
        return
    w = tmpdir("smk_")
    for sub in ("workflow", "config", ".test"):
        shutil.copytree(REPO / sub, w / sub,
                        ignore=shutil.ignore_patterns(".snakemake", "results", "logs", "__pycache__"))
    home = w / "home"
    home.mkdir()
    e = {"HOME": str(home), "XDG_CACHE_HOME": str(home / ".cache")}
    must_pass(sh([PYTHON, ".test/setup_dryrun.py"], extra_env=e, cwd=w))
    p = sh([SNAKEMAKE, "disorder", "-n", "-p", "--nolock", "-s", "workflow/Snakefile",
            "--configfile", ".test/config/config.yaml"], extra_env=e, cwd=w, timeout=600)
    must_pass(p)
    out = p.stdout + p.stderr
    for rule in ("disorder_predict", "disorder_map", "disorder_predict_bact", "disorder_stats",
                 "plot_disorder"):
        assert re.search(rf"^{rule}\s+1\s*$", out, re.M), f"{rule} missing from the DAG"
    assert '--tier-sets "A;A,B"' in out and "--metapredict-version V3 --device cpu" in out
    assert "--alignment .test/results/evolution/align/family.aln " in out
    REAL_RUNS.append("snakemake dry run")


if __name__ == "__main__":
    atexit.register(lambda: print(f"real-tool runs: {', '.join(REAL_RUNS) or 'none'}"))
    suite.main("DISORDER TESTS PASSED")
