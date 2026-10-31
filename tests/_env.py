"""Shared helpers for the per-stage synthetic test suites (tests/<stage>/run_tests.py).

Tool locations come from environment variables so the suites stay portable:
  PROM_PYTHON     python with pandas/numpy/scipy/matplotlib/dendropy/biopython (default: this interpreter)
  PROM_TOOLS_BIN  directory holding hmmer/mafft/famsa/clipkit/iqtree2/mmseqs/prodigal/pal2nal/hyphy (prepended to PATH)
  PROM_RSCRIPT    Rscript with ape/phytools/phylolm/caper/ggtree/optparse (default: Rscript on PATH)
  R_LIBS          extra R library path(s), passed through unchanged
tests/run_leaf.sh sources an optional untracked `.unlazy/toolchain.env` that sets these.
"""
from __future__ import annotations
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS_PY = REPO / "workflow" / "scripts" / "py"
SCRIPTS_R = REPO / "workflow" / "scripts" / "R"
FIXTURES = REPO / "tests" / "fixtures"
PYTHON = os.environ.get("PROM_PYTHON") or sys.executable
TOOLS_BIN = os.environ.get("PROM_TOOLS_BIN", "")
RSCRIPT = os.environ.get("PROM_RSCRIPT") or shutil.which("Rscript") or "Rscript"


def env() -> dict:
    e = dict(os.environ)
    if TOOLS_BIN:
        e["PATH"] = TOOLS_BIN + os.pathsep + e.get("PATH", "")
    e.setdefault("MPLBACKEND", "Agg")
    return e


def have(tool: str) -> bool:
    return shutil.which(tool, path=env()["PATH"]) is not None


def run(cmd, cwd=None, check=True, timeout=1800):
    """Run a command (list or shell string) with the toolchain PATH; return CompletedProcess.
    On failure with check=True raise AssertionError carrying the tail of stdout/stderr."""
    shell = isinstance(cmd, str)
    p = subprocess.run(cmd, cwd=cwd or REPO, env=env(), shell=shell, text=True,
                       capture_output=True, timeout=timeout,
                       executable="/bin/bash" if shell else None)
    if check and p.returncode != 0:
        raise AssertionError(f"command failed ({p.returncode}): {cmd}\n--- stdout ---\n"
                             f"{p.stdout[-3000:]}\n--- stderr ---\n{p.stderr[-3000:]}")
    return p


def py(script: str, *args, **kw):
    return run([PYTHON, str(SCRIPTS_PY / script), *map(str, args)], **kw)


def rscript(script: str, *args, **kw):
    return run([RSCRIPT, str(SCRIPTS_R / script), *map(str, args)], **kw)


def _scratch_base() -> str:
    """An ABSOLUTE writable base for test scratch. TMPDIR may be relative (some sandboxes
    set './.tmp'), and tempfile/abspath then give paths that break as soon as a
    subprocess runs with another working directory. PROM_TEST_TMP overrides."""
    for base in (os.environ.get("PROM_TEST_TMP"), os.environ.get("TMPDIR"), "/tmp",
                 tempfile.gettempdir()):
        if base and os.path.isabs(base) and os.path.isdir(base) and os.access(base, os.W_OK):
            return base
    raise RuntimeError("no absolute writable scratch directory (set PROM_TEST_TMP)")


def tmpdir(prefix="promtest_") -> Path:
    """A fresh absolute scratch dir, removed at interpreter exit unless PROM_KEEP_TMP=1."""
    d = Path(tempfile.mkdtemp(prefix=prefix, dir=_scratch_base()))
    if os.environ.get("PROM_KEEP_TMP") != "1":
        import atexit
        import shutil
        atexit.register(shutil.rmtree, str(d), True)
    return d


def read_tsv(path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None)


class Suite:
    """Minimal test collector: register functions with @suite.test, call suite.main(token)."""
    def __init__(self, name: str):
        self.name, self.tests = name, []

    def test(self, fn):
        self.tests.append(fn)
        return fn

    def main(self, success_token: str):
        failed = 0
        for fn in self.tests:
            try:
                fn()
                print(f"  ok    {fn.__name__}")
            except Exception as ex:  # noqa: BLE001 - report every failure
                failed += 1
                print(f"  FAIL  {fn.__name__}: {type(ex).__name__}: {str(ex)[:2000]}")
        print(f"{self.name}: {len(self.tests) - failed}/{len(self.tests)} passed")
        if failed or not self.tests:
            sys.exit(1)
        print(success_token)
