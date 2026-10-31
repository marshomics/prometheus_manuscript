"""Methanogenesis-type statistics suite (design 6.2).   bash tests/run_leaf.sh methanogenesis_stats

Covers workflow/scripts/R/methanogenesis_stats.R with real phylolm / phytools / corHMM
runs on simulated data with known answers, the association_tests.R main guard it relies
on, workflow/envs/rstats.post-deploy.sh (the vendored corHMM 2.8 installed into
temporary libraries exactly as Snakemake runs it, plus every failure path), the
methanogenesis_stats rule of ecology.smk (a Snakemake mini-workflow) and the hand-off of
its table to results_registry.py. Known answers:
  * a predictor with many independent origins that tracks presence is detected
    (primary habitat-adjusted p < 0.05) and passes the replication guard;
  * a single-origin predictor perfectly co-distributed with presence is refused by the
    replication guard (p NA, 'replication guard' note, outside the FDR family);
  * corHMM::fitCorrelationTest's four-model AIC table, from the vendored package;
  * the events-per-variable truncation of the habitat set and the VIF, both checked
    against an independent Python re-implementation, and the VIF refusal;
  * every sensitivity row (unresolved-as-absent, complete-only, leave-clade-out, A+B);
  * failure modes (missing column / package / input, unrooted tree, malformed values,
    bad options) and the genuinely empty inputs (header-only, exit 0).
Prints METHANOGENESIS STATS TESTS PASSED only when every test passes.
"""
import atexit
import hashlib
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, run, read_tsv, REPO, SCRIPTS_PY, SCRIPTS_R,  # noqa: E402
                  PYTHON, RSCRIPT, env)

suite = Suite("methanogenesis_stats")
SCRIPT = SCRIPTS_R / "methanogenesis_stats.R"
ASSOC = SCRIPTS_R / "association_tests.R"
POST_DEPLOY = REPO / "workflow" / "envs" / "rstats.post-deploy.sh"
RSTATS_YAML = REPO / "workflow" / "envs" / "rstats.yaml"
TARBALL = REPO / "workflow" / "envs" / "vendor" / "corHMM_2.8.tar.gz"
CORHMM_SHA256 = "d6e46641f5b2b15720c8bb6c51d715220519dd8a6ba44171d13ef75bfa2eb6df"
SMK = os.environ.get("PROM_SNAKEMAKE") or shutil.which("snakemake") or \
    "/Users/jmarsh/.claude-science/conda/envs/prom-smk/bin/snakemake"

PREDICTORS = ["mgen_co2_reducing", "mgen_acetoclastic", "mgen_methylotrophic",
              "mgen_methyl_dismutating", "mgen_methyl_reducing"]
OUT_COLS = ["predictor", "analysis", "role", "estimate", "std_error", "ci_low", "ci_high",
            "p_value", "p_adj", "n", "n_present", "n_absent", "n_level1", "n_changes_mean",
            "n_changes_lo", "n_changes_hi", "replicated", "delta_aic", "best_model",
            "covariates", "vif", "model", "converged", "note"]
ANALYSES = {"marginal": "secondary", "habitat_adjusted": "primary", "simmap": "descriptive",
            "corhmm": "descriptive", "sens_unresolved_as_absent": "sensitivity",
            "sens_complete_only": "sensitivity", "sens_leave_clade_out": "sensitivity",
            "sens_tiers_AB": "sensitivity"}
ORDERS = ["o__Methanobacteriales", "o__Methanococcales", "o__Methanomicrobiales",
          "o__Methanosarcinales", "o__Methanomassiliicoccales", "o__Methanocellales"]
EPV_MIN, MIN_CLASS, MIN_N, COMPLETE = 10.0, 3, 10, 95.0


def _scratch_root():
    """Absolute scratch root (TMPDIR may be relative in this sandbox); removed at exit
    unless MGEN_KEEP_TMP=1."""
    for base in (os.environ.get("MGEN_TEST_TMP"), "/tmp", tempfile.gettempdir()):
        if base and os.path.isabs(base) and os.path.isdir(base) and os.access(base, os.W_OK):
            root = Path(tempfile.mkdtemp(prefix="mgen_stats_tests_", dir=base))
            if os.environ.get("MGEN_KEEP_TMP") != "1":
                atexit.register(shutil.rmtree, root, True)
            return root
    raise RuntimeError("no absolute writable scratch directory (set MGEN_TEST_TMP)")


SCRATCH = _scratch_root()


def tmpdir(prefix="t_"):
    return Path(tempfile.mkdtemp(prefix=prefix, dir=SCRATCH))


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def tsv(path, header, rows):
    return write(path, "\n".join(["\t".join(header)] +
                                 ["\t".join(str(r[c]) for c in header) for r in rows]) + "\n")


def fnum(x):
    return float("nan") if x in ("", "NA", None) else float(x)


def bh(p):
    m = len(p)
    order = sorted(range(m), key=lambda i: p[i])
    adj, prev = [0.0] * m, 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, p[i] * m / rank)
        adj[i] = prev
    return adj


# =============================================================================== R env
_LIB = {}


def r_env(extra_lib=None):
    e = dict(env())
    e["LANGUAGE"] = "en"
    if extra_lib:
        e["R_LIBS"] = os.pathsep.join([str(extra_lib)] + ([e["R_LIBS"]] if e.get("R_LIBS") else []))
    return e


def post_deploy(script, cwd, conda_prefix=None, extra=None, path_prefix=None):
    e = r_env()
    for k in ("CONDA_PREFIX", "PROM_CORHMM_TARBALL", "PROM_R_LIBRARY"):
        e.pop(k, None)
    if conda_prefix is not None:
        e["CONDA_PREFIX"] = str(conda_prefix)
    e["PATH"] = os.pathsep.join([p for p in (path_prefix, str(Path(RSCRIPT).parent)) if p] +
                                [e.get("PATH", "")])
    if extra:
        e.update(extra)
    return subprocess.run(["/bin/sh", str(script)], cwd=cwd, env=e, capture_output=True,
                          text=True, timeout=900)


def corhmm_lib():
    """The vendored corHMM 2.8, installed ONCE into a temporary library by the post-deploy
    script exactly as Snakemake runs it: a COPY of the script outside workflow/envs
    (Snakemake copies it to <env>.post-deploy.sh), the repository root as working
    directory and $CONDA_PREFIX of a fresh prefix."""
    if "lib" not in _LIB:
        root = tmpdir("corhmm_env_")
        prefix = root / "env"
        lib = prefix / "lib" / "R" / "library"
        lib.mkdir(parents=True)
        copy = root / "conda" / "0a1b2c3d_.post-deploy.sh"
        copy.parent.mkdir()
        shutil.copy(POST_DEPLOY, copy)
        p = post_deploy(copy, REPO, conda_prefix=prefix)
        assert p.returncode == 0, f"post-deploy failed ({p.returncode}):\n{p.stdout[-2000:]}\n{p.stderr[-3000:]}"
        _LIB.update(lib=lib, log=p.stdout + p.stderr)
    return _LIB["lib"]


def mstats(d, *args, master=None, tree=None, out="assoc.tsv", corhmm=False, check=True,
           extra_env=None, timeout=2400):
    a = ["--master", master or d / "master.tsv", "--tree", tree or d / "tree.nwk", "--out", d / out]
    if not corhmm:
        a += ["--corhmm", "FALSE"]
    e = r_env(corhmm_lib() if corhmm else None)
    if extra_env:
        e.update(extra_env)
    p = subprocess.run([RSCRIPT, str(SCRIPT), *map(str, a), *map(str, args)], cwd=REPO, env=e,
                       capture_output=True, text=True, timeout=timeout)
    if check and p.returncode != 0:
        raise AssertionError(f"methanogenesis_stats.R failed ({p.returncode}):\n{p.stdout[-1500:]}\n{p.stderr[-4000:]}")
    df = read_tsv(d / out) if p.returncode == 0 and (d / out).exists() else None
    return p, df


def rows_of(df):
    return {(r["predictor"], r["analysis"]): r for _, r in df.iterrows()}


# =============================================================================== simulation
def _coalescent(tips, rng, clades):
    """Random coalescent-shaped subtree; every internal node's tip set goes to `clades`."""
    items = [(t, 0.0, [t]) for t in tips]
    h = 0.0
    while len(items) > 1:
        k = len(items)
        h += rng.expovariate(k * (k - 1) / 2.0)
        i, j = sorted(rng.sample(range(k), 2))
        a, b = items[i], items[j]
        clades.append(a[2] + b[2])
        items = [x for n, x in enumerate(items) if n not in (i, j)] + \
            [(f"({a[0]}:{h - a[1]:.6f},{b[0]}:{h - b[1]:.6f})", h, a[2] + b[2])]
    return items[0]


def _join(parts, rng):
    items, h = list(parts), max(p[1] for p in parts)
    while len(items) > 1:
        h += rng.uniform(0.3, 0.8)
        i, j = sorted(rng.sample(range(len(items)), 2))
        a, b = items[i], items[j]
        items = [x for n, x in enumerate(items) if n not in (i, j)] + \
            [(f"({a[0]}:{h - a[1]:.6f},{b[0]}:{h - b[1]:.6f})", h, a[2] + b[2])]
    return items[0][0] + ";\n"


def _scattered(clades, rng, n_clades, max_size=4):
    """Tips of n_clades disjoint small clades (2..max_size tips) drawn across the whole
    tree: a trait with many independent origins but real phylogenetic structure."""
    pool = [c for c in clades if 2 <= len(c) <= max_size]
    rng.shuffle(pool)
    chosen, used = [], set()
    for c in pool:
        if len(chosen) == n_clades:
            break
        if not used & set(c):
            chosen.append(c)
            used |= set(c)
    assert len(chosen) == n_clades, (len(chosen), n_clades)
    return used


HABITAT = ["group_gut", "niche_digester", "niche_hot", "niche_marine", "niche_rare",
           "niche_sediment", "niche_twin"]
MASTER_COLS = (["genome_id", "representative", "tree_tip", "gtdb_domain", "gtdb_phylum",
                "gtdb_order", "gtdb_genus", "completeness", "contamination", "pass_primary",
                "genome_size", "mcrA", "methanogen", "is_anme", "methanogenesis_subtype",
                "prometheus_present", "prometheus_tier", "n_defense_systems_excl_prometheus"]
               + HABITAT + ["prometheus_status_primary", "prometheus_status_sensitivity"]
               + PREDICTORS + ["mgen_sub_h2", "mgen_sub_methanol", "subtype_status"])


def _status(tier, tiers):
    return "present" if tier in tiers else ("absent" if tier == "" else "excluded")


def build_fixture(kind, seed=11):
    """kind 'multi': mgen_co2_reducing has many independent origins (24 disjoint small
    clades scattered over the tree) and presence tracks it (90 % concordance); group_gut is a habitat confounder of it;
    niche_twin is a near copy of mgen_methylotrophic (VIF); mgen_acetoclastic is rare (2);
    mgen_methyl_dismutating has unresolved (NA) calls; mgen_methyl_reducing is concentrated
    in o__Methanomassiliicoccales; 6 B-only genomes. kind 'single': mgen_acetoclastic = 1
    for exactly one order and Prometheus is present exactly there (one origin, perfect
    co-distribution). Both carry non-methanogen / ANME / off-tree / non-pass_primary rows
    that must not enter the analysis set."""
    rng = random.Random(seed)
    d = tmpdir(f"mgen_{kind}_")
    sizes = [28, 24, 20, 16, 12, 10] if kind == "multi" else [24, 22, 20, 18, 14]
    groups, k = [], 0
    for oi, size in enumerate(sizes):
        groups.append((ORDERS[oi], [f"M{k + i:03d}" for i in range(size)]))
        k += size
    outg = [f"N{i:02d}" for i in range(10)] + ["A00", "A01", "Z00"]
    clades = []
    parts = [_coalescent(t, rng, clades) for _, t in groups]
    meth_clades = list(clades)
    write(d / "tree.nwk", _join(parts + [_coalescent(outg, rng, clades)], rng))
    rows = []

    def base_row(gid, tip, order, meth, anme=False):
        return dict(genome_id=gid, representative="True", tree_tip=tip, gtdb_domain="d__Archaea",
                    gtdb_phylum="p__Methanobacteriota" if meth else "p__Thermoproteota",
                    gtdb_order=order, gtdb_genus="g__X", completeness=f"{rng.uniform(85, 100):.2f}",
                    contamination=f"{rng.uniform(0, 3):.2f}", pass_primary="True",
                    genome_size=str(int(rng.uniform(1.5e6, 4.5e6))), mcrA=str(meth or anme),
                    methanogen=str(meth), is_anme=str(anme), methanogenesis_subtype="NA",
                    n_defense_systems_excl_prometheus=str(rng.randint(0, 8)),
                    mgen_sub_h2="NA", mgen_sub_methanol="NA", subtype_status="NA",
                    **{h: str(int(rng.random() < 0.4)) for h in HABITAT},
                    **{p: "NA" for p in PREDICTORS})

    meth_rows = []
    for order, tips in groups:
        for t in tips:
            r = base_row(t, t, order, True)
            rows.append(r)
            meth_rows.append(r)
    if kind == "multi":
        co2_on = _scattered(meth_clades, rng, 24)                    # 24 origins, ~half the tips
        methylo_on = _scattered(meth_clades, rng, 12)
        dism_on = _scattered(meth_clades, rng, 9)
        for r in meth_rows:
            co2 = int(r["genome_id"] in co2_on)
            y = co2 if rng.random() < 0.9 else 1 - co2
            methylo = int(r["genome_id"] in methylo_on)
            r.update(mgen_co2_reducing=str(co2), mgen_methylotrophic=str(methylo),
                     mgen_acetoclastic="0", mgen_methyl_dismutating=str(int(r["genome_id"] in dism_on)),
                     mgen_methyl_reducing="0", _y=y,
                     group_gut=str(co2 if rng.random() < 0.75 else 1 - co2),
                     niche_marine=str(int(rng.random() < 0.4)), niche_sediment=str(int(rng.random() < 0.3)),
                     niche_hot=str(int(rng.random() < 0.5)), niche_digester=str(int(rng.random() < 0.35)),
                     niche_rare="0", niche_twin=str(methylo))
        for r in rng.sample(meth_rows, 3):                     # niche_twin: 3 flips
            r["niche_twin"] = str(1 - int(r["niche_twin"]))
        for r in rng.sample(meth_rows, 2):
            r["niche_rare"] = "1"
        for r in rng.sample(meth_rows, 2):
            r["mgen_acetoclastic"] = "1"
        for r in rng.sample(meth_rows, 12):
            r["niche_digester"] = "NA"
        for r in rng.sample(meth_rows, 25):
            r["mgen_methyl_dismutating"] = "NA"
        for r in meth_rows:
            if r["gtdb_order"] == ORDERS[4]:
                r["mgen_methyl_reducing"] = "1"
        for r in rng.sample([r for r in meth_rows if r["gtdb_order"] != ORDERS[4]], 3):
            r["mgen_methyl_reducing"] = "1"
    else:
        for r in meth_rows:
            ac = int(r["gtdb_order"] == ORDERS[3])
            r.update(mgen_acetoclastic=str(ac), _y=ac, mgen_co2_reducing="1",
                     mgen_methylotrophic="0", mgen_methyl_dismutating="0", mgen_methyl_reducing="0",
                     niche_twin="NA", niche_rare="0")
    for r in meth_rows:
        r["prometheus_tier"] = "A" if r.pop("_y") else ""
    if kind == "multi":
        for r in rng.sample(meth_rows, 6):                     # B-only: excluded under {A}
            r["prometheus_tier"] = "B"
    # rows that must NOT enter the analysis set
    for t in outg:
        r = base_row(t, t, "o__Nitrososphaerales", t == "Z00", anme=t.startswith(("A", "Z")))
        r["prometheus_tier"] = "A" if rng.random() < 0.4 else ""
        if t == "Z00":
            r.update({p: "1" for p in PREDICTORS})             # methanogen flagged ANME
        rows.append(r)
    for i in range(3):                                         # methanogens not on the tree
        r = base_row(f"T{i:02d}", "", ORDERS[0], True)
        r.update({p: "1" for p in PREDICTORS}, prometheus_tier="A")
        rows.append(r)
    for r in [x for x in meth_rows if x["gtdb_order"] == ORDERS[1]][:2]:
        r["pass_primary"] = "False"
    for r in rows:
        r["prometheus_present"] = str(r["prometheus_tier"] != "")
        r["prometheus_status_primary"] = _status(r["prometheus_tier"], ("A",))
        r["prometheus_status_sensitivity"] = _status(r["prometheus_tier"], ("A", "B"))
    tsv(d / "master.tsv", MASTER_COLS, rows)
    return d, rows


_FIX = {}


def fixture(kind):
    if kind not in _FIX:
        _FIX[kind] = build_fixture(kind)
    return _FIX[kind]


def expected_set(rows, pred, tiers=("A",), x_na_as=None, keep=None):
    """Python re-implementation of the analysis set: methanogen tree tips (not ANME,
    pass_primary), tier rule, predictor NA dropped (or coded x_na_as)."""
    out = []
    for r in rows:
        if r["methanogen"] != "True" or r["is_anme"] == "True" or r["pass_primary"] != "True" \
                or not r["tree_tip"]:
            continue
        if keep is not None and not keep(r):
            continue
        t = r["prometheus_tier"]
        y = 1 if t in tiers else (0 if t == "" else None)
        if y is None:
            continue
        x = r[pred]
        if x == "NA":
            if x_na_as is None:
                continue
            x = str(x_na_as)
        out.append(dict(r, _x=int(x), _y=y))
    return out


def _pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    sa = math.sqrt(sum((v - ma) ** 2 for v in a))
    sb = math.sqrt(sum((v - mb) ** 2 for v in b))
    if sa == 0 or sb == 0:
        return float("nan")
    return sum((u - ma) * (v - mb) for u, v in zip(a, b)) / (sa * sb)


def py_select(rows, n_cov=2, epv=EPV_MIN, min_class=MIN_CLASS, min_n=MIN_N, habitat=HABITAT):
    """Independent re-implementation of the habitat selection (eligibility, |phi| order,
    EPV budget on the complete cases, class minima). Returns (included, rows, eligible)."""
    elig = []
    for h in habitat:
        ok = [r for r in rows if r[h] != "NA"]
        n1 = sum(r[h] == "1" for r in ok)
        if min(n1, len(ok) - n1) < min_class:
            continue
        phi = _pearson([r["_x"] for r in ok], [int(r[h]) for r in ok])
        if math.isfinite(phi):
            elig.append((h, phi))
    elig.sort(key=lambda t: (-abs(t[1]), t[0]))
    cur, inc = list(rows), []
    for h, _ in elig:
        nxt = [r for r in cur if r[h] != "NA"]
        npres = sum(r["_y"] for r in nxt)
        ne = min(npres, len(nxt) - npres)
        if 1 + n_cov + len(inc) + 1 > math.floor(ne / epv):
            continue
        n1 = sum(r["_x"] for r in nxt)
        h1 = sum(r[h] == "1" for r in nxt)
        if len(nxt) < min_n or ne < min_class or min(n1, len(nxt) - n1) < min_class or \
                min(h1, len(nxt) - h1) < min_class:
            continue
        vals = [r[h] for r in nxt]
        if any(vals == [r[g] for r in nxt] or vals == [str(1 - int(r[g])) for r in nxt] for g in inc):
            continue                                  # duplicate / complement: redundant
        cur, inc = nxt, inc + [h]
    return inc, cur, [h for h, _ in elig]


def py_vif(rows, terms):
    """VIF of x on (completeness, genome_size, habitat terms) by least squares; z-scoring
    is affine and leaves R^2 unchanged."""
    import numpy as np
    X = np.array([[1.0, float(r["completeness"]), float(r["genome_size"])] +
                  [float(r[h]) for h in terms] for r in rows])
    x = np.array([float(r["_x"]) for r in rows])
    beta, *_ = np.linalg.lstsq(X, x, rcond=None)
    res = x - X @ beta
    r2 = 1 - float(res @ res) / float(((x - x.mean()) ** 2).sum())
    return 1.0 / (1.0 - r2)


_RUNS = {}


def _warm():
    """The three expensive runs, computed once and concurrently (results do not depend on
    scheduling: every random step is seeded):
      main        all five predictors (the default list) on 'multi', corHMM off, threads 3
      single      the single-origin predictor on 'single', corHMM on (vendored package)
      corhmm_co2  the many-origin predictor on 'multi', corHMM on"""
    if _RUNS:
        return
    dm, rows_m = fixture("multi")
    ds, rows_s = fixture("single")
    corhmm_lib()                                            # install before forking runs
    jobs = {"main": (dm, rows_m, ("--threads", "3"), False, "main.tsv"),
            "single": (ds, rows_s, ("--predictors", "mgen_acetoclastic"), True, "single.tsv"),
            "corhmm_co2": (dm, rows_m, ("--predictors", "mgen_co2_reducing", "--nsim", "20"), True,
                           "corhmm_co2.tsv")}

    def one(item):
        name, (d, rows, args, ch, out) = item
        p, a = mstats(d, *args, corhmm=ch, out=out)
        return name, (d, rows, a, p.stderr)
    with ThreadPoolExecutor(max_workers=3) as pool:
        for name, res in pool.map(one, jobs.items()):
            _RUNS[name] = res


def main_run():
    _warm()
    return _RUNS["main"]


def single_run():
    _warm()
    return _RUNS["single"]


# =============================================================================== static
@suite.test
def test_static_parse_guard_and_hygiene():
    for f in (SCRIPT, ASSOC):
        p = run([RSCRIPT, "-e", f'invisible(parse("{f}")); cat("PARSE_OK")'])
        assert "PARSE_OK" in p.stdout, f
    # association_tests.R has a main guard: source() defines its helpers without running
    # main() (which would stop on the missing --master); so does methanogenesis_stats.R,
    # which locates association_tests.R next to itself even when sourced from elsewhere
    code = (f'source("{ASSOC}"); stopifnot(is.function(fit_phyloglm_ms), is.function(run_phyloglm), '
            f'is.function(read_species_tree), is.function(prepare_tree), is.function(tier_response)); '
            f'source("{SCRIPT}"); stopifnot(is.function(analyse_predictor), is.function(A$fit_phyloglm_ms)); '
            f'cat("GUARD_OK")')
    p = subprocess.run([RSCRIPT, "-e", code], cwd=str(SCRATCH), env=r_env(), capture_output=True, text=True)
    assert p.returncode == 0 and "GUARD_OK" in p.stdout, p.stderr[-1500:]
    assert "--master is required" not in p.stderr
    # helper known answers: Fitch-Hartigan counts (polytomy included), the radix (C-locale)
    # tie-break of the left-out order, VIF with two identical adjustment terms (aliased,
    # still finite) and with a copy of x (infinite), Mk Q assembly, strict indicators
    code = (f'source("{SCRIPT}"); '
            'tr <- ape::read.tree(text = "((a:1,b:1,c:1):1,(d:1,e:1):1);"); '
            'stopifnot(fitch_changes(tr, c(a = 1, b = 1, c = 0, d = 0, e = 1)) == 2L, '
            'fitch_changes(tr, c(a = 1, b = 1, c = 1, d = 1, e = 1)) == 0L, '
            'fitch_changes(ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);"), c(a = 0, b = 1, c = 0, d = 1)) == 2L); '
            'lo <- leave_out_order(data.frame(x = c(1, 1, 1, 1, 0), order = c("o__B", "o__a", "o__B", "o__a", "o__Z"))); '
            'stopifnot(lo$order == "o__B", lo$n == 2L, lo$tied); '
            'set.seed(1); d <- data.frame(x = rep(0:1, 10), h1 = rbinom(20, 1, 0.5), zc1 = rnorm(20)); d$h2 <- d$h1; '
            'v <- vif_of_x(d, c("zc1", "h1", "h2")); stopifnot(is.finite(v), abs(v - vif_of_x(d, c("zc1", "h1"))) < 1e-12); '
            'd$h3 <- d$x; stopifnot(is.infinite(vif_of_x(d, c("zc1", "h3")))); '
            'q <- mk_q(list(states = c("0", "1"), index.matrix = matrix(c(NA, 2, 1, NA), 2), rates = c(0.5, 2))); '
            'stopifnot(q["0", "1"] == 0.5, q["1", "0"] == 2, q["0", "0"] == -0.5); '
            'stopifnot(identical(parse_indicator(c("1", "0", "NA", "", "true"), "c"), c(1L, 0L, NA, NA, 1L))); '
            'e <- tryCatch(parse_indicator(c("1", "2"), "mgen_x"), error = function(e) conditionMessage(e)); '
            'stopifnot(grepl("mgen_x", e)); cat("HELPERS_OK")')
    p = subprocess.run([RSCRIPT, "-e", code], cwd=str(SCRATCH), env=r_env(), capture_output=True, text=True)
    assert p.returncode == 0 and "HELPERS_OK" in p.stdout, p.stderr[-1500:]
    run(["/bin/sh", "-n", str(POST_DEPLOY)])
    text = POST_DEPLOY.read_text()
    assert text.startswith("#!/bin/sh\n") and "|| true" not in text and "set -eu" in text
    code_only = re.sub(r"(?m)^\s*#.*$", "", text)
    bashisms = [b for b in (r"\[\[", r"\bfunction\b", r"\bsource\b", r"<<<", r"\$\{[A-Za-z_]+//",
                            r"\blocal\b", r"\bdeclare\b", r"\bpipefail\b") if re.search(b, code_only)]
    assert not bashisms, f"non-POSIX constructs in the post-deploy script: {bashisms}"
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"), str(REPO / "workflow/rules/ecology.smk"),
             str(REPO / "workflow/rules/synthesis.smk")])
    assert "RULES HYGIENE OK" in p.stdout, p.stdout
    text = (REPO / "workflow/rules/ecology.smk").read_text()
    rule = re.search(r"^rule methanogenesis_stats:\n(.*?)(?=^rule |\Z)", text, flags=re.S | re.M).group(1)
    for flag in ("--predictors", "--habitat", "--covariates", "--primary-tiers", "--sensitivity-tiers",
                 "--min-n", "--min-class", "--min-changes", "--nsim", "--corhmm", "--epv-min",
                 "--vif-max", "--complete-threshold", "--fdr-method", "--btol", "--seed", "--threads"):
        assert flag in rule, flag
    for key in ("predictors", "habitat", "covariates", "primary_tiers", "sensitivity_tiers", "min_n",
                "min_class", "min_changes", "nsim", "corhmm", "epv_min", "vif_max", "complete_threshold"):
        assert re.search(rf'_eco_mgen(_csv|_flag)?\("{key}"', rule), key
    for key in ('cfg("stats/fdr_method"', 'cfg("stats/phyloglm_btol"', 'cfg("seed"'):
        assert key in rule, key
    assert '"../envs/rstats.yaml"' in rule and "# nonfatal:" in rule and "> {log} 2>&1" in rule
    assert 'f"{ECO}/methanogenesis_associations.tsv"' in rule


@suite.test
def test_vendored_tarball_and_pins():
    """The vendored tarball has the recorded SHA-256 (script + README) and every corHMM
    dependency outside base R is pinned exactly in rstats.yaml."""
    assert hashlib.sha256(TARBALL.read_bytes()).hexdigest() == CORHMM_SHA256
    assert f"EXPECTED_SHA256={CORHMM_SHA256}" in POST_DEPLOY.read_text()
    assert CORHMM_SHA256 in (TARBALL.parent / "README.md").read_text()
    with tarfile.open(TARBALL) as tf:
        desc = tf.extractfile("corHMM/DESCRIPTION").read().decode()
    assert re.search(r"^Version: 2\.8$", desc, flags=re.M) and "NeedsCompilation: no" in desc
    deps = set()
    for key in ("Depends", "Imports"):
        m = re.search(rf"^{key}:(.*?)(?=^\S)", desc, flags=re.S | re.M)
        deps |= {re.sub(r"\(.*?\)", "", t).strip() for t in m.group(1).replace("\n", " ").split(",")}
    deps -= {"", "R", "parallel", "stats", "utils", "graphics", "grDevices", "methods"}
    pins = dict(re.findall(r"^\s+- (r-[a-z0-9.]+)=(\S+)", RSTATS_YAML.read_text(), flags=re.M))
    missing = sorted(t for t in deps if f"r-{t.lower()}" not in pins)
    assert not missing, f"corHMM dependencies not pinned in rstats.yaml: {missing}"
    for pkg in ("r-nloptr", "r-gensa", "r-expm", "r-numderiv", "r-corpcor", "r-viridis",
                "r-rmpfr", "r-igraph", "r-mass", "r-nnet"):
        assert re.fullmatch(r"[0-9][0-9._]*", pins[pkg]), (pkg, pins[pkg])
    assert "r-corhmm" not in pins


# =============================================================================== post-deploy
@suite.test
def test_post_deploy_installs_vendored_corhmm():
    lib = corhmm_lib()
    desc = (lib / "corHMM" / "DESCRIPTION").read_text()
    assert re.search(r"^Version: 2\.8$", desc, flags=re.M), desc[:300]
    log = _LIB["log"]
    assert "SHA-256 verified" in log and "corHMM 2.8 installed and loadable" in log, log[-1500:]
    assert "workflow/envs/vendor/corHMM_2.8.tar.gz" in log      # found from the repo root
    q = subprocess.run([RSCRIPT, "-e", 'cat(find.package("corHMM"), as.character(packageVersion("corHMM")), '
                        '"fitCorrelationTest" %in% getNamespaceExports("corHMM"))'],
                       env=r_env(lib), capture_output=True, text=True)
    got = q.stdout.split()
    assert q.returncode == 0 and Path(got[0]).resolve() == (lib / "corHMM").resolve() and \
        got[1:] == ["2.8", "TRUE"], (q.stdout, q.stderr[-800:])


@suite.test
def test_post_deploy_failures():
    root = tmpdir("pd_fail_")

    def layout(name, tar_bytes=None):
        """A copy of the script with its own vendor/ directory (in-place layout)."""
        d = root / name
        (d / "vendor").mkdir(parents=True)
        shutil.copy(POST_DEPLOY, d / "rstats.post-deploy.sh")
        if tar_bytes is not None:
            (d / "vendor" / "corHMM_2.8.tar.gz").write_bytes(tar_bytes)
        prefix = d / "env"
        (prefix / "lib" / "R" / "library").mkdir(parents=True)
        return d, prefix

    good = TARBALL.read_bytes()
    jobs = []
    d, pre = layout("sha", good[:-1] + bytes([good[-1] ^ 0xFF]))       # one flipped byte
    jobs.append(("SHA-256 mismatch", d, pre, None, None))
    d, pre = layout("missing")                                       # no tarball anywhere
    jobs.append(("not found", d, pre, None, None))
    d, pre = layout("noprefix", good)
    jobs.append(("CONDA_PREFIX is not set", d, None, None, None))
    d, pre = layout("nolib", good)
    shutil.rmtree(pre / "lib")
    jobs.append(("R library directory not found", d, pre, None, None))
    d, pre = layout("badr", good)
    fake = d / "fakebin"
    fake.mkdir()
    write(fake / "R", "#!/bin/sh\necho 'simulated R failure' >&2\nexit 3\n").chmod(0o755)
    jobs.append(("R CMD INSTALL of", d, pre, None, str(fake)))
    d, pre = layout("override", good[:1000])                          # the override is checked too
    jobs.append(("SHA-256 mismatch", d, pre, {"PROM_CORHMM_TARBALL": str(d / "vendor" / "corHMM_2.8.tar.gz")},
                 None))

    def one(job):
        msg, d, prefix, extra, pathp = job
        return msg, prefix, post_deploy(d / "rstats.post-deploy.sh", d, conda_prefix=prefix, extra=extra,
                                        path_prefix=pathp)

    with ThreadPoolExecutor(max_workers=3) as pool:
        for msg, prefix, p in pool.map(one, jobs):
            assert p.returncode != 0, (msg, p.stdout[-500:], p.stderr[-500:])
            assert msg in p.stderr, (msg, p.stderr[-800:])
            if prefix is not None and (prefix / "lib" / "R" / "library").exists():
                assert not (prefix / "lib" / "R" / "library" / "corHMM").exists(), msg


# =============================================================================== known answers
@suite.test
def test_output_layout_roles_and_analysis_set():
    d, rows, a, log = main_run()
    assert list(a.columns[:len(OUT_COLS)]) == OUT_COLS, list(a.columns)
    assert len(a) == len(PREDICTORS) * len(ANALYSES)
    assert [(r.predictor, r.analysis) for r in a.itertuples()] == \
        [(p, an) for p in PREDICTORS for an in ANALYSES]
    for r in a.itertuples():
        assert r.role == ANALYSES[r.analysis], (r.predictor, r.analysis, r.role)
    ix = rows_of(a)
    for p in PREDICTORS:
        want = expected_set(rows, p)
        m = ix[(p, "marginal")]
        assert int(m["n"]) == len(want), (p, m["n"], len(want))
        assert int(m["n_present"]) == sum(r["_y"] for r in want)
        assert int(m["n_level1"]) == sum(r["_x"] for r in want)
        assert m["n_excluded_tier"] == "6" and m["tier_set"] == "A", (m["n_excluded_tier"], m["tier_set"])
    n_unres = len(expected_set(rows, "mgen_methyl_dismutating", x_na_as=0)) - \
        len(expected_set(rows, "mgen_methyl_dismutating"))
    assert int(ix[("mgen_methyl_dismutating", "marginal")]["n_unresolved"]) == n_unres > 0
    assert "ANME dropped=1" in log and "not pass_primary dropped=2" in log, log[-2000:]
    # p_adj = BH over the primary rows with a p-value only
    prim = [(k, fnum(r["p_value"])) for k, r in ix.items()
            if r["role"] == "primary" and r["p_value"] != "NA"]
    assert len(prim) >= 2, prim
    for (k, _), q in zip(prim, bh([v for _, v in prim])):
        assert abs(fnum(ix[k]["p_adj"]) - q) < 1e-9, k
    for k, r in ix.items():
        if r["role"] != "primary" or r["p_value"] == "NA":
            assert r["p_adj"] == "NA", k
    for p in PREDICTORS:
        assert ix[(p, "corhmm")]["note"] == "corHMM test disabled (--corhmm FALSE)"


@suite.test
def test_multi_origin_association_detected_and_replicated():
    d, rows, a, _ = main_run()
    ix = rows_of(a)
    s = ix[("mgen_co2_reducing", "simmap")]
    assert s["replicated"] == "TRUE" and fnum(s["n_changes_mean"]) >= 5, dict(s)
    par = int(s["n_changes_parsimony"])
    assert 5 <= par <= 24, par                                   # 24 planted origins
    assert par <= fnum(s["n_changes_lo"]) <= fnum(s["n_changes_mean"]) <= fnum(s["n_changes_hi"])
    assert s["mk_model"] in ("ER", "ARD") and "stochastic maps, seed 12345" in s["note"]
    for an in ("marginal", "habitat_adjusted"):
        r = ix[("mgen_co2_reducing", an)]
        assert fnum(r["p_value"]) < 0.05 and fnum(r["estimate"]) > 1 and fnum(r["ci_low"]) > 0, (an, dict(r))
        assert r["converged"] == "TRUE" and r["replicated"] == "TRUE"
        assert "replication guard" not in r["note"]
    h = ix[("mgen_co2_reducing", "habitat_adjusted")]
    assert fnum(h["p_adj"]) < 0.05
    assert h["habitat_terms"].split(",")[0] == "group_gut", h["habitat_terms"]   # the confounder first
    assert h["covariates"].startswith("z(completeness),z(genome_size),group_gut")
    assert "z(completeness) + z(genome_size) + group_gut" in h["model"]
    # an unrelated many-origin predictor is not 'detected'
    assert fnum(ix[("mgen_methylotrophic", "marginal")]["p_value"]) > 0.01
    # few origins (one order + 3 tips): the replication guard withholds every p-value
    mr = ix[("mgen_methyl_reducing", "simmap")]
    if mr["replicated"] == "FALSE":
        for an in ("marginal", "habitat_adjusted", "sens_tiers_AB"):
            r = ix[("mgen_methyl_reducing", an)]
            assert r["p_value"] == "NA" and r["note"].startswith("replication guard:"), (an, r["note"][:120])


@suite.test
def test_epv_truncation_matches_python():
    d, rows, a, _ = main_run()
    ix = rows_of(a)
    for p in ("mgen_co2_reducing", "mgen_methyl_dismutating", "mgen_methyl_reducing"):
        want = expected_set(rows, p)
        inc, cur, elig = py_select(want)
        r = ix[(p, "habitat_adjusted")]
        got = [t for t in r["habitat_terms"].split(",") if t]
        assert got == inc, (p, got, inc)
        assert int(r["n_habitat_eligible"]) == len(elig), (p, r["n_habitat_eligible"], elig)
        assert "niche_rare" not in elig                         # 2 ones < min_class
        assert int(r["n"]) == len(cur)
        assert int(r["n_terms"]) <= max(3, math.floor(int(r["n_events"]) / EPV_MIN)), dict(r)
        assert int(r["n_terms"]) == 3 + len(inc)
        if len(elig) > len(inc):
            assert "EPV" in r["note"] and ("not added:" in r["note"] or "none of the" in r["note"]), r["note"]
    # stricter EPV: the budget floor(events/25) < x + 2 covariates -> no habitat term
    d2 = tmpdir("mgen_epv_")
    for f in ("master.tsv", "tree.nwk"):
        shutil.copy(d / f, d2 / f)
    _, b = mstats(d2, "--predictors", "mgen_co2_reducing", "--epv-min", "25", "--nsim", "10")
    r = rows_of(b)[("mgen_co2_reducing", "habitat_adjusted")]
    assert r["habitat_terms"] == "" and "already exceed" in r["note"], r["note"]
    assert r["p_value"] != "NA"                                 # still fitted, unadjusted for habitat
    # looser EPV admits more terms, again exactly the Python selection; a habitat column
    # duplicating an included one (niche_gutcopy = group_gut) is skipped as redundant
    m = read_tsv(d2 / "master.tsv")
    m["niche_gutcopy"] = m["group_gut"]
    m.to_csv(d2 / "master.tsv", sep="\t", index=False)
    rows2 = [dict(r, niche_gutcopy=r["group_gut"]) for r in rows]
    _, c = mstats(d2, "--predictors", "mgen_co2_reducing", "--epv-min", "6", "--nsim", "10", out="c.tsv")
    r = rows_of(c)[("mgen_co2_reducing", "habitat_adjusted")]
    inc, _, elig = py_select(expected_set(rows2, "mgen_co2_reducing"), epv=6,
                             habitat=HABITAT + ["niche_gutcopy"])
    got = [t for t in r["habitat_terms"].split(",") if t]
    assert got == inc and len(inc) >= 2 and "group_gut" in inc, (got, inc)
    assert "niche_gutcopy" in elig and "niche_gutcopy" not in got
    assert "niche_gutcopy: redundant with the adjustment terms already included" in r["note"], r["note"]


@suite.test
def test_vif_known_answer_and_refusal():
    d, rows, a, _ = main_run()
    ix = rows_of(a)
    for p in ("mgen_co2_reducing", "mgen_methylotrophic"):
        inc, cur, _ = py_select(expected_set(rows, p))
        r = ix[(p, "habitat_adjusted")]
        v = py_vif(cur, inc)
        assert abs(fnum(r["vif"]) - v) < 1e-6 * max(1.0, v), (p, r["vif"], v)
    m = ix[("mgen_methylotrophic", "habitat_adjusted")]
    assert m["habitat_terms"].split(",")[0] == "niche_twin", m["habitat_terms"]
    assert fnum(m["vif"]) > 5 and m["p_value"] == "NA" and m["estimate"] == "NA" and m["p_adj"] == "NA"
    assert "partial effect not identifiable" in m["note"], m["note"]
    assert ix[("mgen_methylotrophic", "marginal")]["p_value"] != "NA"   # marginal still fitted
    assert fnum(ix[("mgen_co2_reducing", "habitat_adjusted")]["vif"]) < 5


@suite.test
def test_sensitivity_rows():
    d, rows, a, _ = main_run()
    ix = rows_of(a)
    for p in ("mgen_co2_reducing", "mgen_methyl_dismutating"):
        for an in ("sens_unresolved_as_absent", "sens_complete_only", "sens_leave_clade_out",
                   "sens_tiers_AB"):
            r = ix[(p, an)]
            assert r["role"] == "sensitivity" and r["p_adj"] == "NA", (p, an)
            assert r["estimate"] != "NA", (p, an, r["note"])
    # unresolved coded absent: the unresolved genomes join with x = 0
    p = "mgen_methyl_dismutating"
    u, m = ix[(p, "sens_unresolved_as_absent")], ix[(p, "marginal")]
    inc, cur, _ = py_select(expected_set(rows, p, x_na_as=0))
    assert int(u["n"]) == len(cur) and u["habitat_terms"] == ",".join(inc), (u["n"], len(cur), u["habitat_terms"], inc)
    assert int(u["n_level1"]) == sum(r["_x"] for r in cur) and u["n_unresolved"] == "0"
    assert f"({m['n_unresolved']} genome(s))" in u["note"]
    co = ix[("mgen_co2_reducing", "sens_unresolved_as_absent")]
    assert "identical data to habitat_adjusted" in co["note"]
    assert co["estimate"] == ix[("mgen_co2_reducing", "habitat_adjusted")]["estimate"]
    # completeness >= 95 only
    for p in ("mgen_co2_reducing", "mgen_methyl_dismutating"):
        inc, cur, _ = py_select(expected_set(rows, p, keep=lambda r: float(r["completeness"]) >= COMPLETE))
        c = ix[(p, "sens_complete_only")]
        assert int(c["n"]) == len(cur) and c["habitat_terms"] == ",".join(inc), (p, c["n"], len(cur))
        assert "completeness >= 95 only" in c["note"]
    # leave out the gtdb_order with the most x = 1 tips (ties alphabetical)
    for p in ("mgen_co2_reducing", "mgen_methyl_reducing"):
        cnt = {}
        for r in expected_set(rows, p):
            if r["_x"] == 1:
                cnt[r["gtdb_order"]] = cnt.get(r["gtdb_order"], 0) + 1
        top = sorted(o for o, v in cnt.items() if v == max(cnt.values()))[0]
        lo = ix[(p, "sens_leave_clade_out")]
        assert f"left out gtdb_order {top} " in lo["note"], (p, top, lo["note"])
        inc, cur, _ = py_select(expected_set(rows, p, keep=lambda r, top=top: r["gtdb_order"] != top))
        assert int(lo["n"]) == len(cur), (p, lo["n"], len(cur))
    assert "left out gtdb_order o__Methanomassiliicoccales" in \
        ix[("mgen_methyl_reducing", "sens_leave_clade_out")]["note"]
    # response by A+B: the 6 B-only genomes are present, none excluded
    for p in ("mgen_co2_reducing", "mgen_methyl_dismutating"):
        t = ix[(p, "sens_tiers_AB")]
        want = expected_set(rows, p, tiers=("A", "B"))
        inc, cur, _ = py_select(want)
        assert t["tier_set"] == "A,B" and t["n_excluded_tier"] == "0"
        assert int(t["n"]) == len(cur) and int(t["n_present"]) == sum(r["_y"] for r in cur), (p, t["n"], len(cur))
        assert sum(1 for r in want if r["prometheus_tier"] == "B") > 0


@suite.test
def test_rare_predictor_not_fitted():
    _, _, a, _ = main_run()
    ix = rows_of(a)
    for an in ANALYSES:
        r = ix[("mgen_acetoclastic", an)]
        assert r["p_value"] == "NA" and r["estimate"] == "NA", an
    for an in ("marginal", "habitat_adjusted", "simmap"):
        assert "predictor 1=2" in ix[("mgen_acetoclastic", an)]["note"], an
    assert ix[("mgen_acetoclastic", "simmap")]["replicated"] == "NA"


@suite.test
def test_single_origin_refused_by_replication_guard():
    d, rows, a, _ = single_run()
    ix = rows_of(a)
    s = ix[("mgen_acetoclastic", "simmap")]
    assert s["replicated"] == "FALSE" and fnum(s["n_changes_mean"]) < 5, dict(s)
    assert s["n_changes_parsimony"] == "1"                     # one origin
    for an in ("marginal", "habitat_adjusted", "sens_complete_only", "sens_tiers_AB"):
        r = ix[("mgen_acetoclastic", an)]
        assert r["p_value"] == "NA" and r["p_adj"] == "NA", (an, r["p_value"])
        assert r["note"].startswith("replication guard:"), (an, r["note"][:200])
        assert math.isfinite(fnum(r["estimate"])) and math.isfinite(fnum(r["ci_low"])), (an, dict(r))
        assert r["replicated"] == "FALSE"
    assert "quasi-separation" in ix[("mgen_acetoclastic", "marginal")]["note"]
    lo = ix[("mgen_acetoclastic", "sens_leave_clade_out")]
    assert f"left out gtdb_order {ORDERS[3]}" in lo["note"] and "predictor 1=0" in lo["note"]
    assert all(r["p_adj"] == "NA" for r in ix.values())         # nothing enters the FDR family


@suite.test
def test_corhmm_four_model_table():
    """corHMM::fitCorrelationTest from the vendored 2.8: four AICs / weights, delta AIC =
    min(independent) - min(correlated), best = argmin, no p-value. A single origin does
    not support correlated evolution (delta AIC < 0); the many-origin association does
    (delta AIC > 2, a correlated model best)."""
    models = ["independent", "correlated", "hmm_independent", "hmm_correlated"]
    checks = ((single_run()[2], "mgen_acetoclastic", lambda v: v < 0),
              (_RUNS["corhmm_co2"][2], "mgen_co2_reducing", lambda v: v > 2))
    for frame, pred, ok in checks:
        r = rows_of(frame)[(pred, "corhmm")]
        aic = {m: fnum(r[f"aic_{m}"]) for m in models}
        w = {m: fnum(r[f"aicw_{m}"]) for m in models}
        assert all(math.isfinite(v) for v in aic.values()), dict(r)
        assert abs(sum(w.values()) - 1) < 1e-9
        delta = min(aic["independent"], aic["hmm_independent"]) - min(aic["correlated"], aic["hmm_correlated"])
        assert abs(fnum(r["delta_aic"]) - delta) < 1e-6
        assert r["best_model"] == min(aic, key=aic.get) and r["p_value"] == "NA" and r["estimate"] == "NA"
        assert "no p-value" in r["note"] and r["role"] == "descriptive"
        for m in models:
            assert f"{m}: logL=" in r["note"], m
        assert ok(fnum(r["delta_aic"])), (pred, r["delta_aic"], r["best_model"])
    assert rows_of(_RUNS["corhmm_co2"][2])[("mgen_co2_reducing", "corhmm")]["best_model"] in \
        ("correlated", "hmm_correlated")


@suite.test
def test_reproducible_across_threads_and_seeded():
    d, _ = fixture("multi")
    args = ("--predictors", "mgen_co2_reducing,mgen_methylotrophic", "--nsim", "20")
    jobs = [(("--threads", "1"), "t1.tsv"), (("--threads", "2"), "t2.tsv"), (("--seed", "99"), "s99.tsv")]
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(lambda j: mstats(d, *args, *j[0], out=j[1]), jobs))
    assert (d / "t1.tsv").read_text() == (d / "t2.tsv").read_text(), "results depend on --threads"
    a, b = rows_of(read_tsv(d / "t1.tsv")), rows_of(read_tsv(d / "s99.tsv"))
    assert "seed 99" in b[("mgen_co2_reducing", "simmap")]["note"]
    assert a[("mgen_co2_reducing", "simmap")]["n_changes_mean"] != b[("mgen_co2_reducing", "simmap")]["n_changes_mean"]


@suite.test
def test_registry_handoff():
    """The real producer table enters results_registry.py as the exploratory
    methanogenesis family; its BH (primary rows only) equals the producer's p_adj."""
    d, _, a, _ = main_run()
    reg_dir = tmpdir("mgen_reg_")
    headers = {"association": "test\tmodel\tterm\testimate\tstd_error\tstatistic\tp_value\tci_low\tci_high\tn\tnote",
               "cooccurrence": "system_a\tsystem_b\tn_a\tn_b\tn_both\tmethod\tstatistic\tp_value\tp_adj\tdirection",
               "ecology": "predictor\tlevel\tmodel\testimate\tstd_error\tp_value\tp_adj\tn\tnote",
               "selection": "method\tstatistic\tvalue\tp_value\tn_sites_or_branches\tnote",
               "reconciliation": "tool\tevent\tcount\tfraction\tnote", "ancestral": "metric\tvalue\tnote"}
    args = []
    for fam, h in headers.items():
        args += [f"--{fam}", write(reg_dir / f"{fam}.tsv", h + "\n")]
    out = reg_dir / "registry.tsv"
    run([PYTHON, str(SCRIPTS_PY / "results_registry.py"), *map(str, args), "--methanogenesis",
         str(d / "main.tsv"), "--primary-tiers", "A", "--strict-confirmatory", "--out", str(out)])
    reg = read_tsv(out)
    fam = reg[reg.family == "methanogenesis"]
    assert len(fam) == len(a) and set(fam.hypothesis) == {"exploratory"}
    ix = {(r.term, r.test): r for r in fam.itertuples()}
    src = rows_of(a)
    m = sum(1 for r in src.values() if r["role"] == "primary" and r["p_value"] != "NA")
    for (pred, an), r in src.items():
        g = ix[(pred, an)]
        assert g.role == r["role"], (pred, an, g.role)
        if r["p_adj"] != "NA":
            assert abs(float(g.p_adj) - fnum(r["p_adj"])) < 1e-9 and g.fdr_m == str(m), (pred, an)
        else:
            assert g.p_adj == "" and g.fdr_m == "", (pred, an, g.p_adj)
    s = ix[("mgen_co2_reducing", "simmap")]
    assert s.estimate_scale == "count" and \
        abs(float(s.estimate) - fnum(src[("mgen_co2_reducing", "simmap")]["n_changes_mean"])) < 1e-9
    assert ix[("mgen_co2_reducing", "habitat_adjusted")].estimate_scale == "log_odds"
    assert ix[("mgen_co2_reducing", "sens_tiers_AB")].tier_set == "A,B"
    cor = rows_of(_RUNS["corhmm_co2"][2])
    out2 = reg_dir / "registry2.tsv"
    run([PYTHON, str(SCRIPTS_PY / "results_registry.py"), *map(str, args), "--methanogenesis",
         str(_RUNS["corhmm_co2"][0] / "corhmm_co2.tsv"), "--primary-tiers", "A", "--out", str(out2)])
    g = {(r.term, r.test): r for r in read_tsv(out2).itertuples()}[("mgen_co2_reducing", "corhmm")]
    assert g.estimate_scale == "statistic" and abs(float(g.statistic) - fnum(cor[("mgen_co2_reducing", "corhmm")]["delta_aic"])) < 1e-9


# =============================================================================== failures
@suite.test
def test_failure_modes():
    d, rows = fixture("multi")
    f = tmpdir("mgen_fail_")
    m = read_tsv(d / "master.tsv")
    m.drop(columns=["mgen_methyl_reducing"]).to_csv(f / "nocol.tsv", sep="\t", index=False)
    m.drop(columns=["gtdb_order"]).to_csv(f / "noorder.tsv", sep="\t", index=False)
    bad = m.copy()
    bad.loc[bad.index[0], "mgen_co2_reducing"] = "2"
    bad.to_csv(f / "badval.tsv", sep="\t", index=False)
    badh = m.copy()
    badh.loc[badh.index[0], "niche_hot"] = "maybe"
    badh.to_csv(f / "badhab.tsv", sep="\t", index=False)
    dup = m.copy()
    dup.loc[dup.index[1], "genome_id"] = dup.loc[dup.index[0], "genome_id"]
    dup.to_csv(f / "dup.tsv", sep="\t", index=False)
    badtier = m.copy()
    badtier.loc[badtier.index[0], "prometheus_tier"] = "Rejected"
    badtier.to_csv(f / "badtier.tsv", sep="\t", index=False)
    write(f / "zero.tsv", "")
    unrooted = run([RSCRIPT, "-e", f'cat(ape::write.tree(ape::unroot(ape::read.tree("{d / "tree.nwk"}"))))']).stdout
    write(f / "unrooted.nwk", unrooted + "\n")
    write(f / "nobl.nwk", re.sub(r":[0-9.]+", "", (d / "tree.nwk").read_text()))
    cases = [
        ("lacks required column(s): mgen_methyl_reducing", {"--master": f / "nocol.tsv", "--predictors": ",".join(PREDICTORS)}),
        ("lacks required column(s): gtdb_order", {"--master": f / "noorder.tsv"}),
        ("not in {1, 0, NA}", {"--master": f / "badval.tsv"}),
        ("column 'niche_hot'", {"--master": f / "badhab.tsv"}),
        ("duplicated genome_id", {"--master": f / "dup.tsv"}),
        ("prometheus_tier", {"--master": f / "badtier.tsv"}),
        ("0-byte", {"--master": f / "zero.tsv"}),
        ("master table not found", {"--master": f / "nope.tsv"}),
        ("not rooted", {"--tree": f / "unrooted.nwk"}),
        ("branch lengths", {"--tree": f / "nobl.nwk"}),
        ("--min-changes must be > 0", {"--min-changes": "0"}),
        ("--nsim must be >= 1", {"--nsim": "0"}),
        ("is not a p.adjust method", {"--fdr-method": "qvalue"}),
        ("--corhmm must be TRUE or FALSE", {"--corhmm": "maybe"}),
        ("--primary-tiers", {"--primary-tiers": "Z"}),
        ("--habitat column(s) not in the master table", {"--habitat": "niche_nope"}),
        ("lacks required column(s): nope_cov", {"--covariates": "completeness,nope_cov"}),
    ]

    def one(job):
        i, (msg, kv) = job
        args = {"--master": d / "master.tsv", "--tree": d / "tree.nwk", "--corhmm": "FALSE",
                "--predictors": "mgen_co2_reducing", "--nsim": "5"}
        args.update(kv)
        out = f / f"x{i}.tsv"
        flat = [str(x) for kv2 in args.items() for x in kv2]
        p = subprocess.run([RSCRIPT, str(SCRIPT), *flat, "--out", str(out)], cwd=REPO, env=r_env(),
                           capture_output=True, text=True, timeout=600)
        return msg, p, out.exists()
    with ThreadPoolExecutor(max_workers=3) as pool:
        for msg, p, made in pool.map(one, list(enumerate(cases))):
            assert p.returncode != 0 and not made, (msg, p.returncode, made, p.stderr[-400:])
            assert msg in p.stderr, (msg, p.stderr[-600:])
    # --corhmm TRUE without the corHMM package is an error (vendored lib not on R_LIBS)
    probe = subprocess.run([RSCRIPT, "-e", 'cat(requireNamespace("corHMM", quietly=TRUE))'],
                           env=r_env(), capture_output=True, text=True)
    if probe.stdout.strip() == "FALSE":
        p = subprocess.run([RSCRIPT, str(SCRIPT), "--master", str(d / "master.tsv"), "--tree",
                            str(d / "tree.nwk"), "--predictors", "mgen_co2_reducing", "--corhmm", "TRUE",
                            "--out", str(f / "nopkg.tsv")], cwd=REPO, env=r_env(), capture_output=True,
                           text=True, timeout=600)
        assert p.returncode != 0 and "corHMM" in p.stderr and not (f / "nopkg.tsv").exists(), p.stderr[-600:]
    else:
        print("    note: corHMM is installed in the base R library; missing-corHMM case not run")
    # a missing core package (phylolm hidden) is an error, not a skipped test
    e = r_env()
    e["R_LIBS"] = ""
    e["R_LIBS_USER"] = "/nonexistent"
    q = subprocess.run([RSCRIPT, "-e", 'cat(requireNamespace("phylolm", quietly=TRUE))'], env=e,
                       capture_output=True, text=True)
    if q.stdout.strip() == "FALSE":
        p = subprocess.run([RSCRIPT, str(SCRIPT), "--master", str(d / "master.tsv"), "--tree",
                            str(d / "tree.nwk"), "--corhmm", "FALSE", "--out", str(f / "z.tsv")],
                           env=e, capture_output=True, text=True, cwd=REPO)
        assert p.returncode != 0 and not (f / "z.tsv").exists(), p.stderr[-400:]


@suite.test
def test_empty_inputs_header_only():
    d, rows = fixture("multi")
    f = tmpdir("mgen_empty_")
    write(f / "empty.nwk", "();\n")
    m = read_tsv(d / "master.tsv")
    m.iloc[0:0].to_csv(f / "norows.tsv", sep="\t", index=False)
    nm = m.copy()
    nm["methanogen"] = "False"
    nm.to_csv(f / "nometh.tsv", sep="\t", index=False)
    for name, kw, why in (("t.tsv", {"tree": f / "empty.nwk", "master": d / "master.tsv"}, "species tree is empty"),
                          ("r.tsv", {"master": f / "norows.tsv", "tree": d / "tree.nwk"}, "master table has zero genomes"),
                          ("m.tsv", {"master": f / "nometh.tsv", "tree": d / "tree.nwk"}, "no methanogen tip")):
        p, a = mstats(f, out=name, **kw)
        assert a is not None and a.empty and list(a.columns[:len(OUT_COLS)]) == OUT_COLS, name
        assert why in p.stderr, (name, p.stderr[-400:])


# =============================================================================== Snakemake rule
MINI_SNAKEFILE = r"""
import os
configfile: "config.yaml"
RESULTS, LOGS = "results", "logs"
GENOMES = []
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
    input: f"{RESULTS}/ecology/methanogenesis_associations.tsv"
rule build_master_table:
    input: "master_src.tsv"
    output: master=f"{RESULTS}/scaffolds/master_table.tsv"
    shell: "cp {input} {output.master}"
rule species_tree:
    input: "tree_src.nwk"
    output: tree=f"{RESULTS}/scaffolds/species_tree.treefile"
    shell: "cp {input} {output.tree}"
include: "@ECOLOGY_SMK@"
"""


def _mini_workflow(config_extra=""):
    d, _ = fixture("multi")
    w = tmpdir("mgen_smk_")
    (w / "workflow").mkdir()
    os.symlink(REPO / "workflow" / "scripts", w / "workflow" / "scripts")
    shutil.copy(d / "master.tsv", w / "master_src.tsv")
    shutil.copy(d / "tree.nwk", w / "tree_src.nwk")
    write(w / "samples.tsv", "genome_id\n")
    write(w / "config.yaml", "samples: samples.tsv\nseed: 7\nstats:\n  confidence_tiers_primary: [A]\n"
          "  confidence_tiers_sensitivity: [A, B]\n  fdr_method: BH\n" + config_extra)
    # the real ecology.smk with only its conda: directives removed (Snakemake 8.30 touches
    # the conda wrapper of a job even without --use-conda); shell/params/log unchanged
    src = (REPO / "workflow/rules/ecology.smk").read_text()
    stripped, n_conda = re.subn(r'\n[ ]+conda:[ ]*\n[ ]+"\.\./envs/[a-z_]+\.yaml"[ ]*(?=\n)', "", src)
    assert n_conda == 7, f"expected 7 conda directives in ecology.smk, found {n_conda}"
    write(w / "workflow" / "rules" / "ecology.smk", stripped)
    write(w / "Snakefile", MINI_SNAKEFILE.replace("@ECOLOGY_SMK@", str(w / "workflow/rules/ecology.smk")))
    return w


def _smk(w, *extra):
    e = r_env()
    e["PATH"] = os.pathsep.join([str(Path(PYTHON).parent), str(Path(RSCRIPT).parent), e["PATH"]])
    e["HOME"] = str(w / "home")
    e["XDG_CACHE_HOME"] = str(w / "home" / ".cache")
    cmd = [SMK, "all", "-s", str(w / "Snakefile"), "-d", str(w), "--cores", "2", "--nolock",
           "--drop-metadata", *extra]
    return subprocess.run(cmd, cwd=w, env=e, text=True, capture_output=True, timeout=1800)


@suite.test
def test_snakemake_rule_end_to_end():
    assert os.path.exists(SMK) or shutil.which("snakemake"), f"snakemake not found ({SMK}); set PROM_SNAKEMAKE"
    w = _mini_workflow("methanogenesis:\n  stats:\n    predictors: [mgen_co2_reducing, mgen_methylotrophic]\n"
                       "    nsim: 20\n    corhmm: false\n    min_changes: 4\n    complete_threshold: 90\n"
                       "    habitat: [group_gut, niche_hot]\n")
    dry = _smk(w, "-n", "-p")
    assert dry.returncode == 0, dry.stderr[-3000:]
    out = dry.stdout + dry.stderr
    for frag in ("--predictors 'mgen_co2_reducing,mgen_methylotrophic'", "--habitat 'group_gut,niche_hot'",
                 "--covariates 'completeness,genome_size'", "--primary-tiers A --sensitivity-tiers A,B",
                 "--min-n 10 --min-class 3", "--min-changes 4.0 --nsim 20 --corhmm FALSE",
                 "--epv-min 10.0 --vif-max 5.0", "--complete-threshold 90.0",
                 "--fdr-method BH --btol 10 --seed 7", "logs/ecology/methanogenesis_stats.log"):
        assert frag in out, (frag, out[-2500:])
    p = _smk(w)
    assert p.returncode == 0, p.stderr[-4000:]
    a = read_tsv(w / "results/ecology/methanogenesis_associations.tsv")
    assert sorted(set(a.predictor)) == ["mgen_co2_reducing", "mgen_methylotrophic"]
    ix = rows_of(a)
    s = ix[("mgen_co2_reducing", "simmap")]
    assert "nsim=20" in s["model"] and "min_changes=4" in s["note"] and "seed 7" in s["note"]
    assert ix[("mgen_co2_reducing", "corhmm")]["note"] == "corHMM test disabled (--corhmm FALSE)"
    assert "completeness >= 90 only" in ix[("mgen_co2_reducing", "sens_complete_only")]["note"]
    terms = ix[("mgen_co2_reducing", "habitat_adjusted")]["habitat_terms"]
    assert set(t for t in terms.split(",") if t) <= {"group_gut", "niche_hot"}, terms
    log = (w / "logs/ecology/methanogenesis_stats.log").read_text()
    assert "[methanogenesis_stats]" in log and "row(s) for 2 predictor(s)" in log, log[-1500:]
    # a missing required column fails the rule loudly (no output)
    m = read_tsv(w / "master_src.tsv")
    m.drop(columns=["mgen_methylotrophic"]).to_csv(w / "master_src.tsv", sep="\t", index=False)
    shutil.rmtree(w / "results")
    p = _smk(w)
    assert p.returncode != 0 and not (w / "results/ecology/methanogenesis_associations.tsv").exists()
    assert "lacks required column(s): mgen_methylotrophic" in \
        (w / "logs/ecology/methanogenesis_stats.log").read_text()


if __name__ == "__main__":
    suite.main("METHANOGENESIS STATS TESTS PASSED")
