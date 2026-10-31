"""Stage 7 (synthesis) synthetic test suite.

Covers results_registry.py, plot_summary_forest.py, plot_tree_master.R and the
synthesis.smk rules on synthetic inputs with known answers:
  * the FDR code against R's stats::p.adjust (real R) and hand-computed BH values;
  * exact confirmatory tagging against the labels the real producers write;
  * tier-set sensitivity rows, duplicates, strict label-mismatch failure, empty and
    malformed inputs;
  * the always-exploratory methanogenesis family (methanogenesis_stats.R): producer roles
    kept, BH over its primary rows only, count / statistic scales for the stochastic-map
    and corHMM rows, role / p-value consistency failures, tier-set demotion;
  * the per-panel forest (scales, null lines, not-estimable rows, placeholders);
  * the master tree with real ggtree/gheatmap/ggnewscale (tip matching, quoted labels,
    'no data' vs absence, ANME, B-only genome kept as tier B, failure modes);
  * a Snakemake dry run of the synthesis rules in a private copy of the repository.
Run:  bash tests/run_leaf.sh synthesis
"""
from __future__ import annotations

import io
import os
from concurrent.futures import ThreadPoolExecutor
import random
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, rscript, run, tmpdir, read_tsv,  # noqa: E402
                  REPO, SCRIPTS_PY, SCRIPTS_R, PYTHON, RSCRIPT)

# Private, writable matplotlib cache for every subprocess (the toolchain HOME may be
# read-only for this process); no .pyc files written into the repository.
_MPL = tempfile.mkdtemp(prefix="syn_mpl_")
os.environ["MPLCONFIGDIR"] = _MPL
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True
sys.path.insert(0, str(SCRIPTS_PY))

import results_registry as rr  # noqa: E402

S = Suite("synthesis")
OWNED_PY = [SCRIPTS_PY / "results_registry.py", SCRIPTS_PY / "plot_summary_forest.py",
            HERE / "run_tests.py"]
OWNED_R = [SCRIPTS_R / "plot_tree_master.R"]
SCHEMA_COLS = ["family", "test", "term", "estimate", "ci_low", "ci_high", "p_value",
               "p_adj", "hypothesis"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def write(path: Path, header: str, rows=()):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    return path


def close(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


ASSOC_H = "test\tmodel\tterm\testimate\tstd_error\tstatistic\tp_value\tci_low\tci_high\tn\tnote"
ASSOC_ROWS = [  # labels exactly as association_tests.R (revised, six rows) writes them
    ("naive_glm", "glm(present~methanogen, binomial); LRT p, profile CI", "methanogen", 1.5,
     0.5, 9.1, 0.0027, 0.6, 2.6, 40, "assumes independence; baseline only; log-odds"),
    ("naive_fisher", "fisher.test(methanogen x present); conditional MLE odds ratio",
     "methanogen", 1.4, "NA", "NA", 0.004, 0.5, 2.4, 40,
     "assumes independence; baseline only; log conditional-MLE odds ratio"),
    ("phyloglm_methanogen", "phyloglm(present~methanogen, MPLE)", "methanogen", 1.1, 0.6,
     1.83, 0.067, -0.1, 2.3, 40, "phylo logistic (MPLE); bootstrap CI (n_boot=200); log-odds"),
    ("phyloglm_confounders", "phyloglm(present~methanogen+completeness+scale(genome_size)"
     "+n_defense_excl_prometheus, MPLE)", "methanogen", 0.9, 0.4, 2.25, 0.0244, "NA", "NA",
     38, "HEADLINE: methanogen effect adjusted (log-odds)"),
    ("fitPagel", "fitPagel(methanogen, present)", "correlated_evolution", "NA", "NA", 8.2,
     0.016, "NA", "NA", 40, "Pagel94 dependent vs independent"),
    ("D_statistic", "caper::phylo.d(present)", "D", 0.2, "NA", "NA", 0.3, "NA", "NA", 40,
     "Fritz-Purvis D; p_value = P(D_random <= D_obs); caper Pval1=0.3 Pval0=0.2"),
]
COOC_H = "system_a\tsystem_b\tn_a\tn_b\tn_both\tmethod\tstatistic\tp_value\tp_adj\tdirection"
COOC_ROWS = [("Prometheus", "CBASS", 10, 12, 6, "fitPagel", 5.1, 0.04, 0.08,
              "co-occurrence (OR=3.5)"),
             ("Prometheus", "Gabija", 10, 8, 1, "fisher(fallback)", 0.25, 0.2, 0.2,
              "exclusion (OR=0.25)")]
PRED_H = "predictor\tlevel\tmodel\testimate\tstd_error\tp_value\tp_adj\tn\tnote"
ECO_ROWS = [("niche_soil", "1", "phyloglm", 0.8, 0.4, 0.045, 0.045, 40, "binary"),
            ("n_proviruses", "NA", "phyloglm", 0.1, 0.05, 0.045, 0.045, 40, "viral")]
SEL_H = "method\tstatistic\tvalue\tp_value\tn_sites_or_branches\tnote"
SEL_ROWS = [("FEL", "sites_positive", 3, "", 3, "p<=0.1, beta>alpha"),
            ("MEME", "sites_episodic", 2, "", 2, "p<=0.1"),
            ("BUSTED", "episodic_diversifying_p", 0.03, 0.03, "", "gene-wide test"),
            ("aBSREL", "branches_selected", 1, "", 1, "corrected p<=0.05"),
            ("RELAX", "K", 0.45, 0.01, "", "Test=5 terminal methanogen branches")]
REC_H = "tool\tevent\tcount\tfraction\tnote"
EVENTS = ["Duplication", "Transfer", "Loss", "Speciation", "Origination"]
REC_ROWS = ([("generax", e, {"Transfer": 5, "Speciation": 20}.get(e, 2),
              0.8 if e in ("Transfer", "Speciation") else "", "parsed *_events.txt")
             for e in EVENTS] +
            [("ale", e, 0, "", "ALE not on PATH (skipped)") for e in EVENTS])
ANC_H = "metric\tvalue\tnote"
ANC_ROWS = [("n_gains", 4.2, "mean 0->1"), ("n_losses", 7.0, "mean 1->0"),
            ("rate_gain", 0.1, "Q"), ("rate_loss", 0.2, "Q"),
            ("n_independent_acquisitions", 4.2, "mean independent gains"),
            ("gains_on_methanogen_branches", 3.1, "observed"),
            ("expected_gains_methanogen", 1.2, "null mean (Mk_simulation(ARD))"),
            ("enrichment_p", 0.004975, "one-sided P(null >= observed)")]


def producer_set(d: Path, **override):
    """Write a full set of producer tables (overridable per family) and return paths."""
    spec = {
        "association": (ASSOC_H, ASSOC_ROWS), "cooccurrence": (COOC_H, COOC_ROWS),
        "ecology": (PRED_H, ECO_ROWS), "selection": (SEL_H, SEL_ROWS),
        "reconciliation": (REC_H, REC_ROWS), "ancestral": (ANC_H, ANC_ROWS),
    }
    spec.update(override)
    return {fam: write(d / f"{fam}.tsv", h, rows) for fam, (h, rows) in spec.items()}


def registry_cli(paths, out, *extra, check=True, disorder="", methanogenesis=None):
    args = []
    for fam in ("association", "cooccurrence", "ecology", "selection", "reconciliation",
                "ancestral"):
        args += [f"--{fam}", paths[fam]]
    args += ["--disorder", disorder, "--out", out, *extra]
    if methanogenesis is not None:
        args += ["--methanogenesis", methanogenesis]
    return py("results_registry.py", *args, check=check)


# methanogenesis_stats.R columns (SCHEMA ecology/methanogenesis_associations.tsv) + extras
MGEN_COLS = ["predictor", "analysis", "role", "estimate", "std_error", "ci_low", "ci_high",
             "p_value", "p_adj", "n", "n_present", "n_absent", "n_level1", "n_changes_mean",
             "n_changes_lo", "n_changes_hi", "replicated", "delta_aic", "best_model",
             "covariates", "vif", "model", "converged", "note", "tier_set"]


def mgen_row(pred, analysis, role, tier="A", **kw):
    r = {c: "NA" for c in MGEN_COLS}
    r.update(predictor=pred, analysis=analysis, role=role, tier_set=tier, model="m", note="",
             covariates="", best_model="", n="80")
    r.update({k: str(v) for k, v in kw.items()})
    return tuple(r[c] for c in MGEN_COLS)


MGEN_ROWS = [
    mgen_row("mgen_co2_reducing", "marginal", "secondary", estimate=1.2, std_error=0.4,
             ci_low=0.416, ci_high=1.984, p_value=0.003, replicated="TRUE"),
    mgen_row("mgen_co2_reducing", "habitat_adjusted", "primary", estimate=1.1, std_error=0.45,
             ci_low=0.218, ci_high=1.982, p_value=0.0145, p_adj=0.029, replicated="TRUE"),
    mgen_row("mgen_co2_reducing", "simmap", "descriptive", n_changes_mean=12.4,
             n_changes_lo=9, n_changes_hi=16, replicated="TRUE"),
    mgen_row("mgen_co2_reducing", "corhmm", "descriptive", delta_aic=6.5,
             best_model="correlated", replicated="TRUE"),
    mgen_row("mgen_co2_reducing", "sens_unresolved_as_absent", "sensitivity", estimate=1.0,
             std_error=0.4, ci_low=0.216, ci_high=1.784, p_value=0.012),
    mgen_row("mgen_co2_reducing", "sens_tiers_AB", "sensitivity", tier="A,B", estimate=1.0,
             std_error=0.4, ci_low=0.216, ci_high=1.784, p_value=0.0001),
    mgen_row("mgen_methylotrophic", "marginal", "secondary", estimate=0.2, std_error=0.5,
             ci_low=-0.78, ci_high=1.18, p_value=0.69),
    mgen_row("mgen_methylotrophic", "habitat_adjusted", "primary", estimate=0.3, std_error=0.5,
             ci_low=-0.68, ci_high=1.28, p_value=0.55, p_adj=0.55),
    # replication guard: estimate kept, p withheld -> not counted in m
    mgen_row("mgen_acetoclastic", "habitat_adjusted", "primary", estimate=2.0, std_error=1.0,
             ci_low=0.04, ci_high=3.96, replicated="FALSE",
             note="replication guard: mean 1.1 state changes < min_changes=5"),
    mgen_row("mgen_acetoclastic", "simmap", "descriptive", n_changes_mean=1.1,
             n_changes_lo=1, n_changes_hi=2, replicated="FALSE"),
    mgen_row("mgen_acetoclastic", "corhmm", "descriptive",
             note="corHMM test disabled (--corhmm FALSE)"),
]


def reg_index(path):
    df = read_tsv(path)
    return df, {(r.family, r.test, r.term, r.tier_set): r for r in df.itertuples()}


def svg_text(path):
    return open(path, encoding="utf-8").read()


# --------------------------------------------------------------------------- #
# static checks
# --------------------------------------------------------------------------- #
@S.test
def test_py_compile_owned():
    d = tmpdir("syn_pyc_")
    for f in OWNED_PY:
        run([PYTHON, "-c", "import py_compile, sys; py_compile.compile(sys.argv[1], "
             "cfile=sys.argv[2], doraise=True)", str(f), str(d / (f.name + "c"))])


@S.test
def test_r_parse_owned():
    for f in OWNED_R:
        run([RSCRIPT, "-e", f"invisible(parse(file = '{f}'))"])


@S.test
def test_rules_hygiene_and_logs():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"),
             str(REPO / "workflow" / "rules" / "synthesis.smk")])
    assert "RULES HYGIENE OK" in p.stdout, p.stdout
    text = (REPO / "workflow" / "rules" / "synthesis.smk").read_text()
    rules = re.findall(r"^rule (\w+):", text, flags=re.M)
    assert set(rules) == {"plot_validation_tree", "master_tree_figure", "results_registry",
                          "plot_summary_forest"}, rules
    assert text.count("> {log} 2>&1") == len(rules), "every rule must redirect into its log"
    assert "|| true" not in text


# --------------------------------------------------------------------------- #
# FDR
# --------------------------------------------------------------------------- #
@S.test
def test_bh_known_answers_and_na_excluded_from_m():
    # hand-computed BH: sorted 0.005,0.01,0.03,0.04 -> 0.02,0.02,0.04,0.04
    got = rr.p_adjust([0.01, 0.04, 0.03, 0.005], "BH")
    exp = [0.02, 0.04, 0.04, 0.02]
    assert all(close(g, e) for g, e in zip(got, exp)), got
    # NA is excluded from m: m=2 -> [0.02, None, 0.04]. Counting the NA (m=3) would give
    # [0.03, None, 0.06] -- the inflation the registry must avoid.
    got = rr.p_adjust([0.01, None, 0.04], "BH")
    assert got[1] is None and close(got[0], 0.02) and close(got[2], 0.04), got
    got = rr.p_adjust(["0.01", "NA", "", float("nan"), "0.04"], "BH")
    assert close(got[0], 0.02) and close(got[4], 0.04) and got[1:4] == [None] * 3, got
    assert rr.p_adjust([0.3], "BH") == [0.3]                     # single test unchanged
    assert rr.p_adjust([None, None], "BH") == [None, None]
    assert rr.bh_correction([0.01, 0.04]) == rr.p_adjust([0.01, 0.04], "BH")
    for bad in ("qvalue", "storey", ""):
        try:
            rr.p_adjust([0.1, 0.2], bad)
        except ValueError:
            continue
        raise AssertionError(f"unsupported method {bad!r} accepted")


@S.test
def test_p_adjust_matches_R_p_adjust():
    """Real-tool reference: R's stats::p.adjust (NAs, ties, 0 and 1)."""
    rnd = random.Random(42)
    vecs = [[0.01, 0.04, 0.03, 0.005], [0.2, 0.2, 0.2, 0.01, None, 0.5],
            [0.0, 1.0, 0.5, None], [0.049, 0.05, 0.051, 0.049, 0.9]]
    for _ in range(6):
        v = [rnd.random() ** 3 for _ in range(rnd.randint(2, 40))]
        for i in rnd.sample(range(len(v)), k=len(v) // 5):
            v[i] = None
        vecs.append(v)
    d = tmpdir("syn_padj_")
    src = d / "p.tsv"
    with open(src, "w") as fh:
        for k, v in enumerate(vecs):
            for x in v:
                fh.write(f"{k}\t{'NA' if x is None else repr(x)}\n")
    out = d / "r.tsv"
    rcode = (f"x <- read.delim('{src}', header=FALSE, colClasses=c('integer','numeric'));"
             "res <- do.call(rbind, lapply(split(x, x$V1), function(g) data.frame(k=g$V1,"
             " BH=p.adjust(g$V2,'BH'), BY=p.adjust(g$V2,'BY'), holm=p.adjust(g$V2,'holm'),"
             " bonferroni=p.adjust(g$V2,'bonferroni'))));"
             f"write.table(format(res, digits=17), '{out}', sep='\\t', quote=FALSE, "
             "row.names=FALSE)")
    run([RSCRIPT, "-e", rcode])
    lines = open(out).read().splitlines()
    hdr = [h.strip() for h in lines[0].split("\t")]
    ref = {}
    for ln in lines[1:]:
        f = dict(zip(hdr, (x.strip() for x in ln.split("\t"))))
        ref.setdefault(int(f["k"]), []).append(f)
    n_cmp = 0
    for k, v in enumerate(vecs):
        for meth in ("BH", "BY", "holm", "bonferroni"):
            mine = rr.p_adjust(v, meth)
            for a, row in zip(mine, ref[k]):
                r = row[meth]
                if r == "NA":
                    assert a is None, (k, meth, a)
                else:
                    assert a is not None and close(a, float(r), 1e-12), (k, meth, a, r)
                    n_cmp += 1
    assert n_cmp > 200, n_cmp


# --------------------------------------------------------------------------- #
# registry
# --------------------------------------------------------------------------- #
@S.test
def test_registry_realistic_tagging_estimates_and_fdr():
    d = tmpdir("syn_reg_")
    paths = producer_set(d)
    out = d / "registry.tsv"
    p = registry_cli(paths, out, "--fdr", "BH", "--primary-tiers", "A",
                     "--strict-confirmatory")
    df, ix = reg_index(out)
    assert list(df.columns[:9]) == SCHEMA_COLS, list(df.columns)
    for c in rr.EXTRA_COLS:
        assert c in df.columns, c
    conf = sorted((r.family, r.test, r.term) for r in df.itertuples()
                  if r.hypothesis == "confirmatory")
    assert conf == sorted([
        ("ancestral", "methanogen_gain_enrichment",
         "methanogen_independent_acquisition_enrichment"),
        ("association", "phyloglm_confounders", "methanogen"),
        ("reconciliation", "generax", "transfer"),
        ("selection", "RELAX", "K")]), conf

    def a(t, term="methanogen"):
        return ix[("association", t, term, "")]
    # association BH family = phyloglm x2, fitPagel, D (m = 4); hand-computed:
    # sorted 0.016, 0.0244, 0.067, 0.3 -> 0.0488, 0.0488, 0.0893, 0.3
    exp = {("phyloglm_methanogen", "methanogen"): 0.067 * 4 / 3,
           ("phyloglm_confounders", "methanogen"): 0.0244 * 4 / 2,
           ("fitPagel", "correlated_evolution"): 0.0244 * 4 / 2,
           ("D_statistic", "D"): 0.3}
    for (t, term), v in exp.items():
        assert close(a(t, term).p_adj, v), (t, a(t, term).p_adj, v)
        assert a(t, term).fdr_m == "4"
    # the two independence baselines are reported but kept OUT of the BH family: their
    # tiny (invalid) p-values would otherwise LOWER the headline's adjusted p (m=6 -> 0.0366)
    for t in ("naive_glm", "naive_fisher"):
        assert a(t).role == "baseline" and a(t).hypothesis == "exploratory", t
        assert a(t).p_adj == "" and a(t).fdr_m == "", t
    wrong = rr.p_adjust([0.0027, 0.004, 0.067, 0.0244, 0.016, 0.3], "BH")[3]
    assert close(wrong, 0.0366) and not close(a("phyloglm_confounders").p_adj, wrong)
    assert a("phyloglm_methanogen").role == "secondary"
    dst = a("D_statistic", "D")
    assert dst.role == "secondary" and dst.null_value == "1.0"
    assert "D_random" in dst.note and dst.estimate_scale == "D_statistic"
    # headline: Wald CI from SE when the producer gives no CI; producer CI otherwise
    h = a("phyloglm_confounders")
    assert h.ci_method == "wald_z" and h.estimate_scale == "log_odds" and h.null_value == "0.0"
    assert close(h.ci_low, 0.9 - 1.959963984540054 * 0.4)
    assert close(h.ci_high, 0.9 + 1.959963984540054 * 0.4)
    assert a("naive_glm").ci_method == "producer" and close(a("naive_glm").ci_low, 0.6)
    assert h.n == "38" and "HEADLINE" in h.note
    # selection: BUSTED's value is its p (no estimate); RELAX K on a ratio scale
    b = ix[("selection", "BUSTED", "episodic_diversifying_p", "")]
    assert b.estimate == "" and close(b.p_adj, 0.03)
    k = ix[("selection", "RELAX", "K", "")]
    assert close(k.estimate, 0.45) and k.estimate_scale == "ratio" and k.null_value == "1.0"
    assert close(k.p_adj, 0.02) and k.fdr_m == "2"
    assert ix[("selection", "FEL", "sites_positive", "")].role == "descriptive"
    # co-occurrence: estimate = odds ratio (from 'direction'), statistic = Pagel LR
    c = ix[("cooccurrence", "fitPagel", "Prometheus_vs_CBASS", "")]
    assert close(c.estimate, 3.5) and close(c.statistic, 5.1) and close(c.p_adj, 0.08)
    assert c.estimate_scale == "odds_ratio" and "n_both=6" in c.note
    # ecology: tied p-values; Wald CI
    e = ix[("ecology", "phyloglm", "niche_soil_1", "")]
    assert close(e.p_adj, 0.045) and e.ci_method == "wald_z"
    assert ("ecology", "phyloglm", "n_proviruses", "") in ix
    # reconciliation: GeneRax count kept; ALE 'skipped' zero is NOT a zero-transfer result
    g = ix[("reconciliation", "generax", "transfer", "")]
    assert close(g.estimate, 5) and g.p_value == "" and "vertical_fraction" in g.note
    al = ix[("reconciliation", "ale", "transfer", "")]
    assert al.estimate == "" and al.note.startswith("not estimated") and \
        al.hypothesis == "exploratory"
    # ancestral: enrichment p alone in its family (m = 1); observed vs null mean
    en = ix[("ancestral", "methanogen_gain_enrichment",
             "methanogen_independent_acquisition_enrichment", "")]
    assert close(en.estimate, 3.1) and close(en.null_value, 1.2)
    assert close(en.p_adj, 0.004975) and en.fdr_m == "1" and "fold_enrichment" in en.note
    ia = ix[("ancestral", "independent_acquisitions", "independent_acquisition_count", "")]
    assert ia.hypothesis == "exploratory" and close(ia.estimate, 4.2)
    assert "4 confirmatory" in p.stderr, p.stderr
    assert len(df) == 6 + 2 + 2 + 5 + 2 + 2, len(df)
    # original five-row producer format (D without a p-value) is still accepted
    old = [r for r in ASSOC_ROWS if r[0] != "naive_fisher"]
    old[-1] = old[-1][:6] + ("NA",) + old[-1][7:]
    paths = producer_set(d / "old", association=(ASSOC_H, old))
    registry_cli(paths, d / "old.tsv", "--strict-confirmatory")
    _, ixo = reg_index(d / "old.tsv")
    do = ixo[("association", "D_statistic", "D", "")]
    assert do.role == "descriptive" and do.p_adj == ""
    assert ixo[("association", "phyloglm_confounders", "methanogen", "")].fdr_m == "3"


@S.test
def test_registry_tier_sets_are_sensitivity_not_double_counted():
    d = tmpdir("syn_tier_")
    h = ASSOC_H + "\ttier_set"
    rows = [r + ("A",) for r in ASSOC_ROWS]
    rows += [(r[0], r[1], r[2], r[3], r[4], r[5], 0.001 if r[6] != "NA" else "NA", r[7],
              r[8], 45, r[10], "B,A") for r in ASSOC_ROWS]       # A+B sensitivity re-run
    paths = producer_set(d, association=(h, rows))
    out = d / "reg.tsv"
    registry_cli(paths, out, "--primary-tiers", "A", "--strict-confirmatory")
    df, ix = reg_index(out)
    prim = ix[("association", "phyloglm_confounders", "methanogen", "A")]
    sens = ix[("association", "phyloglm_confounders", "methanogen", "A,B")]
    assert prim.hypothesis == "confirmatory" and prim.fdr_m == "4"
    assert close(prim.p_adj, 0.0244 * 4 / 2)                    # m unchanged by A+B rows
    assert sens.role == "sensitivity" and sens.hypothesis == "exploratory"
    assert sens.p_adj == "" and sens.fdr_m == "" and close(sens.p_value, 0.001)
    assert sum(1 for r in df.itertuples() if r.hypothesis == "confirmatory") == 4
    # 'A+B' and 'B,A' canonicalise to the same set: the A+B rows become primary
    out2 = d / "reg2.tsv"
    registry_cli(paths, out2, "--primary-tiers", "A+B")
    _, ix2 = reg_index(out2)
    assert ix2[("association", "phyloglm_confounders", "methanogen", "A,B")].hypothesis \
        == "confirmatory"
    assert ix2[("association", "phyloglm_confounders", "methanogen", "A")].role == \
        "sensitivity"
    # several tier sets and no primary given: ambiguous -> error
    p = registry_cli(paths, d / "reg3.tsv", check=False)
    assert p.returncode != 0 and "several tier sets" in p.stderr, p.stderr[-500:]


@S.test
def test_registry_duplicates():
    d = tmpdir("syn_dup_")
    paths = producer_set(d, cooccurrence=(COOC_H, COOC_ROWS + [COOC_ROWS[0]]))
    out = d / "reg.tsv"
    p = registry_cli(paths, out)
    df, _ = reg_index(out)
    assert sum(1 for r in df.itertuples() if r.family == "cooccurrence") == 2
    assert "identical duplicate row dropped" in p.stderr
    bad = list(COOC_ROWS[0])
    bad[7] = 0.5
    paths = producer_set(d, cooccurrence=(COOC_H, COOC_ROWS + [tuple(bad)]))
    p = registry_cli(paths, d / "reg2.tsv", check=False)
    assert p.returncode != 0 and "conflicting duplicate" in p.stderr, p.stderr[-400:]


@S.test
def test_registry_strict_confirmatory_label_mismatch():
    d = tmpdir("syn_strict_")
    renamed = [tuple("phyloglm_adjusted" if x == "phyloglm_confounders" else x for x in r)
               for r in ASSOC_ROWS]
    paths = producer_set(d, association=(ASSOC_H, renamed))
    p = registry_cli(paths, d / "r.tsv", "--strict-confirmatory", check=False)
    assert p.returncode != 0, "a renamed headline test must fail in strict mode"
    assert "confirmatory label mismatch" in p.stderr
    assert "association:phyloglm_confounders:methanogen" in p.stderr
    assert "phyloglm_adjusted:methanogen" in p.stderr          # labels found are listed
    p = registry_cli(paths, d / "r2.tsv")                       # lenient: warn only
    df, _ = reg_index(d / "r2.tsv")
    assert "WARNING: confirmatory key association:phyloglm_confounders" in p.stderr
    assert not any(r.family == "association" and r.hypothesis == "confirmatory"
                   for r in df.itertuples())
    # RELAX present but its statistic renamed -> required -> strict error
    sel = [r if r[0] != "RELAX" else ("RELAX", "k_intensity", 0.45, 0.01, "", "x")
           for r in SEL_ROWS]
    paths = producer_set(d, selection=(SEL_H, sel))
    p = registry_cli(paths, d / "r3.tsv", "--strict-confirmatory", check=False)
    assert p.returncode != 0 and "selection:relax:k" in p.stderr.lower()
    # ancestral table present but its metrics renamed: no registry row at all -> error
    anc = [("gain_enrichment_pvalue" if m == "enrichment_p" else
            "n_acquisitions" if m == "n_independent_acquisitions" else m, v, n)
           for m, v, n in ANC_ROWS]
    paths = producer_set(d, ancestral=(ANC_H, anc))
    p = registry_cli(paths, d / "r5.tsv", "--strict-confirmatory", check=False)
    assert p.returncode != 0 and "ancestral:methanogen_gain_enrichment" in p.stderr
    # RELAX disabled in config (no RELAX row) and GeneRax disabled: not required
    paths = producer_set(d, selection=(SEL_H, [r for r in SEL_ROWS if r[0] != "RELAX"]),
                         reconciliation=(REC_H, [("generax", "Speciation", 0, "",
                                                  "run_generax disabled")]))
    registry_cli(paths, d / "r4.tsv", "--strict-confirmatory")
    df, _ = reg_index(d / "r4.tsv")
    assert sum(1 for r in df.itertuples() if r.hypothesis == "confirmatory") == 2
    assert not any(r.family == "reconciliation" for r in df.itertuples())


@S.test
def test_registry_confirmatory_override_and_disorder_always_exploratory():
    d = tmpdir("syn_over_")
    paths = producer_set(d)
    dis = write(d / "disorder.tsv", PRED_H, [
        ("methanogen", "1", "phylolm(disorder~pred, lambda)", 0.03, 0.01, 0.01, 0.02, 30,
         "EXPLORATORY"),
        ("niche_soil", "1", "phylolm(disorder~pred, lambda)", -0.01, 0.01, 0.02, 0.02, 30,
         "EXPLORATORY")])
    out = d / "reg.tsv"
    p = registry_cli(paths, out, "--confirmatory",
                     "association:fitPagel:correlated_evolution;"
                     "disorder:phylolm(disorder~pred, lambda):methanogen_1",
                     disorder=dis)
    df, ix = reg_index(out)
    conf = [(r.family, r.test) for r in df.itertuples() if r.hypothesis == "confirmatory"]
    assert conf == [("association", "fitPagel")], conf
    dm = ix[("disorder", "phylolm(disorder~pred, lambda)", "methanogen_1", "")]
    assert dm.hypothesis == "exploratory" and close(dm.p_adj, 0.02) and dm.fdr_m == "2"
    assert dm.ci_method == "none" and dm.estimate_scale == "coefficient"
    assert "always exploratory" in p.stderr
    try:
        rr.parse_key("association:phyloglm")
    except ValueError:
        pass
    else:
        raise AssertionError("malformed key accepted")


@S.test
def test_registry_methanogenesis_family():
    d = tmpdir("syn_mgen_")
    paths = producer_set(d)
    mg = write(d / "methanogenesis.tsv", "\t".join(MGEN_COLS), MGEN_ROWS)
    out = d / "reg.tsv"
    p = registry_cli(paths, out, "--primary-tiers", "A", "--strict-confirmatory",
                     methanogenesis=mg)
    df, ix = reg_index(out)
    fam = df[df.family == "methanogenesis"]
    assert len(fam) == len(MGEN_ROWS) and set(fam.hypothesis) == {"exploratory"}

    def m(test, term, tier="A"):
        return ix[("methanogenesis", test, term, tier)]
    # roles kept from the producer; FDR family = the primary rows WITH a p-value (m = 2)
    assert m("habitat_adjusted", "mgen_co2_reducing").role == "primary"
    assert m("marginal", "mgen_co2_reducing").role == "secondary"
    assert m("simmap", "mgen_co2_reducing").role == "descriptive"
    assert m("sens_tiers_AB", "mgen_co2_reducing", "A,B").role == "sensitivity"
    bh = rr.p_adjust([0.0145, 0.55], "BH")
    assert close(m("habitat_adjusted", "mgen_co2_reducing").p_adj, bh[0])      # 0.029
    assert close(m("habitat_adjusted", "mgen_methylotrophic").p_adj, bh[1])
    for t, term in (("habitat_adjusted", "mgen_co2_reducing"), ("habitat_adjusted", "mgen_methylotrophic")):
        assert m(t, term).fdr_m == "2" and m(t, term).fdr_method == "BH", (t, term)
    guard = m("habitat_adjusted", "mgen_acetoclastic")
    assert guard.p_value == "" and guard.p_adj == "" and guard.fdr_m == "" and close(guard.estimate, 2.0)
    for t, term, tier in (("marginal", "mgen_co2_reducing", "A"), ("marginal", "mgen_methylotrophic", "A"),
                          ("sens_unresolved_as_absent", "mgen_co2_reducing", "A"),
                          ("sens_tiers_AB", "mgen_co2_reducing", "A,B")):
        r = m(t, term, tier)
        assert r.p_adj == "" and r.fdr_m == "" and r.p_value != "", (t, term)
    # a secondary p of 0.003 and a sensitivity p of 1e-4 would change the primary BH if counted
    wrong = rr.p_adjust([0.003, 0.0145, 0.012, 0.0001, 0.69, 0.55], "BH")[1]
    assert not close(m("habitat_adjusted", "mgen_co2_reducing").p_adj, wrong)
    # scales: log-odds (producer Wald CI), count (simmap), statistic (corHMM delta AIC)
    h = m("habitat_adjusted", "mgen_co2_reducing")
    assert h.estimate_scale == "log_odds" and h.null_value == "0.0" and h.ci_method == "producer"
    assert close(h.ci_low, 0.218) and close(h.std_error, 0.45)
    s = m("simmap", "mgen_co2_reducing")
    assert s.estimate_scale == "count" and close(s.estimate, 12.4) and close(s.ci_low, 9) and \
        close(s.ci_high, 16) and s.p_value == "" and "replicated=TRUE" in s.note
    assert "replicated=FALSE" in m("simmap", "mgen_acetoclastic").note
    c = m("corhmm", "mgen_co2_reducing")
    assert c.estimate_scale == "statistic" and close(c.estimate, 6.5) and close(c.statistic, 6.5)
    assert "best_model=correlated" in c.note and c.p_value == ""
    assert m("corhmm", "mgen_acetoclastic").estimate == ""
    # the other families are untouched (association m = 4 as without the new family)
    assert ix[("association", "phyloglm_confounders", "methanogen", "")].fdr_m == "4"
    assert "methanogenesis=2" in p.stderr, p.stderr[-600:]
    # never confirmatory, even when named
    p = registry_cli(paths, d / "reg_conf.tsv", "--primary-tiers", "A", "--confirmatory",
                     "association:phyloglm_confounders:methanogen;"
                     "methanogenesis:habitat_adjusted:mgen_co2_reducing", methanogenesis=mg)
    df2, ix2 = reg_index(d / "reg_conf.tsv")
    assert ix2[("methanogenesis", "habitat_adjusted", "mgen_co2_reducing", "A")].hypothesis == "exploratory"
    assert "always exploratory" in p.stderr
    # optional input: without --methanogenesis the family is absent; header-only = empty family
    registry_cli(paths, d / "reg_none.tsv", "--primary-tiers", "A")
    assert not (read_tsv(d / "reg_none.tsv").family == "methanogenesis").any()
    empty = write(d / "mg_empty.tsv", "\t".join(MGEN_COLS))
    registry_cli(paths, d / "reg_empty.tsv", "--primary-tiers", "A", methanogenesis=empty)
    assert not (read_tsv(d / "reg_empty.tsv").family == "methanogenesis").any()
    # roles are explicit, so several tier sets need no --primary-tiers
    registry_cli(paths, d / "reg_noprim.tsv", methanogenesis=mg)
    _, ix3 = reg_index(d / "reg_noprim.tsv")
    assert ix3[("methanogenesis", "habitat_adjusted", "mgen_co2_reducing", "A")].fdr_m == "2"
    # a registry primary tier set different from the producer's demotes every row
    p = registry_cli(paths, d / "reg_ab.tsv", "--primary-tiers", "A,B", methanogenesis=mg)
    df4, _ = reg_index(d / "reg_ab.tsv")
    f4 = df4[df4.family == "methanogenesis"]
    assert set(f4.role) == {"sensitivity"} and (f4.p_adj == "").all()
    assert "demoted to sensitivity" in p.stderr


@S.test
def test_registry_methanogenesis_malformed_inputs_fail():
    d = tmpdir("syn_mgen_bad_")
    paths = producer_set(d)
    bad_role = [mgen_row("mgen_co2_reducing", "habitat_adjusted", "secondary", estimate=1, p_value=0.1)]
    p_on_desc = [mgen_row("mgen_co2_reducing", "simmap", "descriptive", n_changes_mean=3, p_value=0.2)]
    unknown = [mgen_row("mgen_co2_reducing", "pagel", "secondary", estimate=1, p_value=0.1)]
    big_p = [mgen_row("mgen_co2_reducing", "habitat_adjusted", "primary", estimate=1, p_value=1.7)]
    bad_label = [mgen_row("mgen_co2_reducing", "marginal", "headline", estimate=1, p_value=0.1)]
    for name, rows, needle in (("role", bad_role, "must have role 'primary'"),
                               ("pdesc", p_on_desc, "carries no p-value"),
                               ("unknown", unknown, "unknown analysis 'pagel'"),
                               ("bigp", big_p, "outside [0, 1]"),
                               ("label", bad_label, "role 'headline' not in")):
        mg = write(d / f"{name}.tsv", "\t".join(MGEN_COLS), rows)
        out = d / f"{name}.out"
        p = registry_cli(paths, out, "--primary-tiers", "A", methanogenesis=mg, check=False)
        assert p.returncode != 0 and not out.exists(), name
        assert needle in p.stderr, (name, p.stderr[-600:])
    nocol = write(d / "nocol.tsv", "predictor\trole\testimate\tp_value",
                  [("mgen_co2_reducing", "primary", 1, 0.1)])
    p = registry_cli(paths, d / "nocol.out", methanogenesis=nocol, check=False)
    assert p.returncode != 0 and "lacks SCHEMA column(s) ['analysis']" in p.stderr, p.stderr[-400:]
    p = registry_cli(paths, d / "missing.out", methanogenesis=str(d / "nope.tsv"), check=False)
    assert p.returncode != 0 and not (d / "missing.out").exists()


@S.test
def test_registry_empty_inputs_give_header_only():
    d = tmpdir("syn_empty_")
    paths = producer_set(d, association=(ASSOC_H, []), cooccurrence=(COOC_H, []),
                         ecology=(PRED_H, []), selection=(SEL_H, []),
                         reconciliation=(REC_H, [("generax", "Speciation", 0, "",
                                                  "trivial species/gene tree")]),
                         ancestral=(ANC_H, []))
    out = d / "reg.tsv"
    p = registry_cli(paths, out, "--strict-confirmatory")
    lines = open(out).read().splitlines()
    assert lines == [rr.HEADER], lines
    assert "wrote 0 row(s)" in p.stderr


@S.test
def test_registry_failures_are_loud():
    d = tmpdir("syn_fail_")
    base = producer_set(d)
    cases = []
    p1 = dict(base)
    p1["selection"] = d / "missing.tsv"
    cases.append(("missing", p1, ()))
    z = d / "zero.tsv"
    z.write_text("")
    p2 = dict(base)
    p2["ancestral"] = z
    cases.append(("0-byte", p2, ()))
    nocol = write(d / "nocol.tsv", "test\tterm\testimate", [("naive_glm", "methanogen", 1)])
    p3 = dict(base)
    p3["association"] = nocol
    cases.append(("schema", p3, ()))
    badp = write(d / "badp.tsv", SEL_H, [("RELAX", "K", 0.4, 1.5, "", "x")])
    p4 = dict(base)
    p4["selection"] = badp
    cases.append(("p>1", p4, ()))
    cases.append(("fdr", base, ("--fdr", "qvalue")))
    cases.append(("disorder-missing", base, ("--disorder", str(d / "nope.tsv"))))
    for name, paths, extra in cases:
        out = d / f"{name}.out"
        p = registry_cli(paths, out, *extra, check=False)
        assert p.returncode != 0, f"{name}: expected a failure, got exit 0"
        assert not out.exists(), f"{name}: output written despite failure"


# --------------------------------------------------------------------------- #
# summary forest
# --------------------------------------------------------------------------- #
@S.test
def test_forest_one_panel_per_confirmatory_row():
    d = tmpdir("syn_forest_")
    paths = producer_set(d)
    reg = d / "registry.tsv"
    registry_cli(paths, reg, "--primary-tiers", "A")
    fig = d / "figs"
    p = py("plot_summary_forest.py", "--registry", reg, "--outdir", fig,
           "--formats", "png,svg,pdf")
    for ext in ("png", "svg", "pdf"):
        f = fig / f"summary_forest.{ext}"
        assert f.exists() and f.stat().st_size > 1000, f
    assert ("4 confirmatory panel(s): significant, significant, estimate_only, significant"
            in p.stderr), p.stderr
    s = svg_text(fig / "summary_forest.svg")
    for frag in ("Methanogen association (phyloglm,", "Horizontal transfers (GeneRax DTL)",
                 "branches (RELAX K)", "p_adj = 0.049 (BH, m=4)", "p_adj = 0.02 (BH, m=2)",
                 "estimate only", "log-odds (null 0)", "ratio, log axis (null 1)"):
        assert frag in s, frag
    for frag in ("naive", "Gabija", "niche_soil"):       # exploratory rows not plotted
        assert frag not in s, frag


@S.test
def test_forest_panel_scales_and_nulls():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import plot_summary_forest as fp

    def row(**kw):
        base = {"note": "", "p_value": "", "p_adj": "", "ci_low": "", "ci_high": "",
                "null_value": ""}
        base.update(kw)
        return pd.Series(base)
    fig, axes = plt.subplots(4, 1)
    st = fp.draw_panel(axes[0], row(estimate="0.45", estimate_scale="ratio",
                                    null_value="1.0", p_value="0.01", p_adj="0.02"), 0.05)
    lo, hi = axes[0].get_xlim()
    assert st == "significant" and axes[0].get_xscale() == "log" and lo < 0.45 < 1.0 < hi
    st = fp.draw_panel(axes[1], row(estimate="5", estimate_scale="count"), 0.05)
    assert st == "estimate_only" and axes[1].get_xlim()[0] == 0.0 and \
        axes[1].get_xscale() == "linear"
    st = fp.draw_panel(axes[2], row(estimate="0.9", ci_low="-0.2", ci_high="2.0",
                                    estimate_scale="log_odds", null_value="0.0",
                                    p_value="0.2", p_adj="0.3"), 0.05)
    lo, hi = axes[2].get_xlim()
    assert st == "tested" and lo < -0.2 and hi > 2.0
    st = fp.draw_panel(axes[3], row(estimate="", estimate_scale="log_odds",
                                    note="failed: fewer than 5 genomes"), 0.05)
    assert st == "not_estimable"
    plt.close(fig)


@S.test
def test_forest_not_estimable_placeholders_and_errors():
    d = tmpdir("syn_forest2_")
    failed = [r if r[0] != "phyloglm_confounders" else
              (r[0], r[1], r[2], "NA", "NA", "NA", "NA", "NA", "NA", 40,
               "failed: fewer than 5 genomes with complete confounders") for r in ASSOC_ROWS]
    paths = producer_set(d, association=(ASSOC_H, failed))
    reg = d / "reg.tsv"
    registry_cli(paths, reg, "--strict-confirmatory")
    p = py("plot_summary_forest.py", "--registry", reg, "--outdir", d / "f1")
    assert "not_estimable" in p.stderr and "4 confirmatory panel(s)" in p.stderr, p.stderr
    assert "not estimable: failed: fewer than 5 genomes" in svg_text(
        d / "f1" / "summary_forest.svg")
    # header-only registry -> labelled placeholder (genuinely empty input)
    empty = d / "empty.tsv"
    empty.write_text(rr.HEADER + "\n")
    py("plot_summary_forest.py", "--registry", empty, "--outdir", d / "f2")
    assert "results registry has no tests" in svg_text(d / "f2" / "summary_forest.svg")
    # only exploratory rows -> placeholder naming the reason
    df = read_tsv(reg)
    df = df[df.hypothesis == "exploratory"]
    expl = d / "expl.tsv"
    df.to_csv(expl, sep="\t", index=False)
    py("plot_summary_forest.py", "--registry", expl, "--outdir", d / "f3")
    assert "no confirmatory result" in svg_text(d / "f3" / "summary_forest.svg")
    # missing registry / missing columns / bad alpha -> non-zero exit
    for args in (["--registry", d / "nope.tsv"],
                 ["--registry", write(d / "bad.tsv", "family\ttest", [("a", "b")])],
                 ["--registry", reg, "--alpha", "1.5"]):
        q = py("plot_summary_forest.py", *args, "--outdir", d / "f4", check=False)
        assert q.returncode != 0, args


# --------------------------------------------------------------------------- #
# master tree (real ggtree / gheatmap / ggnewscale)
# --------------------------------------------------------------------------- #
MASTER_H = ("genome_id\trepresentative\ttree_tip\tprometheus_present\tprometheus_tier\t"
            "methanogen\tis_anme\tmethanogenesis_subtype\tniche_soil\tniche_marine\tgroup_Marine")
MASTER_ROWS = [
    ("G1", "True", "G1", "True", "A", "True", "False", "hydrogenotrophic", 1, 0, 0),
    ("G2", "True", "G2", "True", "B", "True", "False", "methylotrophic", 1, 1, 1),  # B-only
    ("G3", "True", "G3", "False", "", "False", "True", "NA", 0, 1, 1),             # ANME
    ("G4", "True", "G4", "True", "C", "False", "False", "NA", 0, 0, 0),
    ("G5", "True", "G5", "False", "", "False", "False", "NA", 1, "", 0),   # marine unknown
    ("G6", "True", "G6", "True", "A", "True", "False", "NA", 0, 1, 1),
    ("G7", "True", "G7", "False", "", "False", "False", "NA", "", "", ""),  # all niches unknown
    ("GENOME8", "True", "G8", "True", "A", "True", "False", "acetoclastic", 1, 1, 1),
]  # G9 has no master row
TREE9 = ("(((('G1':0.1,G2:0.2):0.1,(G3:0.15,G4:0.1):0.2):0.1,(G5:0.3,G6:0.1):0.2):0.05,"
         "(G7:0.2,(G8:0.1,G9:0.12):0.1):0.1);")
VIRAL_H = ("genome_id\tn_proviruses\ttotal_provirus_kb\tmean_checkv_completeness\t"
           "n_rloop_proviruses\tfrac_rloop_proviruses")
VIRAL_ROWS = [("G1", 2, 30, 80, 1, 0.5), ("G2", 0, 0, "", 0, 0), ("G3", 5, 100, 90, 2, 0.4)]


def tree_inputs(d: Path, tree=TREE9, master_rows=MASTER_ROWS, master_h=MASTER_H,
                viral_rows=VIRAL_ROWS):
    (d / "tree.nwk").write_text(tree + "\n")
    write(d / "master.tsv", master_h, master_rows)
    write(d / "viral.tsv", VIRAL_H, viral_rows)
    write(d / "gl.tsv", ANC_H, ANC_ROWS[:4] + [("n_independent_acquisitions", 3.2, "mean")])
    return d


def master_tree(d: Path, check=True, **kw):
    a = dict(tree=d / "tree.nwk", master=d / "master.tsv", viral=d / "viral.tsv",
             gains=d / "gl.tsv")
    a.update(kw)
    return rscript("plot_tree_master.R", "--tree", a["tree"], "--master", a["master"],
                   "--viral", a["viral"], "--gains", a["gains"],
                   "--plot-prefix", d / "fig" / "master_tree",
                   "--annotation-out", d / "fig" / "ann.tsv", check=check, timeout=600)


@S.test
def test_master_tree_rings_and_tip_matching():
    run([RSCRIPT, "-e", "for (p in c('ape','ggtree','ggnewscale','ggplot2','optparse',"
         "'svglite','ragg')) if (!requireNamespace(p, quietly=TRUE)) stop('missing R "
         "package: ', p)"])
    d = tree_inputs(tmpdir("syn_mt_"))
    p = master_tree(d)
    for ext in ("png", "svg"):
        f = d / "fig" / f"master_tree.{ext}"
        assert f.exists() and f.stat().st_size > 5000, f
    assert "9 tips (8 matched), rectangular layout, 7 rings" in p.stderr, p.stderr
    assert "1/9 tip(s) have no master row" in p.stderr
    ann = {r.label: r for r in read_tsv(d / "fig" / "ann.tsv").itertuples()}
    assert set(ann) == {f"G{i}" for i in range(1, 10)}          # quoted 'G1' matched as G1
    exp = {  # label: (tier, metabolism, subtype, niche_richness, group_Marine, n_prov, frac)
        "G1": ("A", "methanogen", "hydrogenotrophic", "1", "no", "2", "0.5"),
        "G2": ("B", "methanogen", "methylotrophic", "2", "yes", "0", ""),   # frac undefined
        "G3": ("absent", "ANME (methanotroph)", "not methanogen", "1", "yes", "5", "0.4"),
        "G4": ("C", "neither", "not methanogen", "0", "no", "", ""),
        "G5": ("absent", "neither", "not methanogen", "1", "no", "", ""),   # NA adds nothing
        "G6": ("A", "methanogen", "unresolved", "1", "yes", "", ""),
        "G7": ("absent", "neither", "not methanogen", "", "no data", "", ""),  # unknown != 0
        "G8": ("A", "methanogen", "acetoclastic", "2", "yes", "", ""),     # via tree_tip
        "G9": ("no data", "no data", "no data", "", "no data", "", ""),    # not an absence
    }
    for lab, want in exp.items():
        r = ann[lab]
        got = (r.tier, r.metabolism, r.subtype, r.niche_richness, r.group_Marine,
               r.n_proviruses, r.frac_rloop)
        assert got == want, (lab, got)
    assert ann["G9"].matched == "FALSE" and ann["G8"].matched == "TRUE"
    s = svg_text(d / "fig" / "master_tree.svg")
    assert "mean independent acquisitions = 3.20" in s and "ANME (methanotroph)" in s
    # one legend per ring (ggnewscale): without it the fill scales would collapse into one
    for title in ("Prometheus tier", "mcr metabolism", "Methanogenesis", "Niche richness",
                  "In Marine", "Proviral load", "R-loop-prone"):
        assert title in s, title


def random_newick(labels, rnd):
    nodes = [f"{lab}:{rnd.uniform(0.05, 0.5):.3f}" for lab in labels]
    while len(nodes) > 1:
        i, j = sorted(rnd.sample(range(len(nodes)), 2))
        b = nodes.pop(j)
        a = nodes.pop(i)
        nodes.append(f"({a},{b}):{rnd.uniform(0.01, 0.3):.3f}")
    return nodes[0].rsplit(":", 1)[0] + ";"


@S.test
def test_master_tree_circular_quoted_long_labels():
    rnd = random.Random(3)
    ids = [f"GCA_{rnd.randrange(10**9):09d}.1" + ("_long_strain_name" if i % 5 == 0 else "")
           for i in range(60)]
    labels = [f"'{x}'" if i % 2 else x for i, x in enumerate(ids)]
    tiers = ["A", "B", "C", ""]
    rows = [(g, "True", g, "True" if tiers[i % 4] else "False", tiers[i % 4],
             "True" if i % 3 == 0 else "False", "False",
             "hydrogenotrophic" if i % 3 == 0 else "NA", i % 2, (i + 1) % 2, i % 2)
            for i, g in enumerate(ids)]
    d = tree_inputs(tmpdir("syn_mtc_"), tree=random_newick(labels, rnd),
                    master_rows=rows[2:],
                    viral_rows=[(g, i % 4, 1, 1, 0, 0) for i, g in enumerate(ids)])
    p = master_tree(d)
    assert "60 tips (58 matched), circular layout" in p.stderr, p.stderr
    ann = read_tsv(d / "fig" / "ann.tsv")
    assert len(ann) == 60 and (ann.matched == "TRUE").sum() == 58
    assert set(ann[ann.matched == "FALSE"].label) == set(ids[:2])
    assert (ann[ann.label.isin(ids[2:])].tier == "no data").sum() == 0
    want_b = {g for i, g in enumerate(ids[2:], start=2) if tiers[i % 4] == "B"}
    assert set(ann[ann.tier == "B"].label) == want_b             # B stays B, not absent


@S.test
def test_master_tree_failures_and_empty_input():
    base = tree_inputs(tmpdir("syn_mtf_"))
    cases = {
        "no species-tree tip matches": dict(master_rows=[
            ("X" + r[0], "True", "") + r[3:] for r in MASTER_ROWS]),
        "lacks column": dict(master_h=MASTER_H.replace("prometheus_tier", "tier_x")),
        "unexpected prometheus_tier": dict(master_rows=MASTER_ROWS[:-1] + [
            ("GENOME8", "True", "G8", "True", "Rejected", "True", "False", "NA", 0, 0, 0)]),
        "malformed tree": dict(tree="((G1:0.1,G2:0.2),(G3:0.1,G4:0.1;"),
    }
    zero = base / "zero.tsv"
    zero.write_text("")
    jobs = [(msg, tree_inputs(tmpdir("syn_mtf_"), **kw), {}) for msg, kw in cases.items()]
    jobs += [("viral not found", base, {"viral": base / "nope.tsv"}),
             ("0-byte", base, {"master": zero})]
    with ThreadPoolExecutor(max_workers=2) as ex:      # two R processes at a time
        res = list(ex.map(lambda j: (j[0], master_tree(j[1], check=False, **j[2])), jobs))
    for msg, p in res:
        assert p.returncode != 0, f"{msg}: expected failure"
        if msg != "malformed tree":
            assert msg in p.stderr, (msg, p.stderr[-600:])
    # header-only viral summary: rings degrade to 'no data' with a warning (exit 0)
    hv = tree_inputs(tmpdir("syn_mtv_"), viral_rows=[])
    p = master_tree(hv)
    assert "viral summary is header-only" in p.stderr
    ann = read_tsv(hv / "fig" / "ann.tsv")
    assert (ann.n_proviruses == "").all() and (ann.frac_rloop == "").all()
    # genuinely empty input: '();' tree -> placeholder + header-only annotation, exit 0
    e = tree_inputs(tmpdir("syn_mte_"), tree="();")
    master_tree(e)
    assert (e / "fig" / "master_tree.png").exists() and (e / "fig" / "master_tree.svg").exists()
    ann = open(e / "fig" / "ann.tsv").read().splitlines()
    assert len(ann) == 1 and ann[0].startswith("label\tmatched\ttier"), ann


@S.test
def test_master_tree_pinned_rstats_stack():
    """Same figure under the versions workflow/envs/rstats.yaml pins (ggplot2 3.5.1,
    ggnewscale 0.4.10, ggtree 3.10.0, patchwork 1.2.0), when that private env exists."""
    rs = os.environ.get("PROM_PINNED_RSCRIPT") or str(
        Path(PYTHON).resolve().parents[2] / "prom-synthesis-r" / "bin" / "Rscript")
    if not os.path.exists(rs):
        print("    note: pinned rstats env not found (set PROM_PINNED_RSCRIPT); not run")
        return
    env = dict(os.environ)
    env.pop("R_LIBS", None)
    ver = subprocess.run([rs, "-e", "cat(as.character(packageVersion('ggplot2')), "
                          "as.character(packageVersion('ggnewscale')), "
                          "as.character(packageVersion('ggtree')))"],
                         capture_output=True, text=True, env=env, check=True).stdout.split()
    assert ver == ["3.5.1", "0.4.10", "3.10.0"], ver
    def render(tag_r):
        tag, r = tag_r
        d = tree_inputs(tmpdir(f"syn_pin_{tag}_"))
        e = env if tag == "pinned" else dict(os.environ)
        p = subprocess.run([r, str(SCRIPTS_R / "plot_tree_master.R"), "--tree", d / "tree.nwk",
                            "--master", d / "master.tsv", "--viral", d / "viral.tsv",
                            "--gains", d / "gl.tsv", "--plot-prefix", d / "fig" / "mt",
                            "--annotation-out", d / "fig" / "ann.tsv"],
                           capture_output=True, text=True, env=e, timeout=600)
        assert p.returncode == 0, (tag, p.stderr[-1500:])
        assert (d / "fig" / "mt.png").stat().st_size > 5000
        return open(d / "fig" / "ann.tsv").read()
    with ThreadPoolExecutor(max_workers=2) as ex:
        outs = list(ex.map(render, (("pinned", rs), ("toolchain", RSCRIPT))))
    assert outs[0] == outs[1], "ring values differ between the pinned and toolchain stacks"


@S.test
def test_save_plot_render_error_fails():
    """plotting_theme.R save_plot(): a render error FAILS the call (non-zero exit), names the
    error and leaves no partial figure. (Review: it used to draw a blank 'Figure unavailable'
    panel and return normally, so a broken figure could pass as a finished output; callers
    such as plot_tree_master.R had to intercept the message to notice.)"""
    d = tmpdir("syn_sp_")
    code = (f"source('{SCRIPTS_R / 'plotting_theme.R'}'); "
            f"save_plot(function() stop('boom'), '{d / 'x'}'); cat('NOT REACHED')")
    p = run([RSCRIPT, "-e", code], check=False)
    out = (p.stdout or "") + (p.stderr or "")
    assert p.returncode != 0 and "NOT REACHED" not in out, out
    assert "render error: boom" in out and "failed to render" in out, out
    assert not list(d.glob("x.*")), sorted(d.iterdir())


# --------------------------------------------------------------------------- #
# Snakemake dry run of the synthesis rules (private copy of the repository)
# --------------------------------------------------------------------------- #
def _snakemake_bin():
    cand = [os.environ.get("PROM_SNAKEMAKE", ""), shutil.which("snakemake") or "",
            str(Path(PYTHON).resolve().parents[2] / "prom-smk" / "bin" / "snakemake")]
    return next((c for c in cand if c and os.path.exists(c)), None)


def _copy_repo(dst: Path, pristine_rules=None):
    buf = io.BytesIO()
    skip = {".snakemake", ".git", "results", "locks", ".tmp"}
    with tarfile.open(fileobj=buf, mode="w") as tf:
        tf.add(REPO, arcname=".",
               filter=lambda ti: None if Path(ti.name).name in skip else ti)
    buf.seek(0)
    with tarfile.open(fileobj=buf, mode="r") as tf:
        tf.extractall(dst)
    if pristine_rules is not None:          # isolate synthesis.smk from concurrent edits
        for f in pristine_rules.glob("*.smk"):
            if f.name != "synthesis.smk":
                shutil.copy(f, dst / "workflow" / "rules" / f.name)
        shutil.copy(pristine_rules.parent / "Snakefile", dst / "workflow" / "Snakefile")


def _dry_run(smk, root: Path):
    env = dict(os.environ)
    home = Path(tempfile.mkdtemp(prefix="syn_home_"))
    env.update(HOME=str(home), XDG_CACHE_HOME=str(home / ".cache"))
    subprocess.run([PYTHON, ".test/setup_dryrun.py"], cwd=root, env=env, check=True,
                   capture_output=True, text=True)
    t = ".test/results"
    return subprocess.run(
        [smk, f"{t}/synthesis/results_registry.tsv", f"{t}/synthesis/figures/master_tree.png",
         f"{t}/synthesis/figures/summary_forest.svg",
         f"{t}/prometheus/plots/prometheus_validation_tree.png",
         "-n", "-p", "--nolock", "-s", "workflow/Snakefile",
         "--configfile", ".test/config/config.yaml"],
        cwd=root, env=env, capture_output=True, text=True, timeout=600)


@S.test
def test_snakemake_dry_run_synthesis_rules():
    smk = _snakemake_bin()
    if smk is None:
        print("    note: snakemake not found (set PROM_SNAKEMAKE); dry run not executed")
        return
    root = Path(tempfile.mkdtemp(prefix="syn_dag_"))
    _copy_repo(root)
    p = _dry_run(smk, root)
    out = p.stdout + p.stderr
    if p.returncode != 0:
        if "synthesis.smk" in out:
            raise AssertionError("dry run failed in synthesis.smk:\n" + out[-3000:])
        pristine = REPO.parent / (REPO.name + "_pristine_56f5437") / "workflow" / "rules"
        if not pristine.is_dir():
            raise AssertionError("dry run failed outside synthesis.smk and no pristine "
                                 "copy is available to isolate it:\n" + out[-2000:])
        print("    note: current DAG fails outside synthesis.smk (concurrent edits); "
              "re-checking synthesis.smk against the pristine stage rule files")
        root = Path(tempfile.mkdtemp(prefix="syn_dag2_"))
        _copy_repo(root, pristine_rules=pristine)
        p = _dry_run(smk, root)
        out = p.stdout + p.stderr
        assert p.returncode == 0, out[-3000:]
    for rule in ("results_registry", "plot_summary_forest", "master_tree_figure",
                 "plot_validation_tree"):
        m = re.search(rf"^rule {rule}:\n(.*?)\n\n", out, flags=re.S | re.M)
        assert m, f"rule {rule} not scheduled"
        assert "log: .test/logs/synthesis/" in m.group(1), (rule, m.group(1))
    assert "--strict-confirmatory" in out and '--primary-tiers "A"' in out
    assert '--disorder ".test/results/disorder/disorder_associations.tsv"' in out
    # the within-methanogen pathway tests feed the registry (always-exploratory family)
    assert ('--methanogenesis ".test/results/ecology/methanogenesis_associations.tsv"' in out
            or "--methanogenesis .test/results/ecology/methanogenesis_associations.tsv" in out), \
        out[-3000:]
    m = re.search(r"^rule methanogenesis_stats:\n(.*?)\n\n", out, flags=re.S | re.M)
    assert m and "log: .test/logs/ecology/methanogenesis_stats.log" in m.group(1), out[-3000:]
    assert "workflow/scripts/R/methanogenesis_stats.R" in out and "--corhmm TRUE" in out
    assert "master_tree_annotation.tsv" in out and "--formats png,svg" in out


if __name__ == "__main__":
    S.main("SYNTHESIS TESTS PASSED")
