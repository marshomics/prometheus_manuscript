"""Stage 4 (distribution) synthetic test suite.

Run from the repo root:  bash tests/run_leaf.sh distribution
Prints DISTRIBUTION TESTS PASSED only when every test passes.

20 tests. Covers build_defense_matrix.py, compute_prevalence.py, association_tests.R,
cooccurrence.R, plot_distribution.py and plot_distribution_extra.py on inputs with KNOWN
answers (hand-built tables; trees and traits simulated by sim_data.R with a known
generative model), their empty-input and failure behaviour, and real runs of the R
tools (phylolm, phytools, caper) checked against independent direct calls
(verify_direct.R). At most 2 threads per tool.
"""
from __future__ import annotations

import sys
sys.dont_write_bytecode = True   # importing pipeline modules must not write into workflow/scripts/py/__pycache__

import json
import math
import os
import sys
import functools
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _env import (Suite, py, rscript, run, tmpdir, read_tsv, REPO, SCRIPTS_PY, SCRIPTS_R,  # noqa: E402
                  PYTHON, RSCRIPT)

HERE = Path(__file__).resolve().parent
SIM = HERE / "sim_data.R"
VERIFY = HERE / "verify_direct.R"
suite = Suite("distribution")

OWNED_PY = ["build_defense_matrix.py", "compute_prevalence.py", "plot_distribution.py",
            "plot_distribution_extra.py"]
OWNED_R = ["association_tests.R", "cooccurrence.R"]
ASSOC_COLS = ["test", "model", "term", "estimate", "std_error", "statistic", "p_value",
              "ci_low", "ci_high", "n", "note"]
ASSOC_ORDER = ["naive_glm", "naive_fisher", "phyloglm_methanogen", "phyloglm_confounders",
               "fitPagel", "D_statistic"]
COOC_BASE = ["system_a", "system_b", "n_a", "n_b", "n_both", "method", "statistic", "p_value",
             "p_adj", "direction"]
PREV_BASE = ["rank", "taxon", "n_genomes", "n_present", "prevalence", "ci_low", "ci_high",
             "median_completeness", "tier_set"]


def fails(res, *needles):
    """Assert a CompletedProcess failed (non-zero exit) and mentions every needle."""
    assert res.returncode != 0, f"expected failure, got exit 0\n{res.stdout[-800:]}\n{res.stderr[-800:]}"
    blob = (res.stdout + res.stderr).lower()
    for n in needles:
        assert n.lower() in blob, f"failure message lacks {n!r}:\n{res.stderr[-1500:]}"


def close(a, b, tol, what):
    assert abs(float(a) - float(b)) <= tol, f"{what}: {a} vs {b} (tol {tol})"


def fnum(x):
    return float("nan") if x in ("", "NA") else float(x)


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(v) for v in r) + "\n")


@functools.lru_cache(maxsize=None)
def sim(scenario: str, n: int, seed: int) -> Path:
    d = tmpdir(f"dist_sim_{scenario}_")
    run([RSCRIPT, str(SIM), str(d), scenario, str(n), str(seed)])
    return d


def expected_sets(master: Path, tree: Path, tiers):
    """Independent (Python) implementation of the analysis set + Contract v2 item 6 rule."""
    import re
    txt = open(tree).read()
    tips = set(re.findall(r"[(,]\s*'?([^():,;'\[\]]+)'?\s*:", txt))
    m = read_tsv(master)
    cand = m[(m.pass_primary == "True") & (m.tree_tip != "") & m.tree_tip.isin(tips)]
    excl = cand[(cand.prometheus_tier != "") & ~cand.prometheus_tier.isin(tiers)]
    use = cand[(cand.prometheus_tier == "") | cand.prometheus_tier.isin(tiers)]
    return use, excl


# --------------------------------------------------------------------------- #
# static checks
# --------------------------------------------------------------------------- #
@suite.test
def test_syntax_python_and_r():
    for f in OWNED_PY:   # syntax check in memory (no .pyc written outside the owned paths)
        compile(open(SCRIPTS_PY / f).read(), str(SCRIPTS_PY / f), "exec")
    for f in OWNED_R + [str(SIM), str(VERIFY)]:
        p = f if os.path.isabs(f) else str(SCRIPTS_R / f)
        run([RSCRIPT, "-e", f"invisible(parse('{p}'))"])


@suite.test
def test_rules_hygiene_and_logs():
    r = run([PYTHON, str(REPO / "tests" / "check_rules.py"), str(REPO / "workflow/rules/distribution.smk")])
    assert "RULES HYGIENE OK" in r.stdout
    text = open(REPO / "workflow/rules/distribution.smk").read()
    assert "|| true" not in text and "2>/dev/null" not in text
    assert "--require-all" in text, "archaeal defense_matrix must require every DefenseFinder table"
    assert text.count("> {log} 2>&1") == 5, "every rule must send stdout+stderr to its log"


# --------------------------------------------------------------------------- #
# compute_prevalence.py
# --------------------------------------------------------------------------- #
def prevalence_master(path):
    """P1: 3xA, 1xB, 1xC, 5 absent (pass) + 2 tier-A genomes failing pass_primary.
    P2: 4 absent. Two P1 genomes have no genus. Completeness known per genome."""
    hdr = ["genome_id", "tree_tip", "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
           "gtdb_family", "gtdb_genus", "gtdb_species", "completeness", "contamination",
           "pass_primary", "prometheus_present", "prometheus_tier", "prometheus_n_loci"]
    rows = []
    spec = ([("A", "True")] * 3 + [("B", "True"), ("C", "True")] + [("", "True")] * 5
            + [("A", "False")] * 2)
    for i, (tier, pp) in enumerate(spec):
        genus = "" if i in (5, 6) else "G1"
        rows.append([f"0{i:03d}", f"0{i:03d}", "Archaea", "P1", "c", "o", "F1", genus, "s",
                     90 + i, 1, pp, "True" if tier else "False", tier, int(bool(tier))])
    for j in range(4):
        rows.append([f"1{j:03d}", f"1{j:03d}", "Archaea", "P2", "c", "o", "F2", "G2", "s",
                     80, 1, "True", "False", "", 0])
    write_tsv(path, hdr, rows)
    return hdr


@suite.test
def test_wilson_matches_scipy():
    sys.path.insert(0, str(SCRIPTS_PY))
    import compute_prevalence as cp
    from scipy.stats import binomtest
    for k, n in [(0, 7), (3, 8), (4, 9), (8, 8), (1, 1), (17, 230)]:
        lo, hi = cp.wilson_ci(n, k)
        ref = binomtest(k, n).proportion_ci(confidence_level=0.95, method="wilson")
        close(lo, ref.low, 1e-12, f"wilson low {k}/{n}")
        close(hi, ref.high, 1e-12, f"wilson high {k}/{n}")
    assert all(math.isnan(v) for v in cp.wilson_ci(0, 0))


@suite.test
def test_prevalence_tier_rule_known_answer():
    from scipy.stats import binomtest
    d = tmpdir("dist_prev_")
    prevalence_master(d / "m.tsv")
    py("compute_prevalence.py", "--master", d / "m.tsv", "--primary-tiers", "A", "--out", d / "p.tsv")
    p = read_tsv(d / "p.tsv")
    assert list(p.columns[:9]) == PREV_BASE, list(p.columns)
    p1 = p[(p["rank"] == "gtdb_phylum") & (p.taxon == "P1")].iloc[0]
    # tier A: B-only and C-only genomes are EXCLUDED (not absences); non-pass genomes ignored
    assert (int(p1.n_genomes), int(p1.n_present), int(p1.n_excluded), int(p1.n_genomes_total)) == (8, 3, 2, 10)
    close(p1.prevalence, 3 / 8, 1e-6, "P1 prevalence tier A")
    ref = binomtest(3, 8).proportion_ci(method="wilson")
    close(p1.ci_low, ref.low, 1e-6, "P1 ci_low"); close(p1.ci_high, ref.high, 1e-6, "P1 ci_high")
    # median completeness over the 8 denominator genomes (90,91,92,95..99)
    close(p1.median_completeness, 95.5, 1e-9, "P1 median completeness")
    # sensitivity A+B: the B-only genome is PRESENT, the C-only one still excluded
    assert (int(p1.n_genomes_sensitivity), int(p1.n_present_sensitivity),
            int(p1.n_excluded_sensitivity)) == (9, 4, 1)
    assert p1.tier_set == "A" and p1.tier_set_sensitivity == "A,B"
    ov = p[p["rank"] == "overall"].iloc[0]
    assert (int(ov.n_genomes), int(ov.n_present), int(ov.n_excluded)) == (12, 3, 2)
    p2 = p[(p["rank"] == "gtdb_phylum") & (p.taxon == "P2")].iloc[0]
    assert int(p2.n_present) == 0 and fnum(p2.ci_low) == 0.0 and fnum(p2.prevalence) == 0.0
    g1 = p[(p["rank"] == "gtdb_genus") & (p.taxon == "G1")].iloc[0]
    assert int(g1.n_genomes_total) == 8, "genomes without a genus must be left out of genus rows"
    assert set(p["rank"]) == {"overall", "gtdb_phylum", "gtdb_class", "gtdb_order", "gtdb_family", "gtdb_genus"}


@suite.test
def test_prevalence_primary_AB_and_bacterial_cli():
    d = tmpdir("dist_prev_ab_")
    prevalence_master(d / "m.tsv")
    py("compute_prevalence.py", "--master", d / "m.tsv", "--primary-tiers", "A,B", "--out", d / "p.tsv")
    p1 = read_tsv(d / "p.tsv").query("taxon == 'P1'").iloc[0]
    assert (int(p1.n_genomes), int(p1.n_present), int(p1.n_excluded)) == (9, 4, 1)
    # bacterial slim master (build_bacterial_table.py columns), bact_prevalence CLI
    hdr = ["genome_id", "tree_tip", "gtdb_domain", "gtdb_phylum", "gtdb_class", "gtdb_order",
           "gtdb_family", "gtdb_genus", "gtdb_species", "completeness", "contamination", "genome_size",
           "pass_primary", "prometheus_present", "prometheus_tier", "prometheus_n_loci",
           "n_defense_systems", "n_defense_systems_excl_prometheus"]
    rows = [[f"GCF_{i}", f"GCF_{i}", "Bacteria", "Pseudomonadota", "c", "o", "f", "g", "s", 95, 1, 4e6,
             "True" if i != 2 else "False", "True" if i < 3 else "False", "A" if i < 3 else "", int(i < 3), 3, 2]
            for i in range(6)]
    write_tsv(d / "b.tsv", hdr, rows)
    py("compute_prevalence.py", "--master", d / "b.tsv", "--primary-tiers", "A", "--out", d / "bp.tsv")
    ov = read_tsv(d / "bp.tsv").query("rank == 'overall'").iloc[0]
    assert (int(ov.n_genomes), int(ov.n_present)) == (5, 2)


@suite.test
def test_prevalence_empty_and_failures():
    d = tmpdir("dist_prev_bad_")
    hdr = prevalence_master(d / "m.tsv")
    write_tsv(d / "empty.tsv", hdr, [])
    py("compute_prevalence.py", "--master", d / "empty.tsv", "--primary-tiers", "A", "--out", d / "e.tsv")
    first = open(d / "e.tsv").read().splitlines()
    assert len(first) == 1 and first[0].split("\t")[:9] == PREV_BASE
    m = read_tsv(d / "m.tsv")
    m.drop(columns=["prometheus_tier"]).to_csv(d / "nocol.tsv", sep="\t", index=False)
    fails(py("compute_prevalence.py", "--master", d / "nocol.tsv", "--primary-tiers", "A",
             "--out", d / "x.tsv", check=False), "prometheus_tier")
    bad = m.copy(); bad.loc[0, "prometheus_tier"] = ""          # present=True but no tier
    bad.to_csv(d / "incons.tsv", sep="\t", index=False)
    fails(py("compute_prevalence.py", "--master", d / "incons.tsv", "--primary-tiers", "A",
             "--out", d / "x.tsv", check=False), "disagrees")
    bad = m.copy(); bad.loc[0, "prometheus_tier"] = "D"
    bad.to_csv(d / "tierD.tsv", sep="\t", index=False)
    fails(py("compute_prevalence.py", "--master", d / "tierD.tsv", "--primary-tiers", "A",
             "--out", d / "x.tsv", check=False), "unexpected")
    fails(py("compute_prevalence.py", "--master", d / "m.tsv", "--primary-tiers", "X",
             "--out", d / "x.tsv", check=False), "subset of a,b,c")
    bad = m.copy(); bad.loc[0, "pass_primary"] = "maybe"
    bad.to_csv(d / "pp.tsv", sep="\t", index=False)
    fails(py("compute_prevalence.py", "--master", d / "pp.tsv", "--primary-tiers", "A",
             "--out", d / "x.tsv", check=False), "non-boolean")
    fails(py("compute_prevalence.py", "--master", d / "missing.tsv", "--primary-tiers", "A",
             "--out", d / "x.tsv", check=False))


# --------------------------------------------------------------------------- #
# build_defense_matrix.py
# --------------------------------------------------------------------------- #
DF2_HDR = ["sys_id", "type", "subtype", "activity", "sys_beg", "sys_end", "protein_in_syst",
           "genes_count", "name_of_profiles_in_sys"]
PLACEHOLDER_HDR = ["sys_id", "type", "subtype", "sys_beg", "sys_end", "protein_in_syst",
                   "name_of_profiles_in_sys"]


def df_dir_fixture(root: Path):
    def put(g, hdr, rows):
        (root / g).mkdir(parents=True, exist_ok=True)
        write_tsv(root / g / f"{g}_defense_finder_systems.tsv", hdr, rows)
    put("000123", DF2_HDR, [
        ["000123_1_1_RM", "RM", "RM_Type_I", "Defense", "000123_1_1", "000123_1_3", "a,b,c", 3, "x"],
        ["000123_1_4_Gab", "Gabija", "Gabija", "Defense", "000123_1_4", "000123_1_5", "d,e", 2, "y"],
        ["000123_2_7_Pro", "Prometheus", "Prometheus", "Defense", "000123_2_7", "000123_2_7", "f", 1, "Prometheus__ProA"],
        ["000123_3_1_ADF", "Anti_CBASS", "Acb1", "Antidefense", "000123_3_1", "000123_3_1", "g", 1, "z"]])
    put("g2", DF2_HDR, [])                                   # DefenseFinder found nothing
    put("g3", PLACEHOLDER_HDR, [["g3_1_2_CB", "CBASS", "CBASS_I", "g3_1_2", "g3_1_4", "h,i", "w"]])
    write_tsv(root.parent / "master.tsv", ["genome_id", "pass_primary"],
              [["000123", "True"], ["g2", "True"], ["g3", "False"]])


@suite.test
def test_defense_matrix_known_answer():
    d = tmpdir("dist_dm_")
    df_dir_fixture(d / "df")
    py("build_defense_matrix.py", "--master", d / "master.tsv", "--df-dir", d / "df",
       "--require-all", "--out", d / "mat.tsv")
    mat = read_tsv(d / "mat.tsv")
    assert list(mat.columns) == ["genome_id", "CBASS", "Gabija", "Prometheus", "RM"], list(mat.columns)
    assert list(mat.genome_id) == ["000123", "g2", "g3"], "ids must stay strings (leading zeros)"
    got = {r.genome_id: (r.CBASS, r.Gabija, r.Prometheus, r.RM) for r in mat.itertuples()}
    assert got == {"000123": ("0", "1", "1", "1"), "g2": ("0", "0", "0", "0"), "g3": ("1", "0", "0", "0")}


@suite.test
def test_defense_matrix_missing_empty_and_failures():
    d = tmpdir("dist_dm_bad_")
    df_dir_fixture(d / "df")
    write_tsv(d / "m4.tsv", ["genome_id"], [["000123"], ["g2"], ["g3"], ["g4"]])
    fails(py("build_defense_matrix.py", "--master", d / "m4.tsv", "--df-dir", d / "df",
             "--require-all", "--out", d / "x.tsv", check=False), "missing")
    r = py("build_defense_matrix.py", "--master", d / "m4.tsv", "--df-dir", d / "df", "--out", d / "b.tsv")
    assert "1 genome(s) without a systems table" in r.stderr
    assert read_tsv(d / "b.tsv").set_index("genome_id").loc["g4"].astype(int).sum() == 0
    (d / "df" / "g2" / "g2_defense_finder_systems.tsv").write_text("")      # 0-byte = crashed
    fails(py("build_defense_matrix.py", "--master", d / "master.tsv", "--df-dir", d / "df",
             "--out", d / "x.tsv", check=False))
    write_tsv(d / "df" / "g2" / "g2_defense_finder_systems.tsv", ["sys_id", "kind"], [["a", "RM"]])
    fails(py("build_defense_matrix.py", "--master", d / "master.tsv", "--df-dir", d / "df",
             "--out", d / "x.tsv", check=False), "no 'type' column")
    write_tsv(d / "dup.tsv", ["genome_id"], [["g3"], ["g3"]])
    fails(py("build_defense_matrix.py", "--master", d / "dup.tsv", "--df-dir", d / "df",
             "--out", d / "x.tsv", check=False), "duplicated")
    write_tsv(d / "empty.tsv", ["genome_id", "pass_primary"], [])
    py("build_defense_matrix.py", "--master", d / "empty.tsv", "--df-dir", d / "df", "--out", d / "e.tsv")
    assert open(d / "e.tsv").read() == "genome_id\tPrometheus\n"


# --------------------------------------------------------------------------- #
# association_tests.R (real phylolm / phytools / caper runs)
# --------------------------------------------------------------------------- #
ASSOC_ARGS = ["--boot", "20", "--d-permutations", "199", "--seed", "7", "--threads", "2"]


def run_assoc(sd: Path, out: Path, tiers="A", extra=(), check=True, tree=None, master=None):
    return rscript("association_tests.R", "--master", master or sd / "master.tsv",
                   "--tree", tree or sd / "tree.nwk", "--primary-tiers", tiers,
                   *ASSOC_ARGS, *extra, "--out", out, check=check)


def verify(sd: Path, tiers="A", seed=7, dperm=199) -> dict:
    out = sd / f"verify_{tiers}.json"
    run([RSCRIPT, str(VERIFY), "assoc", str(sd / "master.tsv"), str(sd / "tree.nwk"), tiers,
         str(seed), str(dperm), str(out)])
    return json.load(open(out))


@functools.lru_cache(maxsize=None)
def assoc_signal(tiers="A"):
    sd = sim("signal", 100, 21)
    out = sd / f"assoc_{tiers}.tsv"
    run_assoc(sd, out, tiers, extra=() if tiers == "A" else ("--boot", "5"))
    return sd, read_tsv(out)


@suite.test
def test_association_signal_against_direct_calls():
    sd, a = assoc_signal("A")
    assert list(a.columns) == ASSOC_COLS and list(a.test) == ASSOC_ORDER, (list(a.columns), list(a.test))
    ref = verify(sd, "A")
    use, excl = expected_sets(sd / "master.tsv", sd / "tree.nwk", ["A"])
    assert ref["n"] == len(use) and ref["n_excluded"] == len(excl) > 0
    r = a.set_index("test")
    for t in ["naive_glm", "naive_fisher", "phyloglm_methanogen", "fitPagel", "D_statistic"]:
        assert int(r.loc[t, "n"]) == len(use), (t, r.loc[t, "n"], len(use))
        assert f"n_excluded_tier={len(excl)}" in r.loc[t, "note"], r.loc[t, "note"]
    close(r.loc["naive_glm", "estimate"], ref["glm_est"], 1e-6, "glm log-odds")
    close(r.loc["naive_glm", "p_value"], ref["glm_lrt_p"], 1e-9, "glm LRT p")
    close(r.loc["naive_fisher", "estimate"], ref["fisher_logor"], 1e-6, "Fisher log OR")
    close(r.loc["naive_fisher", "p_value"], ref["fisher_p"], 1e-12, "Fisher p")
    # Pagel: script (fitMk 'pruning') == direct phytools::fitPagel (default 'lik' engine)
    close(r.loc["fitPagel", "statistic"], ref["pagel_lr"], 1e-3, "Pagel LR")
    close(r.loc["fitPagel", "p_value"], ref["pagel_p"], 1e-4, "Pagel p")
    # independent model factorises into two 2-state chains (guard used by pagel_test)
    close(ref["sep_logL_sum"], ref["pagel_indep_logL"], 1e-3, "independent logL factorisation")
    # D: same seed and permutation count -> identical estimate; (k+1)/(N+1) p-value
    close(r.loc["D_statistic", "estimate"], ref["D"], 1e-9, "D estimate")
    close(r.loc["D_statistic", "p_value"], ref["D_p_random"], 1e-12, "D p (vs random)")
    assert abs(fnum(r.loc["D_statistic", "p_value"]) * 200 - round(fnum(r.loc["D_statistic", "p_value"]) * 200)) < 1e-9
    # phyloglm: multi-start can only improve on the default single fit's penalised logLik
    import re
    pen = float(re.search(r"penlogLik (-?[0-9.eE+-]+)", r.loc["phyloglm_methanogen", "note"]).group(1))
    assert pen >= ref["pglm_default_pen"] - 1e-6, (pen, ref["pglm_default_pen"])
    if ref["pglm_default_conv"] == 0 and abs(pen - ref["pglm_default_pen"]) < 1e-3:
        close(r.loc["phyloglm_methanogen", "estimate"], ref["pglm_default_est"], 0.05, "phyloglm estimate")
    # known generative signal (beta_methanogen = 1.5 > 0): detected by the phylogenetic models
    for t in ["phyloglm_methanogen", "phyloglm_confounders"]:
        est, lo, p = fnum(r.loc[t, "estimate"]), fnum(r.loc[t, "ci_low"]), fnum(r.loc[t, "p_value"])
        assert est > 0 and p < 0.05 and lo < est < fnum(r.loc[t, "ci_high"]), (t, est, lo, p)
        assert "replicates fitted" in r.loc[t, "note"] and "best of 4 starts" in r.loc[t, "note"]
    assert r.loc["phyloglm_confounders", "note"].startswith("HEADLINE")


@suite.test
def test_association_B_only_excluded_under_A_present_under_AB():
    sd, a = assoc_signal("A")
    _, ab = assoc_signal("A,B")
    use_a, excl_a = expected_sets(sd / "master.tsv", sd / "tree.nwk", ["A"])
    use_ab, excl_ab = expected_sets(sd / "master.tsv", sd / "tree.nwk", ["A", "B"])
    b_only = set(use_ab[use_ab.prometheus_tier == "B"].genome_id)
    assert b_only and not (b_only & set(use_a.genome_id)), "B-only genomes must be excluded under tier A"
    n_a = int(a.set_index("test").loc["fitPagel", "n"]); n_ab = int(ab.set_index("test").loc["fitPagel", "n"])
    assert (n_a, n_ab) == (len(use_a), len(use_ab)) and n_ab - n_a == len(b_only)
    n_present_ab = int(use_ab.prometheus_tier.isin(["A", "B"]).sum())
    assert f"n_present={n_present_ab}" in ab.set_index("test").loc["fitPagel", "note"]
    assert f"n_excluded_tier={len(excl_ab)}" in ab.set_index("test").loc["naive_glm", "note"]


@suite.test
def test_association_clumped_D_and_node_labels():
    sd = sim("clumped", 80, 3)
    out = sd / "assoc.tsv"
    m = read_tsv(sd / "master.tsv").assign(methanogen="False")   # only D is estimable -> fast
    m.to_csv(sd / "master_dOnly.tsv", sep="\t", index=False)
    # duplicated quoted '100/100' support labels must not break caper::comparative.data
    run_assoc(sd, out, master=sd / "master_dOnly.tsv", extra=("--boot", "0"))
    r = read_tsv(out).set_index("test")
    D, p = fnum(r.loc["D_statistic", "estimate"]), fnum(r.loc["D_statistic", "p_value"])
    assert D < 0.3 and p <= 0.01, (D, p)


@suite.test
def test_association_empty_not_estimable_and_failures():
    sd = sim("signal", 100, 21)
    d = tmpdir("dist_assoc_bad_")
    (d / "empty.nwk").write_text("();\n")
    run_assoc(sd, d / "e.tsv", tree=d / "empty.nwk")
    assert open(d / "e.tsv").read().rstrip("\n").split("\t") == ASSOC_COLS
    m = read_tsv(sd / "master.tsv")
    m.assign(methanogen="False").to_csv(d / "nometh.tsv", sep="\t", index=False)
    run_assoc(sd, d / "nm.tsv", master=d / "nometh.tsv", extra=("--boot", "0"))
    r = read_tsv(d / "nm.tsv").set_index("test")
    assert list(r.index) == ASSOC_ORDER
    for t in ASSOC_ORDER[:5]:
        assert r.loc[t, "p_value"] == "NA" and "no variation in methanogen" in r.loc[t, "note"]
    assert fnum(r.loc["D_statistic", "p_value"]) < 1, "D needs only presence variation"
    m.assign(prometheus_tier="", prometheus_present="False").to_csv(d / "noprom.tsv", sep="\t", index=False)
    run_assoc(sd, d / "np.tsv", master=d / "noprom.tsv", extra=("--boot", "0"))
    r = read_tsv(d / "np.tsv")
    assert (r.p_value == "NA").all() and r.note.str.contains("no variation in Prometheus").all()
    # failures: unrooted tree, missing column, inconsistent master, no branch lengths
    (d / "unrooted.nwk").write_text("(000001:1,000002:1,(000003:1,000004:1):1);\n")
    fails(run_assoc(sd, d / "x.tsv", tree=d / "unrooted.nwk", check=False), "not rooted")
    m.drop(columns=["n_defense_systems_excl_prometheus"]).to_csv(d / "nocol.tsv", sep="\t", index=False)
    fails(run_assoc(sd, d / "x.tsv", master=d / "nocol.tsv", check=False), "n_defense_systems_excl_prometheus")
    bad = m.copy(); bad.loc[bad.prometheus_tier == "A", "prometheus_present"] = "False"
    bad.to_csv(d / "incons.tsv", sep="\t", index=False)
    fails(run_assoc(sd, d / "x.tsv", master=d / "incons.tsv", check=False), "disagrees")
    (d / "nobl.nwk").write_text("((000001,000002),(000003,000004));\n")
    fails(run_assoc(sd, d / "x.tsv", tree=d / "nobl.nwk", check=False), "branch lengths")
    fails(run_assoc(sd, d / "x.tsv", tiers="Z", check=False), "subset of a,b,c")


@suite.test
def test_association_complete_separation():
    """methanogen == presence: Wald GLM is meaningless (old output: 25.4 +/- 7475, p=0.997)."""
    sd = sim("signal", 100, 21)
    d = tmpdir("dist_assoc_sep_")
    m = read_tsv(sd / "master.tsv")
    m["methanogen"] = (m.prometheus_tier == "A").map({True: "True", False: "False"})
    m.to_csv(d / "sep.tsv", sep="\t", index=False)
    run_assoc(sd, d / "s.tsv", master=d / "sep.tsv", extra=["--boot", "0"])
    r = read_tsv(d / "s.tsv").set_index("test")
    assert r.loc["naive_glm", "estimate"] == "NA" and "separation" in r.loc["naive_glm", "note"]
    assert fnum(r.loc["naive_glm", "p_value"]) < 1e-6, "LRT p stays valid under separation"
    assert r.loc["naive_fisher", "estimate"] == "NA" and fnum(r.loc["naive_fisher", "ci_low"]) > 0
    assert fnum(r.loc["naive_fisher", "p_value"]) < 1e-6
    est = fnum(r.loc["phyloglm_methanogen", "estimate"])
    assert math.isfinite(est) and est > 0, "Firth-penalised MPLE stays finite under separation"
    assert "Wald CI" in r.loc["phyloglm_methanogen", "note"]


@suite.test
def test_association_zero_length_branch():
    sd = sim("clumped", 60, 5)
    d = tmpdir("dist_assoc_zero_")
    txt = open(sd / "tree.nwk").read()
    import re
    m = re.search(r"(\(\s*'?(\d{6})'?:)([0-9.eE-]+)(,\s*'?(\d{6})'?:)([0-9.eE-]+)\)", txt)
    assert m, "no cherry found"
    txt = txt[:m.start(3)] + "0" + txt[m.end(3):m.start(6)] + "0" + txt[m.end(6):]   # zero-length cherry
    (d / "zero.nwk").write_text(txt)
    read_tsv(sd / "master.tsv").assign(methanogen="False").to_csv(d / "m.tsv", sep="\t", index=False)
    r = run_assoc(sd, d / "z.tsv", tree=d / "zero.nwk", master=d / "m.tsv", extra=("--boot", "0"))
    assert "zero-length branch" in r.stderr
    t = read_tsv(d / "z.tsv").set_index("test")
    assert fnum(t.loc["D_statistic", "estimate"]) == fnum(t.loc["D_statistic", "estimate"])  # not NaN


@suite.test
def test_association_reproducible():
    """Same inputs and seed, --threads 2 (cached run) vs --threads 1: byte-identical output,
    including the parametric-bootstrap refits (5 replicates) and the seeded IG10 starts."""
    sd, _ = assoc_signal("A,B")
    d = tmpdir("dist_assoc_rep_")
    run_assoc(sd, d / "t1.tsv", "A,B", extra=("--boot", "5", "--threads", "1"))
    assert open(d / "t1.tsv").read() == open(sd / "assoc_A,B.tsv").read(), "not reproducible across thread counts"
    assert "replicates fitted" in open(d / "t1.tsv").read() or "bootstrap CI not reported" in open(d / "t1.tsv").read()


@suite.test
def test_pagel_guard_helpers():
    sd = sim("signal", 60, 5)
    code = f"""
suppressPackageStartupMessages({{library(ape); library(phytools)}})
source('{SCRIPTS_R / "association_tests.R"}')
m <- read.delim('{sd / "master.tsv"}', colClasses = 'character', na.strings = character(0))
tr <- prepare_tree(read_species_tree('{sd / "tree.nwk"}'))
k <- m$pass_primary == 'True' & m$tree_tip %in% tr$tip.label & m$prometheus_tier %in% c('', 'A')
tr <- keep.tip(tr, m$tree_tip[k]); m <- m[k, ][match(tr$tip.label, m$tree_tip[k]), ]
x <- setNames(as.integer(m$methanogen == 'True'), tr$tip.label)
y <- setNames(as.integer(m$prometheus_tier == 'A'), tr$tip.label)
pg <- pagel_test(tr, x, y)
fp <- fitPagel(tr, factor(x), factor(y), lik.func = 'pruning')
ld2 <- pagel_refit_dependent(tr, factor(x, levels = c(0, 1)), factor(y, levels = c(0, 1)), fp$independent.Q)
stopifnot(abs(pg$lr - fp$lik.ratio) < 1e-3, pg$lr >= 0, ld2 >= as.numeric(fp$independent.logL) - 1e-6)
stopifnot(identical(tier_response(c('', 'A', 'B', 'C'), 'A'), c(0L, 1L, NA, NA)))
stopifnot(identical(tier_response(c('', 'A', 'B', 'C'), c('A', 'B')), c(0L, 1L, 1L, NA)))
cat('GUARD_OK', pg$lr, ld2, as.numeric(fp$independent.logL), '\\n')
"""
    r = run([RSCRIPT, "-e", code])
    assert "GUARD_OK" in r.stdout, r.stdout


# --------------------------------------------------------------------------- #
# cooccurrence.R
# --------------------------------------------------------------------------- #
def bh(p):
    import numpy as np
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    adj = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(adj, 1); return out


@functools.lru_cache(maxsize=None)
def cooc_run(threads: int):
    sd = sim("signal", 60, 5)
    out = sd / f"cooc_t{threads}.tsv"
    rscript("cooccurrence.R", "--matrix", sd / "matrix.tsv", "--master", sd / "master.tsv",
            "--tree", sd / "tree.nwk", "--fdr", "BH", "--seed", "11", "--primary-tiers", "A",
            "--threads", threads, "--out", out)
    return sd, read_tsv(out)


@suite.test
def test_cooccurrence_known_answer_against_direct():
    sd, c = cooc_run(2)
    assert list(c.columns[:10]) == COOC_BASE
    assert list(c.system_b) == ["CoSys", "ExSys", "RandSys"], "RareSys (2 genomes) must be outside the window"
    use, excl = expected_sets(sd / "master.tsv", sd / "tree.nwk", ["A"])
    n_a = int((use.prometheus_tier == "A").sum())
    assert (c.n.astype(int) == len(use)).all() and (c.n_a.astype(int) == n_a).all()
    assert c.note.str.contains(f"n_excluded_tier={len(excl)}").all()
    r = c.set_index("system_b")
    assert r.loc["CoSys", "direction"].startswith("co-occurrence") and fnum(r.loc["CoSys", "p_adj"]) < 0.05
    assert r.loc["ExSys", "direction"].startswith("exclusion")
    assert fnum(r.loc["CoSys", "odds_ratio"]) > 1 and fnum(r.loc["ExSys", "odds_ratio"]) < 1
    import numpy as np
    np.testing.assert_allclose(c.p_adj.astype(float), bh(c.p_value.astype(float)), rtol=1e-10)
    out = sd / "verify_cooc.json"
    run([RSCRIPT, str(VERIFY), "cooc", str(sd / "master.tsv"), str(sd / "tree.nwk"), "A",
         str(sd / "matrix.tsv"), "CoSys,ExSys,RandSys", str(out)])
    ref = json.load(open(out))
    for s in ["CoSys", "ExSys", "RandSys"]:
        close(r.loc[s, "statistic"], ref["systems"][s]["lr"], 1e-3, f"{s} Pagel LR")
        close(r.loc[s, "p_value"], ref["systems"][s]["p"], 1e-4, f"{s} Pagel p")
        assert int(r.loc[s, "n_b"]) == ref["systems"][s]["n_b"]
    assert (c.method == "fitPagel").all(), "no silent fallback to another test"
    assert c.phyloglm_logodds.map(fnum).notna().all(), "genome-size-adjusted phyloglm must be fitted"


@suite.test
def test_cooccurrence_thread_count_does_not_change_results():
    _, c1 = cooc_run(1)
    _, c2 = cooc_run(2)
    assert c1.equals(c2), "results differ between --threads 1 and 2"


@suite.test
def test_cooccurrence_bacterial_cli_empty_and_failures():
    sd = sim("signal", 60, 5)
    d = tmpdir("dist_cooc_bad_")
    # bact_cooccurrence CLI (no new flags), slim master where tree_tip = genome_id for every genome
    m = read_tsv(sd / "master.tsv")
    slim = m[["genome_id", "pass_primary", "prometheus_present", "prometheus_tier", "genome_size"]].assign(tree_tip=m.genome_id)
    slim.to_csv(d / "bm.tsv", sep="\t", index=False)
    read_tsv(sd / "matrix.tsv")[["genome_id", "Prometheus", "CoSys"]].to_csv(d / "m1.tsv", sep="\t", index=False)
    rscript("cooccurrence.R", "--matrix", d / "m1.tsv", "--master", d / "bm.tsv",
            "--tree", sd / "tree.nwk", "--fdr", "BH", "--seed", "12345", "--out", d / "b.tsv")
    b = read_tsv(d / "b.tsv")
    assert list(b.system_b) == ["CoSys"] and (b.method == "fitPagel").all()
    # genuinely empty inputs -> header-only
    (d / "empty.nwk").write_text("();\n")
    rscript("cooccurrence.R", "--matrix", d / "m1.tsv", "--master", sd / "master.tsv",
            "--tree", d / "empty.nwk", "--out", d / "e1.tsv")
    assert open(d / "e1.tsv").read().split("\n")[0].split("\t")[:10] == COOC_BASE
    assert open(d / "e1.tsv").read().count("\n") == 1
    m.assign(prometheus_tier="", prometheus_present="False").to_csv(d / "noprom.tsv", sep="\t", index=False)
    rscript("cooccurrence.R", "--matrix", d / "m1.tsv", "--master", d / "noprom.tsv",
            "--tree", sd / "tree.nwk", "--out", d / "e2.tsv")
    assert open(d / "e2.tsv").read().count("\n") == 1
    # failures
    mat = read_tsv(sd / "matrix.tsv"); bad = mat.copy(); bad.loc[0, "CoSys"] = "2"
    bad.to_csv(d / "bad.tsv", sep="\t", index=False)
    fails(rscript("cooccurrence.R", "--matrix", d / "bad.tsv", "--master", sd / "master.tsv",
                  "--tree", sd / "tree.nwk", "--out", d / "x.tsv", check=False), "not strictly 0/1")
    fails(rscript("cooccurrence.R", "--matrix", sd / "matrix.tsv", "--master", sd / "master.tsv",
                  "--tree", sd / "tree.nwk", "--fdr", "xyz", "--out", d / "x.tsv", check=False), "--fdr")
    m.iloc[1:].to_csv(d / "short.tsv", sep="\t", index=False)
    fails(rscript("cooccurrence.R", "--matrix", sd / "matrix.tsv", "--master", d / "short.tsv",
                  "--tree", sd / "tree.nwk", "--out", d / "x.tsv", check=False), "absent from the master")
    (d / "unrooted.nwk").write_text("(000001:1,000002:1,(000003:1,000004:1):1);\n")
    fails(rscript("cooccurrence.R", "--matrix", sd / "matrix.tsv", "--master", sd / "master.tsv",
                  "--tree", d / "unrooted.nwk", "--out", d / "x.tsv", check=False), "not rooted")


# --------------------------------------------------------------------------- #
# plot_distribution.py / plot_distribution_extra.py
# --------------------------------------------------------------------------- #
@suite.test
def test_plots_real_inputs_placeholders_and_helpers():
    sd, a = assoc_signal("A")
    cd, _c = cooc_run(2)
    d = tmpdir("dist_plots_")
    py("compute_prevalence.py", "--master", sd / "master.tsv", "--primary-tiers", "A", "--out", d / "prev.tsv")
    py("plot_distribution.py", "--prevalence", d / "prev.tsv", "--tests", sd / "assoc_A.tsv",
       "--cooccurrence", cd / "cooc_t2.tsv", "--outdir", d / "p1")
    for stem in ["prevalence_by_taxon", "association_forest", "cooccurrence_heatmap"]:
        for ext in ["png", "svg"]:
            f = d / "p1" / f"{stem}.{ext}"
            assert f.exists() and f.stat().st_size > 1000, f
    py("plot_distribution_extra.py", "--prevalence", d / "prev.tsv", "--cooccurrence", cd / "cooc_t2.tsv",
       "--master", sd / "master.tsv", "--outdir", d / "p2", "--formats", "png")
    assert sorted(p.name for p in (d / "p2").iterdir()) == ["cooccurrence_oddsratio.png", "prevalence_caterpillar.png"]
    # header-only inputs -> labelled placeholders, exit 0
    (d / "h").mkdir()
    write_tsv(d / "h" / "prev.tsv", PREV_BASE, []); write_tsv(d / "h" / "a.tsv", ASSOC_COLS, [])
    write_tsv(d / "h" / "c.tsv", COOC_BASE, [])
    py("plot_distribution.py", "--prevalence", d / "h/prev.tsv", "--tests", d / "h/a.tsv",
       "--cooccurrence", d / "h/c.tsv", "--outdir", d / "p3")
    py("plot_distribution_extra.py", "--prevalence", d / "h/prev.tsv", "--cooccurrence", d / "h/c.tsv",
       "--master", sd / "master.tsv", "--outdir", d / "p3")
    assert len(list((d / "p3").glob("*.svg"))) == 5
    write_tsv(d / "h" / "badc.tsv", ["system_b", "p_value"], [["X", 0.1]])
    fails(py("plot_distribution.py", "--prevalence", d / "h/prev.tsv", "--tests", d / "h/a.tsv",
             "--cooccurrence", d / "h/badc.tsv", "--outdir", d / "p4", check=False), "missing required")
    sys.path.insert(0, str(SCRIPTS_PY))
    import plot_distribution as pdist
    import plot_distribution_extra as pext
    assert [pdist.direction_sign(s) for s in ["exclusion (OR=0.2)", "co-occurrence (OR=Inf)", "neutral (OR=1)",
                                              "ambiguous (OR=NA)"]] == [-1, 1, 0, 0]
    close(pext.log2_or_from_counts(10, 4, 4, 4), math.log2(4.5 * 6.5 / (0.5 * 0.5)), 1e-12, "Haldane log2 OR")
    close(pext.log2_or_from_counts(20, 6, 8, 3), math.log2(3 * 9 / (3 * 5)), 1e-12, "sample log2 OR")
    import pandas as pd
    prev = pd.DataFrame({"rank": ["gtdb_phylum"], "taxon": ["P"], "prevalence": [0.5], "tier_set": ["A"]})
    mast = pd.DataFrame({"pass_primary": ["True"] * 4 + ["False"], "prometheus_tier": ["A", "B", "", "", "A"]})
    close(pext.overall_prevalence(prev, mast), 1 / 3, 1e-12, "overall prevalence with B excluded")


if __name__ == "__main__":
    suite.main("DISTRIBUTION TESTS PASSED")
