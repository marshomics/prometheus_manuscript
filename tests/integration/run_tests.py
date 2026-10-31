#!/usr/bin/env python3
"""End-to-end integration test (node gate N2): the REAL Snakemake workflow, Phase 1
(`preprocess`) then Phase 2 (`all`), on the synthetic project of make_world.py.

Real tools run wherever they are installed here: Prodigal, DefenseFinder 3.0.0 (pinned
models 3.1.0 + CasFinder 3.1.0), HMMER, MAFFT, ClipKIT, IQ-TREE 2, MMseqs2, PAL2NAL,
HyPhy 2.5.62 and every in-house Python/R script. Stand-ins replace only the tools whose
reference data cannot be installed here (GTDB-Tk, CheckM2, BBMap statswrapper, geNomad,
CheckV); they reproduce the pinned versions' file layouts and id rules. GeneRax is
disabled (its 2.1.3 builds crash on this macOS host; see tests/evolution/REVIEW.md).
The workflow is copied with its `conda:` directives removed (tools come from PATH).

Checks: every stage consumes its producer's real output (the run itself) and the
SCHEMA.md column contract holds; the planted signal survives end to end; the review's
fixes hold in the assembled pipeline (GTDB-Tk id guard, real family gene tree, DUF3320
optional for tier A, tier rule columns, configured figure formats, provenance).
Prints INTEGRATION TESTS PASSED.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import REPO, PYTHON, RSCRIPT, TOOLS_BIN, Suite, tmpdir, read_tsv  # noqa: E402
sys.path.insert(0, str(HERE))
import make_world as mw  # noqa: E402

import yaml  # noqa: E402

suite = Suite("integration")
SMK = os.environ.get("PROM_SNAKEMAKE") or shutil.which("snakemake")
ENVS = Path(PYTHON).resolve().parents[1].parent          # .../conda/envs
STATE = {}


def df_models():
    p = os.environ.get("PROM_DF_MODELS")
    cands = [p] if p else []
    cands += [os.path.join(os.path.abspath(tempfile.gettempdir()),
                           "prometheus-tests-dfmodels-3.1.0_cas3.1.0"),
              "/tmp/prometheus-tests-dfmodels-3.1.0_cas3.1.0"]
    for c in cands:
        if c and all(os.path.isfile(os.path.join(c, k, "metadata.yml"))
                     for k in ("defense-finder-models", "CasFinder")):
            return os.path.abspath(c)
    return None


def tool(name, *envs):
    for e in envs:
        p = ENVS / e / "bin" / name
        if p.exists():
            return p
    if TOOLS_BIN and (Path(TOOLS_BIN) / name).exists():
        return Path(TOOLS_BIN) / name
    return None


def toolbin(root):
    """One directory with exactly the binaries the rules call (no env's python can shadow
    the analysis python)."""
    tb = root / "toolbin"
    tb.mkdir()
    links = {"python": Path(PYTHON), "python3": Path(PYTHON), "Rscript": Path(RSCRIPT)}
    for n in ("prodigal", "hmmscan", "hmmsearch", "hmmpress", "hmmbuild", "hmmemit", "hmmalign",
              "mafft", "iqtree2", "mmseqs", "pal2nal.pl", "seqkit", "perl"):
        p = tool(n, "prom-tools")
        if p:
            links[n] = p
    for n in ("clipkit", "famsa"):
        links[n] = tool(n, "prom-evolution-clipkit", "prom-tools")
    links["hyphy"] = tool("hyphy", "prom-evolution-hyphy", "prom-tools")
    for n in ("defense-finder", "macsyfinder", "macsydata"):
        links[n] = tool(n, "prom-detection_diversity-defensefinder")
    missing = [n for n, p in links.items() if not p]
    if missing:
        raise AssertionError(f"integration test needs these tools installed: {missing}")
    for n, p in links.items():
        os.symlink(p, tb / n)
    return tb


def copy_workflow(root):
    """Private copy of workflow/ with the conda: directives removed (tools come from PATH;
    Snakemake 8.30 otherwise calls `conda` for every rule with a directive)."""
    shutil.copytree(REPO / "workflow", root / "workflow",
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    pat2 = re.compile(r'^[ \t]+conda:[ \t]*\n[ \t]+f?"[^"]+\.yaml"[ \t]*\n', re.M)
    pat1 = re.compile(r'^[ \t]+conda:[ \t]*f?"[^"]+\.yaml"[ \t]*\n', re.M)
    for f in list((root / "workflow" / "rules").glob("*.smk")) + [root / "workflow" / "Snakefile"]:
        t = pat1.sub("", pat2.sub("", f.read_text()))
        assert not re.search(r"^[ \t]+conda:", t, re.M), f"conda directive left in {f.name}"
        f.write_text(t)


def write_config(root, models):
    c = yaml.safe_load(open(REPO / "config" / "config.yaml"))
    c.update({"samples": "results/preprocess/species_metadata.tsv",
              "metadata": "results/preprocess/species_metadata.tsv",
              "genomes_dir": "results/preprocess/genomes", "genome_extension": "fna",
              "results_dir": "results", "logs_dir": "logs", "seed": 12345})
    c["input"].update({"metadata": "data/input_metadata.tsv", "genomes_dir": "data/genomes",
                       "genome_extension": "fasta"})
    c["preprocess"]["gtdbtk_batch_size"] = 10                       # 17 genomes -> 2 batches
    c["db"] = {"checkm2": "db/checkm2.dmnd", "gtdbtk_data": "db/gtdbtk",
               "pfam_hmm": "db/pfam/Pfam-A.hmm", "kofam_ko_list": "db/kofam/ko_list",
               "kofam_profiles": "db/kofam/methanogenesis_profiles",
               "genomad": "db/genomad_db", "checkv": "db/checkv-db"}
    c["refs"] = {"prometheus_faa": "refs/prometheus_reference.faa",
                 "paralog_faa": "refs/paralog_reference.faa", "mcr_hmm_dir": "refs/mcr_hmm",
                 "mcrA_ref_aln": "", "mcrA_ref_labels": ""}
    c["prometheus"].update({"defensefinder_models_dir": models, "validation_tree_model": "LG+G4"})
    c["diversity"]["rarefaction_permutations"] = 20
    c["tree"].update({"scope": "focal", "pmsf": False, "iqtree_model": "LG+G4", "outgroup": [],
                      "outgroup_taxon": "p__Synth_O", "ufboot": 1000, "alrt": 1000})
    c["methanogen"]["anme_via_mcrA_tree"] = False
    c["methanogenesis"] = {
        "enabled": True,
        "rules": "workflow/resources/methanogenesis/rules/methanogenesis_rules.yaml",
        "readthrough": True, "absence_alpha": 0.05,
        "stats": {"min_n": 4, "min_class": 2, "min_changes": 2, "nsim": 20, "corhmm": False,
                  "epv_min": 5, "vif_max": 5, "complete_threshold": 95},
    }
    ev = c["evolution"]
    ev.update({"gene_tree_model": "LG+G4", "domain_tree_model": "LG+G4", "domain_min_taxa": 6})
    ev["reconciliation"]["run_generax"] = False
    ev["ancestral"].update({"nsim": 50, "nnull": 50, "prevalence_tol": 0.25})
    ev["selection"].update({"methods": ["FEL", "RELAX"], "gard": False})
    c["stats"].update({"phyloglm_boot": 10, "phylo_d_permutations": 99, "cooccurrence_min_count": 2})
    c["plotting"]["formats"] = ["png", "pdf"]                       # non-default on purpose
    c["threads"] = {"small": 1, "medium": 2, "large": 4}
    (root / "config").mkdir()
    with open(root / "config" / "config.yaml", "w") as fh:
        yaml.safe_dump(c, fh, sort_keys=False)


def smk(root, *targets, check=True, timeout=5400):
    home = root / "home"
    home.mkdir(exist_ok=True)
    e = {k: v for k, v in os.environ.items() if k in ("R_LIBS", "LANG", "LC_ALL", "USER", "LOGNAME")}
    e.update({"PATH": os.pathsep.join([str(root / "toolbin"), str(root / "fakebin"), TOOLS_BIN or "",
                                       "/usr/bin", "/bin", "/usr/sbin", "/sbin"]),
              "HOME": str(home), "XDG_CACHE_HOME": str(home / ".cache"), "MPLCONFIGDIR": str(home / ".mpl"),
              "TMPDIR": str(root / "tmp"), "PYTHONPATH": str(root / "fakepy"), "MPLBACKEND": "Agg",
              "FAKE_GTDBTK_TAX": str(root / "data" / "fake_tax.tsv"),
              "FAKE_CHECKM2_TABLE": str(root / "data" / "fake_checkm2.tsv"),
              "STUB_SPECIES_MSA": str(root / "data" / "species_msa.faa"),
              # provenance probes need an active conda env; any real env satisfies the rule
              "CONDA_PREFIX": str(Path(PYTHON).resolve().parents[1])})
    (root / "tmp").mkdir(exist_ok=True)
    cmd = [SMK, *targets, "--cores", "6", "--nolock", "--drop-metadata", "-s", "workflow/Snakefile",
           "--configfile", "config/config.yaml", "--show-failed-logs"]
    t0 = time.time()
    p = subprocess.run(cmd, cwd=root, env=e, text=True, capture_output=True, timeout=timeout)
    STATE.setdefault("timing", []).append((" ".join(targets), round(time.time() - t0)))
    if check and p.returncode != 0:
        raise AssertionError(f"snakemake {' '.join(targets)} failed ({p.returncode}):\n"
                             f"{p.stdout[-4000:]}\n--- stderr ---\n{p.stderr[-6000:]}")
    return p


def world():
    if "root" in STATE:
        return STATE["root"]
    if not SMK:
        raise AssertionError("snakemake not found (set PROM_SNAKEMAKE)")
    models = df_models()
    if not models:
        raise AssertionError("pinned DefenseFinder models not found (set PROM_DF_MODELS to a "
                             "macsydata --target dir with defense-finder-models 3.1.0 + CasFinder 3.1.0)")
    root = tmpdir("integration_")
    os.environ["PATH"] = os.pathsep.join([TOOLS_BIN or "", os.environ.get("PATH", "")])  # hmmbuild etc.
    STATE["truth"] = mw.build(root, PYTHON)
    toolbin(root)
    copy_workflow(root)
    write_config(root, models)
    STATE["root"] = root
    return root


def cols(path):
    with open(path) as fh:
        return fh.readline().rstrip("\n").split("\t")


def require_cols(path, required):
    have = cols(path)
    miss = [c for c in required if c not in have]
    assert not miss, f"{path.name}: missing SCHEMA columns {miss} (have {have[:20]}...)"


# --------------------------------------------------------------------------- #
NICHE_COLS = ("group_Marine", "niche_anoxic", "niche_thermophile")


def expected_species_niche():
    """{rep label: {col: '1'|'0'|''}} from make_world's own NICHE/TAX: the species value is 1
    if any genome of the species is 1, else 0 if any is 0, else '' (unknown is not absence)."""
    def species(lab):
        t = mw.TAX[lab]
        return t.split(";s__")[-1] if ";s__" in t else None
    out = {}
    for lab in mw.REPS:
        members = [m for m in mw.LABELS if species(m) and species(m) == species(lab)]
        out[lab] = {}
        for c in NICHE_COLS:
            vals = [mw.NICHE[m][c] for m in members]
            out[lab][c] = "1" if "1" in vals else "0" if "0" in vals else ""
    return out


@suite.test
def phase1_preprocess():
    root = world()
    smk(root, "preprocess")
    R = root / "results" / "preprocess"
    acc, truth = mw.ACC, STATE["truth"]
    sp = read_tsv(R / "species_metadata.tsv")
    require_cols(R / "species_metadata.tsv", [
        "genome_id", "filename", "source", "classification", "gtdb_domain", "gtdb_phylum",
        "gtdb_species", "assembly_path", "completeness", "contamination", "gc", "genome_size",
        "n_contigs", "n50", "quality_score", "n_genomes_in_species", "n_passing_in_species",
        "group_Marine", "niche_anoxic", "niche_thermophile"])
    assert sorted(sp["genome_id"]) == sorted(truth["reps"]), sorted(sp["genome_id"])
    m1 = sp.set_index("genome_id").loc[acc["M1"]]
    assert m1["n_genomes_in_species"] == "2" and m1["niche_thermophile"] == "1", dict(m1)
    spi = sp.set_index("genome_id")
    for lab, want in expected_species_niche().items():
        for c, w in want.items():
            assert spi.loc[acc[lab], c] == w, (lab, c, spi.loc[acc[lab], c], w)
    rep = read_tsv(R / "filter_report.tsv").set_index("genome_id")
    st = rep["status"].to_dict()
    assert st[acc["M1"]] == "representative" and st[acc["M1b"]] != "representative", st
    assert st[acc["Q1"]] == "failed_quality" and st[acc["B1"]] == "excluded_domain", st
    assert st[acc["U1"]] != "representative", st
    for g in truth["reps"]:
        assert (root / "results" / "preprocess" / "genomes" / f"{g}.fna").exists(), g


@suite.test
def phase2_all_runs():
    root = world()
    if not (root / "results" / "preprocess" / "species_metadata.tsv").exists():
        smk(root, "preprocess")
    smk(root, "all")


@suite.test
def gene_calls_and_detection_contract():
    root = world()
    R = root / "results"
    for g in STATE["truth"]["reps"]:
        gt = R / "genes" / g / f"{g}.genes.tsv"
        require_cols(gt, ["protein_id", "contig", "contig_index", "gene_index", "start", "end",
                          "strand", "partial", "contig_length", "length_aa"])
        t = read_tsv(gt)
        assert t["protein_id"].str.fullmatch(re.escape(g) + r"_[12]_\d+").all(), t["protein_id"][:3]
    assert (R / "defense" / ".defensefinder_ok").exists()
    require_cols(R / "prometheus" / "candidate_sources.tsv",
                 ["protein_id", "genome_id", "detect_defensefinder", "df_score", "df_subtype"])
    cand = read_tsv(R / "prometheus" / "candidate_sources.tsv")
    called = {g for g in STATE["truth"]["df_called"] if g in STATE["truth"]["reps"]}
    assert set(cand["genome_id"]) == called, (sorted(set(cand["genome_id"])), sorted(called))
    assert mw.ACC["M4"] not in set(cand["genome_id"])     # split-domain homologue, below GA2


@suite.test
def catalog_tiers_and_duf3320_optional():
    root = world()
    p = root / "results" / "prometheus" / "prometheus_catalog.tsv"
    require_cols(p, ["locus_id", "genome_id", "protein_id", "length_aa", "detect_defensefinder",
                     "architecture", "n_domains", "catalytic_status", "defense_context",
                     "tree_assignment", "fragmented", "confidence_tier", "reason",
                     "missing_required", "canonical_order", "has_DUF3320", "placement_support"])
    cat = read_tsv(p).set_index("genome_id")
    assert not (cat["tree_assignment"] == "Paralog").any(), cat["tree_assignment"].to_dict()
    assert (cat["tree_assignment"] == "Prometheus").sum() >= 6, cat["tree_assignment"].to_dict()
    n2 = cat.loc[mw.ACC["N2"]]            # the DUF3320-less homologue
    assert n2["has_DUF3320"] in ("False", "0") and n2["missing_required"] == "", dict(n2)
    assert "DUF3320" not in n2["reason"], n2["reason"]
    four = [mw.ACC[k] for k in ("M1", "M2", "M3", "M5", "N1")]   # M4: not called by DefenseFinder
    tiers = cat.loc[four, "confidence_tier"].to_dict()
    assert set(tiers.values()) <= {"A", "B"}, tiers
    if n2["tree_assignment"] == "Prometheus" and n2["catalytic_status"] == "intact":
        assert n2["confidence_tier"] == "A", dict(n2)   # DUF3320 absence alone never blocks A
    STATE["tiers"] = cat["confidence_tier"].to_dict()


@suite.test
def scaffolds_methanogenesis_pathways():
    """Planted pathway capacity is recovered end to end: KO search -> pathway calls -> the
    master table's mgen_* columns and composite subtype."""
    root = world()
    M = root / "results" / "scaffolds" / "metabolism"
    truth = STATE["truth"]["mgen_pathways"]
    sub_truth = STATE["truth"]["mgen_subtype"]
    hits = read_tsv(M / "methanogenesis_ko_hits.tsv")
    assert set(hits.columns) >= {"genome_id", "protein_id", "ko", "call", "contig_index",
                                 "gene_index", "start", "end", "strand"}
    planted = STATE["truth"]["mgen_planted_kos"]
    for gid, kos in planted.items():
        if gid not in set(STATE["truth"]["reps"]):
            continue
        got = set(hits[hits.genome_id == gid].ko)
        assert set(kos) <= got, (gid, sorted(set(kos) - got))
    assert not set(hits.genome_id) - set(STATE["truth"]["methanogens"]), \
        "only methanogens are searched"

    pw = read_tsv(M / "methanogenesis_pathways.tsv").set_index("genome_id")
    for gid, want in truth.items():
        if gid not in pw.index:
            continue
        for path, exp in want.items():
            assert pw.loc[gid, path] == exp, (gid, path, pw.loc[gid, path], "expected", exp)
        assert pw.loc[gid, "methanogenesis_subtype"] == sub_truth[gid], \
            (gid, pw.loc[gid, "methanogenesis_subtype"])
    non_meth = [g for g in pw.index if g not in set(STATE["truth"]["methanogens"])]
    assert non_meth, "the world has non-methanogens"
    for g in non_meth:
        assert pw.loc[g, "co2_reducing"] == "not_applicable", (g, pw.loc[g].to_dict())
        assert pw.loc[g, "methanogenesis_subtype"] == "NA"

    m = read_tsv(root / "results" / "scaffolds" / "master_table.tsv").set_index("genome_id")
    for col in ["mgen_co2_reducing", "mgen_acetoclastic", "mgen_methylotrophic",
                "mgen_methyl_dismutating", "mgen_methyl_reducing", "subtype_status",
                "subtype_taxonomy_concordance"]:
        assert col in m.columns, col
    code = {"present": "1", "absent": "0"}
    for gid, want in truth.items():
        if gid not in m.index:
            continue
        for path, exp in want.items():
            assert m.loc[gid, f"mgen_{path}"] == code[exp], (gid, path, m.loc[gid, f"mgen_{path}"])
        assert m.loc[gid, "methanogenesis_subtype"] == sub_truth[gid]
    for g in non_meth:
        if g in m.index:
            assert m.loc[g, "mgen_co2_reducing"] == "NA", (g, m.loc[g, "mgen_co2_reducing"])


@suite.test
def scaffolds_master_and_species_tree():
    root = world()
    R = root / "results" / "scaffolds"
    rm = (R / "species_tree" / "rooting_method.txt").read_text()
    assert rm.startswith("method=outgroup"), rm[:200]
    tree = (R / "species_tree.treefile").read_text().strip()
    tips = set(re.findall(r"[(,]([A-Za-z0-9_.\-]+):", tree))
    assert tips == set(STATE["truth"]["reps"]), sorted(tips)
    gd = (R / "gtdbtk_align" / "gtdbtk_dropped.tsv").read_text()
    assert "user_" not in gd
    require_cols(R / "master_table.tsv", [
        "genome_id", "representative", "tree_tip", "gtdb_phylum", "completeness", "pass_primary",
        "genome_size", "mcrA", "methanogen", "is_anme", "prometheus_present", "prometheus_tier",
        "prometheus_n_loci", "n_defense_systems", "n_defense_systems_excl_prometheus",
        "group_Marine", "niche_anoxic", "niche_thermophile", "prometheus_status_primary",
        "prometheus_status_sensitivity", "mcr_group"])
    m = read_tsv(R / "master_table.tsv").set_index("genome_id")
    assert len(m) == 13
    meth = {g for g, v in m["methanogen"].items() if v in ("True", "1")}
    assert meth == set(STATE["truth"]["methanogens"]) & set(STATE["truth"]["reps"]), sorted(meth)
    assert set(m["prometheus_status_primary"]) <= {"present", "absent", "excluded"}
    for g, t in STATE.get("tiers", {}).items():
        st = m.loc[g, "prometheus_status_primary"]
        assert st == ("present" if t == "A" else "excluded" if t in ("B", "C") else st), (g, t, st)
    exp = expected_species_niche()
    for lab, want in exp.items():
        for c, w in want.items():
            assert m.loc[mw.ACC[lab], c] == w, (lab, c, m.loc[mw.ACC[lab], c], w)
    miss = read_tsv(R / "ecology_missingness.tsv").set_index("column")
    for c in NICHE_COLS:
        vals = [exp[lab][c] for lab in exp]
        got = tuple(miss.loc[c, k] for k in ("n_genomes", "n_1", "n_0", "n_NA"))
        want = tuple(str(x) for x in (len(vals), vals.count("1"), vals.count("0"), vals.count("")))
        assert got == want, (c, got, want)
    assert int(miss.loc["niche_thermophile", "n_NA"]) > 0      # the world has unknowns


@suite.test
def distribution_contract():
    root = world()
    R = root / "results" / "distribution"
    require_cols(R / "association_tests.tsv", ["test", "model", "term", "estimate", "std_error",
                                               "statistic", "p_value", "ci_low", "ci_high", "n", "note"])
    a = read_tsv(R / "association_tests.tsv")
    assert set(a["test"]) == {"naive_glm", "naive_fisher", "phyloglm_methanogen",
                              "phyloglm_confounders", "fitPagel", "D_statistic"}, set(a["test"])
    require_cols(R / "prevalence_by_taxon.tsv", ["rank", "taxon", "n_genomes", "n_present", "prevalence",
                                                 "ci_low", "ci_high", "tier_set", "n_excluded"])
    pv = read_tsv(R / "prevalence_by_taxon.tsv")
    assert (pv.iloc[0]["rank"], pv.iloc[0]["taxon"]) == ("overall", "all")
    require_cols(R / "cooccurrence.tsv", ["system_a", "system_b", "n_a", "n_b", "n_both", "method",
                                          "statistic", "p_value", "p_adj", "direction"])
    require_cols(R / "defense_matrix.tsv", ["genome_id"])


@suite.test
def evolution_contract_and_real_gene_tree():
    root = world()
    R = root / "results" / "evolution"
    tree = (R / "tree" / "family.treefile").read_text().strip()
    n_tips = len(re.findall(r"[(,]([^(),:;]+):", tree))
    cat = read_tsv(root / "results" / "prometheus" / "prometheus_catalog.tsv")
    n_valid = int((cat["confidence_tier"] != "Rejected").sum())
    assert tree != "();" and n_valid >= 4 and n_tips == n_valid, (n_tips, n_valid, tree[:200])  # EVO-34
    assert len((R / "tree" / "family.ufboot").read_text().splitlines()) == 1000
    assert (R / "align" / "family.codon.fasta").stat().st_size > 0
    require_cols(R / "align" / "family.codon_report.tsv", ["id", "status"])
    ss = read_tsv(R / "selection" / "selection_summary.tsv")
    require_cols(R / "selection" / "selection_summary.tsv",
                 ["method", "statistic", "value", "p_value", "n_sites_or_branches", "note"])
    relax = ss[(ss["method"] == "RELAX") & (ss["statistic"] == "K")]
    assert len(relax) == 1 and not relax.iloc[0]["note"].startswith("not run"), ss.to_dict("records")
    require_cols(R / "selection" / "site_table.tsv", ["site", "alignment_col", "domain",
                                                      "beta_alpha_or_dnds", "p_value", "classification"])
    rec = read_tsv(R / "reconciliation" / "reconciliation_events.tsv")
    assert rec["note"].str.contains("disabled").all()
    gl = read_tsv(R / "ancestral" / "gains_losses.tsv")
    require_cols(R / "ancestral" / "gains_losses.tsv", ["metric", "value", "note", "tier_set"])
    assert "enrichment_p" in set(gl["metric"]) and "n_tips_excluded_tier" in set(gl["metric"])
    require_cols(R / "domains" / "rf_distances.tsv", ["domain", "n_tips", "rf_to_family",
                                                     "normalized_rf", "note"])


@suite.test
def ecology_contract():
    root = world()
    R = root / "results" / "ecology"
    require_cols(R / "viral" / "viral_summary.tsv", ["genome_id", "n_proviruses", "total_provirus_kb",
                                                    "mean_checkv_completeness", "n_rloop_proviruses",
                                                    "frac_rloop_proviruses"])
    v = read_tsv(R / "viral" / "viral_summary.tsv")
    assert len(v) == 13 and set(v["n_proviruses"]) <= {"0", "1", "2"}, v["n_proviruses"].tolist()
    assert (v.loc[v["n_proviruses"] == "0", "frac_rloop_proviruses"] == "").all()
    require_cols(R / "ecology_associations.tsv", ["predictor", "level", "model", "estimate",
                                                  "std_error", "p_value", "p_adj", "n", "note"])
    e = read_tsv(R / "ecology_associations.tsv")
    assert {"group_Marine", "niche_anoxic"} <= set(e["predictor"]), set(e["predictor"])
    # n = species with a KNOWN niche value and a primary response (present / absent); unknown
    # niches are dropped (never fitted as 0) and tier-excluded genomes have no response
    m = read_tsv(root / "results" / "scaffolds" / "master_table.tsv")
    usable = m[(m["niche_thermophile"] != "") & m["prometheus_status_primary"].isin(["present", "absent"])]
    n_all = int((m["prometheus_status_primary"].isin(["present", "absent"])).sum())
    assert len(usable) < n_all, (len(usable), n_all)   # the world's unknowns must matter
    for n in e.loc[e["predictor"] == "niche_thermophile", "n"]:
        assert n in ("", str(len(usable))), (n, len(usable), n_all)
    assert (R / "ecology_associations_sensitivity.tsv").exists()


@suite.test
def registry_confirmatory_rows():
    root = world()
    p = root / "results" / "synthesis" / "results_registry.tsv"
    require_cols(p, ["family", "test", "term", "estimate", "ci_low", "ci_high", "p_value", "p_adj",
                     "hypothesis", "role"])
    r = read_tsv(p)
    conf = {(f, t, m) for f, t, m, h in zip(r["family"], r["test"], r["term"], r["hypothesis"])
            if h == "confirmatory"}
    for key in [("association", "phyloglm_confounders", "methanogen"),
                ("selection", "RELAX", "K")]:
        assert key in conf, (key, sorted(conf))
    assert any(f == "ancestral" for f, _, _ in conf), sorted(conf)
    base = r[r["test"].isin(["naive_glm", "naive_fisher"])]
    assert (base["p_adj"] == "").all() and (base["hypothesis"] == "exploratory").all()


@suite.test
def figures_follow_plotting_formats_and_provenance():
    root = world()
    R = root / "results"
    figs = [p for p in R.rglob("*") if p.suffix in (".png", ".pdf", ".svg")
            and ("plots" in p.parts or p.parent.name == "figures")]
    assert figs, "no figures"
    stems = {p.with_suffix("") for p in figs}
    bad = [str(s.relative_to(R)) for s in stems
           if not (s.with_suffix(".png").exists() and s.with_suffix(".pdf").exists())]
    assert not bad, bad[:10]
    svg = [str(p.relative_to(R)) for p in figs if p.suffix == ".svg"]
    assert not svg, svg[:10]
    pv = (R / "provenance" / "versions.txt").read_text()
    assert "defense-finder" in pv.lower() or "defensefinder" in pv.lower()


if __name__ == "__main__":
    try:
        suite.main("INTEGRATION TESTS PASSED")
    finally:
        print("timing (s):", STATE.get("timing"))
