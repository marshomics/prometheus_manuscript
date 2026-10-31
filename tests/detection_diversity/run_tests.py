#!/usr/bin/env python3
"""Leaf 1.1.2 (detection_diversity): synthetic tests with known answers + real-tool runs.

Covers check_defensefinder, parse_defensefinder, combine_proteomes, validate_prometheus,
architecture_typology, cluster_diversity, plot_prometheus_validation, plot_diversity,
plot_diversity_extra, plot_validation_tree.R and the SHELL BLOCKS of the rules in
workflow/rules/prometheus.smk / diversity.smk (formatted exactly as Snakemake does and run
under `bash -euo pipefail`), including their empty-input and failure behaviour.

Real tools: hmmscan/hmmpress, mafft, clipkit, iqtree2, mmseqs (PROM_TOOLS_BIN) and, when
installed, DefenseFinder (PROM_DF_BIN, else the sibling conda env
prom-detection_diversity-defensefinder). DefenseFinder's pinned models
(defense-finder-models 3.1.0 + CasFinder 3.1.0) come from PROM_DF_MODELS or are fetched
from codeload.github.com and installed with `macsydata install --target` (the offline
procedure documented in workflow/envs/defensefinder.yaml). Every tool uses <= 2 threads.
Prints DETECTION_DIVERSITY TESTS PASSED only when every test passes.
"""
from __future__ import annotations

import json
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path
from types import SimpleNamespace

sys.dont_write_bytecode = True   # never leave __pycache__ in shared script dirs
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, rscript, run, tmpdir, FIXTURES, read_tsv,  # noqa: E402
                  REPO, SCRIPTS_PY, SCRIPTS_R, PYTHON, TOOLS_BIN, env)

sys.path.insert(0, str(SCRIPTS_PY))
import utils  # noqa: E402
import validate_prometheus as vp  # noqa: E402

suite = Suite("detection_diversity")
TOKEN = "DETECTION_DIVERSITY TESTS PASSED"
THREADS = 2
OWNED_PY = ["check_defensefinder.py", "combine_proteomes.py", "parse_defensefinder.py",
            "validate_prometheus.py", "plot_prometheus_validation.py", "cluster_diversity.py",
            "architecture_typology.py", "plot_diversity.py", "plot_diversity_extra.py"]
OWNED_R = ["plot_validation_tree.R"]
RULES = {"prometheus": REPO / "workflow/rules/prometheus.smk",
         "diversity": REPO / "workflow/rules/diversity.smk"}


def config_spec_json() -> str:
    import yaml
    cfg = yaml.safe_load(open(REPO / "config/config.yaml"))
    return json.dumps(cfg["prometheus"]["domains"])


SPEC_JSON = config_spec_json()
SPEC = utils.load_domain_spec(SPEC_JSON)
_STATE = {}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def expect_fail(res, *needles):
    assert res.returncode != 0, "expected a non-zero exit"
    text = res.stdout + res.stderr
    for n in needles:
        assert n in text, f"expected '{n}' in output; got: {text[-1500:]}"


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def fixture_seqs():
    return dict(utils.read_fasta(str(FIXTURES / "proA_homologs.faa")))


def acc_of(name):
    return name.split("|")[1] if "|" in name else name


def by_acc(tag):
    return [v for k, v in fixture_seqs().items() if tag in k][0]


def python_shim() -> Path:
    """Directory with `python` -> PROM_PYTHON (rule shells call `python ...`)."""
    d = tmpdir("pyshim_")
    shim = d / "python"
    shim.write_text(f"#!/bin/bash\nexec {PYTHON} \"$@\"\n")
    shim.chmod(0o755)
    return d


def rule_shell(smk: Path, rule: str) -> str:
    text = smk.read_text()
    m = re.search(rf"^rule {re.escape(rule)}:\n(.*?)(?=^rule |\Z)", text, flags=re.S | re.M)
    assert m, f"rule {rule} not found in {smk}"
    sm = re.search(r'shell:\s*\n\s*r"""(.*?)"""', m.group(1), flags=re.S)
    assert sm, f"no shell block in rule {rule}"
    return sm.group(1)


def format_rule(smk_key, rule, ns) -> str:
    fmt = {k: (SimpleNamespace(**v) if isinstance(v, dict) else v) for k, v in ns.items()}
    return "set -euo pipefail\n" + rule_shell(RULES[smk_key], rule).format(**fmt)


def run_rule(smk_key, rule, base_env=None, timeout=1200, **ns):
    """Run a rule's shell block formatted like Snakemake (str.format on namespaces)."""
    e = dict(base_env or env())
    e["PATH"] = str(python_shim()) + os.pathsep + e["PATH"]
    return subprocess.run(["bash", "-c", format_rule(smk_key, rule, ns)], cwd=REPO, env=e,
                          text=True, capture_output=True, timeout=timeout)


def df_bin():
    p = os.environ.get("PROM_DF_BIN")
    if p and (Path(p) / "defense-finder").exists():
        return Path(p)
    if TOOLS_BIN:
        cand = Path(TOOLS_BIN).parent.parent / "prom-detection_diversity-defensefinder" / "bin"
        if (cand / "defense-finder").exists():
            return cand
    return None


def df_env():
    e = env()
    e["PATH"] = str(df_bin()) + os.pathsep + e["PATH"]
    e["HOME"] = str(tmpdir("dfhome_"))
    return e


def df_models():
    """Pinned models root (PROM_DF_MODELS, or download + macsydata install), else None.

    The installed models are cached under the temp dir (<tmp>/prometheus-tests-dfmodels-
    3.1.0_cas3.1.0) so reruns skip the ~67 MB download; the preflight re-verifies the
    versions on every run, so a corrupt cache fails loudly rather than being trusted."""
    if "models" in _STATE:
        return _STATE["models"]
    p = os.environ.get("PROM_DF_MODELS")
    if p and (Path(p) / "defense-finder-models").is_dir():
        _STATE["models"] = Path(p)
        return _STATE["models"]
    cache = Path(tempfile.gettempdir()) / "prometheus-tests-dfmodels-3.1.0_cas3.1.0"
    if all((cache / pkg / "metadata.yml").is_file() for pkg in ("defense-finder-models", "CasFinder")):
        _STATE["models"] = cache
        return cache
    root, dl = tmpdir("dfmodels_"), tmpdir("dfmodels_dl_")
    urls = {"defense-finder-models-3.1.0.tar.gz":
            "https://codeload.github.com/mdmparis/defense-finder-models/tar.gz/refs/tags/3.1.0",
            "CasFinder-3.1.0.tar.gz":
            "https://codeload.github.com/macsy-models/CasFinder/tar.gz/refs/tags/3.1.0"}
    try:
        for name, url in urls.items():
            with urllib.request.urlopen(url, timeout=180) as r, open(dl / name, "wb") as fh:
                shutil.copyfileobj(r, fh)
        for name in urls:
            subprocess.run([str(df_bin() / "macsydata"), "install", "--target", str(root), name],
                           cwd=dl, env=df_env(), check=True, capture_output=True, text=True,
                           timeout=600)
    except Exception as exc:  # noqa: BLE001 - offline: the real-DF test reports a skip
        print(f"  note  pinned DefenseFinder models unavailable ({exc})")
        _STATE["models"] = None
        return None
    try:                      # publish atomically; a concurrent run may have won the race
        os.replace(root, cache)
        root = cache
    except OSError:
        pass
    _STATE["models"] = root
    return root


# DefenseFinder v2/v3 table headers (defense_finder_posttreat + macsypy TsvSystemSerializer)
SYS_HDR = ["sys_id", "type", "subtype", "activity", "sys_beg", "sys_end", "protein_in_syst",
           "genes_count", "name_of_profiles_in_sys"]
GEN_HDR = ["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id", "sys_loci",
           "locus_num", "sys_wholeness", "sys_score", "sys_occ", "hit_gene_ref", "hit_status",
           "hit_seq_len", "hit_i_eval", "hit_score", "hit_profile_cov", "hit_seq_cov",
           "hit_begin_match", "hit_end_match", "counterpart", "used_in", "type", "subtype",
           "activity"]


def df_tables(gdir: Path, genome: str, systems: list):
    """systems: [(type, subtype, [(hit_id, gene_name, hit_pos, score)])]."""
    srows, grows = [], []
    for k, (typ, sub, genes) in enumerate(systems, 1):
        sid = f"{genome}_{sub}_{k}"
        ids = [g[0] for g in genes]
        srows.append([sid, typ, sub, "Defense", ids[0], ids[-1], ",".join(ids), str(len(ids)),
                      ",".join(g[1] for g in genes)])
        for hid, gname, pos, score in genes:
            row = dict.fromkeys(GEN_HDR, "")
            row.update(replicon=hid.rsplit("_", 1)[0], hit_id=hid, gene_name=gname,
                       hit_pos=str(pos), model_fqn=f"defense-finder-models/DefenseFinder/{typ}/{sub}",
                       sys_id=sid, hit_score=str(score), hit_i_eval="1e-100",
                       hit_profile_cov="0.95", hit_seq_cov="0.9", type=typ, subtype=sub,
                       activity="Defense", hit_status="mandatory")
            grows.append([row[c] for c in GEN_HDR])
    write(gdir / f"{genome}_defense_finder_systems.tsv",
          "\t".join(SYS_HDR) + "\n" + "".join("\t".join(r) + "\n" for r in srows))
    write(gdir / f"{genome}_defense_finder_genes.tsv",
          "\t".join(GEN_HDR) + "\n" + "".join("\t".join(r) + "\n" for r in grows))


def genes_tsv(gdir: Path, genome: str, genes: list):
    """genes: [(protein_id, contig_index, gene_index, partial, length_aa)] (Contract v2)."""
    cols = ["protein_id", "contig", "contig_index", "gene_index", "start", "end", "strand",
            "partial", "contig_length", "length_aa"]
    lines = ["\t".join(cols)] + ["\t".join([p, f"ctg{c}", str(c), str(g), str(g * 3000),
                                            str(g * 3000 + 3 * n + 2), "+", part, "900000", str(n)])
                                 for p, c, g, part, n in genes]
    write(gdir / genome / f"{genome}.genes.tsv", "\n".join(lines) + "\n")


def domline(target, acc, query, s, e):
    """One 22-column hmmscan --domtblout line (+ description)."""
    return (f"{target} {acc} 100 {query} - 1900 1e-10 100.0 0.0 1 1 1e-10 1e-10 100.0 0.0 "
            f"1 100 {s} {e} {s} {e} 0.90 desc\n")


def real_hmmscan():
    """Real HMMER 3.4 hmmscan --cut_ga --notextw of the fixture proteins vs the fixture
    Pfam subset; returns (domtbl, main output)."""
    if "hmm" not in _STATE:
        d = tmpdir("hmm_")
        shutil.copy(FIXTURES / "pfam_proA_subset.hmm", d / "pfam.hmm")
        run(["hmmpress", "-f", str(d / "pfam.hmm")])
        run(["hmmscan", "--cpu", str(THREADS), "--cut_ga", "--notextw", "-o", str(d / "out.txt"),
             "--domtblout", str(d / "x.domtbl"), str(d / "pfam.hmm"),
             str(FIXTURES / "proA_homologs.faa")])
        _STATE["hmm"] = (d / "x.domtbl", d / "out.txt")
    return _STATE["hmm"]


# --------------------------------------------------------------------------- #
# 0. static checks
# --------------------------------------------------------------------------- #
@suite.test
def compile_owned_python():
    for s in OWNED_PY:
        src = (SCRIPTS_PY / s).read_text()
        compile(src, str(SCRIPTS_PY / s), "exec")   # SyntaxError propagates; no .pyc written


@suite.test
def parse_owned_r():
    for s in OWNED_R:
        r = run([os.environ.get("PROM_RSCRIPT") or "Rscript", "-e",
                 f"invisible(parse(file='{SCRIPTS_R / s}')); cat('PARSED')"])
        assert "PARSED" in r.stdout


@suite.test
def rules_hygiene_and_no_hardcoded_domain_lists():
    r = run([PYTHON, str(REPO / "tests/check_rules.py"), str(RULES["prometheus"]),
             str(RULES["diversity"])])
    assert "RULES HYGIENE OK" in r.stdout
    for s in OWNED_PY:
        t = (SCRIPTS_PY / s).read_text()
        assert not re.search(r"\bCANON\s*=\s*\[", t), f"{s}: hard-coded CANON list"
        assert '"nuclease" in present' not in t, f"{s}: old nuclease-category logic"


# --------------------------------------------------------------------------- #
# 1. check_defensefinder (preflight, ingest copy, output verification)
# --------------------------------------------------------------------------- #
def fake_models(root: Path, df_vers="3.1.0", cas_vers="3.1.0", gene="Prometheus__ProA",
                presence="mandatory", grammar="2.0", with_profile=True):
    pkg = root / "defense-finder-models"
    write(pkg / "metadata.yml", f"---\nmaintainer:\n  name: x\nvers: {df_vers}\n")
    write(pkg / "definitions/DefenseFinder/Prometheus/Prometheus.xml",
          f'<model inter_gene_max_space="5" min_mandatory_genes_required="1" '
          f'min_genes_required="1" vers="{grammar}">\n<gene name="{gene}" '
          f'presence="{presence}"/>\n</model>\n')
    (pkg / "profiles").mkdir(exist_ok=True)
    if with_profile:
        write(pkg / f"profiles/{gene}.hmm", "HMMER3/f [3.4]\nNAME  Prometheus__ProA\n"
              "LENG  1668\nNSEQ  505\nGA    600 600;\nHMM   A C\n//\n")
    write(root / "CasFinder/metadata.yml", f"vers: {cas_vers}\n")
    write(root / "CasFinder/definitions/CAS_Class1.xml",
          f'<model vers="{grammar}"><gene name="cas3" presence="mandatory"/></model>\n')
    return root


def fake_df_bin(d: Path, version="3.0.0") -> Path:
    write(d / "defense-finder", "#!/bin/bash\nif [ \"$1\" = version ]; then echo "
          f"\"Using DefenseFinder version {version}\"; exit 0; fi\nexit 2\n").chmod(0o755)
    return d


def preflight(bin_dir, models, *extra):
    e = env()
    e["PATH"] = str(bin_dir) + os.pathsep + e["PATH"]
    out = tmpdir("pre_") / "ok.tsv"
    p = subprocess.run([PYTHON, str(SCRIPTS_PY / "check_defensefinder.py"), "--run", "True",
                        "--models-dir", str(models), "--out", str(out), *extra],
                       env=e, text=True, capture_output=True)
    return p, out


@suite.test
def preflight_version_and_model_checks():
    d = tmpdir("pf_")
    ok_models, b = fake_models(d / "m_ok"), fake_df_bin(d / "bin")
    p, out = preflight(b, ok_models, "--software-version", "3.0.0", "--models-pin", "3.1.0",
                       "--casfinder-pin", "3.1.0")
    assert p.returncode == 0, p.stderr
    kv = dict(l.split("\t", 1) for l in out.read_text().splitlines()[1:])
    assert (kv["defense_finder_version"], kv["defense_finder_models_version"],
            kv["casfinder_version"]) == ("3.0.0", "3.1.0", "3.1.0")
    assert kv["prometheus_genes"] == "Prometheus__ProA:mandatory"
    assert "LENG=1668" in kv["prometheus_profile_header"] and "GA=600 600;" in kv["prometheus_profile_header"]
    expect_fail(preflight(b, ok_models, "--software-version", "2.0.1")[0], "!= pinned 2.0.1")
    expect_fail(preflight(fake_df_bin(d / "old", "1.3.0"), ok_models)[0], "< required 2.0.0")
    expect_fail(preflight(b, ok_models, "--models-pin", "3.0.0")[0], "!= pinned 3.0.0")
    expect_fail(preflight(b, fake_models(d / "m_v1", df_vers="1.2.4"))[0], "< 2.0.0")
    expect_fail(preflight(b, fake_models(d / "m_cas", cas_vers="3.1.1"), "--casfinder-pin",
                          "3.1.0")[0], "CasFinder 3.1.1 != pinned 3.1.0")
    expect_fail(preflight(b, fake_models(d / "m_gen", gene="Prometheus__ProX"))[0],
                "does not declare gene 'Prometheus__ProA'")
    expect_fail(preflight(b, fake_models(d / "m_acc", presence="accessory"))[0], "expected 'mandatory'")
    expect_fail(preflight(b, fake_models(d / "m_prof", with_profile=False))[0], "is missing")
    nocas = fake_models(d / "m_nocas")
    shutil.rmtree(nocas / "CasFinder")
    expect_fail(preflight(b, nocas)[0], "'CasFinder' is not installed")
    expect_fail(preflight(tmpdir("nobin_"), ok_models)[0], "cannot determine the defense-finder version")


@suite.test
def preflight_ingest_mode():
    d = tmpdir("ing_")
    expect_fail(py("check_defensefinder.py", "--run", "False", "--provided-dir", d / "none",
                   "--out", d / "ok.tsv", check=False), "does not exist")
    df_tables(d / "prov" / "G1", "G1", [])
    py("check_defensefinder.py", "--run", "False", "--provided-dir", d / "prov", "--out", d / "ok.tsv")
    assert "UNVERIFIED" in (d / "ok.tsv").read_text()


@suite.test
def check_outputs_and_ingest_copy():
    d = tmpdir("chk_")
    faa = write(d / "G1.faa", ">G1_1_1\nMKV\n>G1_1_2\nMKV\n>G1_1_3\nMKV\n")
    df_tables(d / "ok", "G1", [("Septu", "Septu", [("G1_1_1", "Septu__PtuA", 1, 100),
                                                   ("G1_1_2", "Septu__PtuB", 2, 90)])])
    chk = lambda sub, *x, check=True: py(
        "check_defensefinder.py", "--check-outputs", *x, "--systems", d / sub / "G1_defense_finder_systems.tsv",
        "--genes", d / sub / "G1_defense_finder_genes.tsv", "--faa", faa, "--genome", "G1", check=check)
    chk("ok", "--expect-gembase")
    df_tables(d / "bad", "G1", [("Septu", "Septu", [("G1_9_9", "Septu__PtuA", 1, 100)])])
    expect_fail(chk("bad", check=False), "absent from", "G1_9_9")
    for sub, old, new, flag, msg in (("v1", "DefenseFinder/Septu", "DefenseFinder_3/Septu", (),
                                      "defense-finder-models v1"),
                                     ("or", "\nG1_1\t", "\nG1\t", ("--expect-gembase",), "--db-type gembase")):
        write(d / sub / "G1_defense_finder_genes.tsv",
              (d / "ok/G1_defense_finder_genes.tsv").read_text().replace(old, new))
        shutil.copy(d / "ok/G1_defense_finder_systems.tsv", d / sub / "G1_defense_finder_systems.tsv")
        expect_fail(chk(sub, *flag, check=False), msg)
    write(d / "mc/G1_defense_finder_systems.tsv", "sys_id\ttype\n")
    shutil.copy(d / "ok/G1_defense_finder_genes.tsv", d / "mc/G1_defense_finder_genes.tsv")
    expect_fail(chk("mc", check=False), "lacks DefenseFinder v2 columns")
    # ingest copy: exact name, a single glob match, ambiguous, missing
    prov = d / "prov"
    df_tables(prov / "G1", "G1", [])
    cp = lambda g, o, check=True: py("check_defensefinder.py", "--ingest-copy", "--provided-dir", prov,
                                     "--genome", g, "--systems", d / o / "s.tsv", "--genes",
                                     d / o / "g.tsv", check=check)
    cp("G1", "o1")
    assert (d / "o1/s.tsv").read_text().startswith("sys_id\ttype")
    (prov / "G2").mkdir()
    for k in ("systems", "genes"):
        shutil.copy(prov / f"G1/G1_defense_finder_{k}.tsv", prov / f"G2/G2.prt_defense_finder_{k}.tsv")
    cp("G2", "o2")
    shutil.copy(prov / "G1/G1_defense_finder_genes.tsv", prov / "G2/other_defense_finder_genes.tsv")
    expect_fail(cp("G2", "o3", check=False), "found 2")
    expect_fail(cp("G9", "o4", check=False), "found 0")


@suite.test
def defensefinder_rule_ingest_mode():
    """Per-genome rule in ingest mode: copy + verification; foreign protein ids -> FAIL."""
    d = tmpdir("dfrule_")
    faa = write(d / "genes/G1/G1.faa", ">G1_1_1\nMKV\n>G1_1_2\nMKV\n")
    df_tables(d / "prov/G1", "G1", [("Septu", "Septu", [("G1_1_1", "Septu__PtuA", 1, 100)])])
    out = d / "out/G1"
    ns = dict(input={"faa": faa, "models_ok": d / "ok"}, wildcards={"genome": "G1"},
              output={"systems": out / "G1_defense_finder_systems.tsv",
                      "genes": out / "G1_defense_finder_genes.tsv"},
              params={"run": "False", "provided": d / "prov", "outdir": out, "models_opt": "",
                      "sys_header": "\t".join(SYS_HDR), "genes_header": "\t".join(GEN_HDR)},
              threads=1, log=d / "log.txt")
    p = run_rule("prometheus", "defensefinder", **ns)
    assert p.returncode == 0, (d / "log.txt").read_text()
    assert read_tsv(out / "G1_defense_finder_genes.tsv")["hit_id"].tolist() == ["G1_1_1"]
    write(faa, ">G1_7_1\nMKV\n")          # the proteome no longer matches the provided output
    p = run_rule("prometheus", "defensefinder", **ns)
    assert p.returncode != 0 and "absent from" in (d / "log.txt").read_text()


# --------------------------------------------------------------------------- #
# 2. combine_proteomes
# --------------------------------------------------------------------------- #
@suite.test
def combine_proteomes_proteome_mode():
    d = tmpdir("cp_")
    write(d / "samples.tsv", "genome_id\tx\nGA\t1\nGB\t2\n")
    write(d / "genes/GA/GA.faa", ">GA_1_1\nMKV*\n>GA_1_2\nMLL\n")
    write(d / "genes/GB/GB.faa", "")
    py("combine_proteomes.py", "--samples", d / "samples.tsv", "--genes-dir", d / "genes",
       "--out-faa", d / "all.faa", "--out-map", d / "map.tsv")
    assert utils.fasta_ids(str(d / "all.faa")) == ["GA__GA_1_1", "GA__GA_1_2"]
    m = read_tsv(d / "map.tsv")
    assert m.columns.tolist() == ["protein_id", "genome_id", "orig_id"]
    assert m["orig_id"].tolist() == ["GA_1_1", "GA_1_2"] and set(m["genome_id"]) == {"GA"}
    # CDS mode used by evolution.smk combined_cds (--suffix fna) stays compatible
    write(d / "genes/GA/GA.fna", ">GA_1_1\nATG\n")
    write(d / "genes/GB/GB.fna", ">GB_1_1\nATG\n")
    py("combine_proteomes.py", "--samples", d / "samples.tsv", "--genes-dir", d / "genes",
       "--suffix", "fna", "--out-faa", d / "all.ffn", "--out-map", d / "cds.tsv")
    assert utils.fasta_ids(str(d / "all.ffn")) == ["GA__GA_1_1", "GB__GB_1_1"]
    write(d / "s2.tsv", "genome_id\nGA\nGC\n")
    expect_fail(py("combine_proteomes.py", "--samples", d / "s2.tsv", "--genes-dir", d / "genes",
                   "--out-faa", d / "x.faa", "--out-map", d / "x.tsv", check=False), "missing")
    expect_fail(py("combine_proteomes.py", "--samples", d / "nope.tsv", "--genes-dir", d / "genes",
                   "--out-faa", d / "x.faa", "--out-map", d / "x.tsv", check=False), "cannot read --samples")


@suite.test
def combine_proteomes_reference_mode():
    d = tmpdir("cpref_")
    write(d / "cand.faa", ">GA__GA_1_5\nMKVL\n")
    write(d / "prom.faa", ">sp|Q9HXX1|PROA_PSEAE some description\nMKV\n>lcl|x:1-2\nMKA\n")
    write(d / "para.faa", ">tr|P38859|DNA2_YEAST\nMKL\n")
    ref = lambda prom, para: ["--ref", f"Prometheus_ref={prom}", "--ref", f"Paralog_ref={para}"]
    py("combine_proteomes.py", "--candidates", d / "cand.faa", *ref(d / "prom.faa", d / "para.faa"),
       "--out-faa", d / "all.faa", "--out-labels", d / "labels.tsv", "--out-map", d / "map.tsv")
    assert utils.fasta_ids(str(d / "all.faa")) == ["GA__GA_1_5", "sp_Q9HXX1_PROA_PSEAE",
                                                   "lcl_x_1-2", "tr_P38859_DNA2_YEAST"]
    assert (d / "labels.tsv").read_text().splitlines() == [
        "sp_Q9HXX1_PROA_PSEAE\tPrometheus_ref", "lcl_x_1-2\tPrometheus_ref",
        "tr_P38859_DNA2_YEAST\tParalog_ref"]
    write(d / "dup.faa", ">sp|A|B\nMK\n>sp:A:B\nMK\n")      # both sanitise to sp_A_B
    expect_fail(py("combine_proteomes.py", "--candidates", d / "cand.faa", *ref(d / "dup.faa", d / "para.faa"),
                   "--out-faa", d / "x.faa", "--out-labels", d / "x.tsv", check=False), "collides")
    write(d / "empty.faa", "")
    expect_fail(py("combine_proteomes.py", "--candidates", d / "cand.faa", *ref(d / "prom.faa", d / "empty.faa"),
                   "--out-faa", d / "x.faa", "--out-labels", d / "x.tsv", check=False), "are empty")


# --------------------------------------------------------------------------- #
# 3. parse_defensefinder
# --------------------------------------------------------------------------- #
def parse_case(d: Path):
    write(d / "samples.tsv", "genome_id\nG1\nG2\n")
    write(d / "all.faa", ">G1__G1_1_1\nMAAA\n>G1__G1_1_2\nMCCC\n>G1__G1_1_3\nMDDD\n>G2__G2_1_1\nMEEE\n")
    write(d / "map.tsv", "protein_id\tgenome_id\torig_id\nG1__G1_1_1\tG1\tG1_1_1\n"
          "G1__G1_1_2\tG1\tG1_1_2\nG1__G1_1_3\tG1\tG1_1_3\nG2__G2_1_1\tG2\tG2_1_1\n")
    return ["--df-dir", d / "df", "--proteome", d / "all.faa", "--pmap", d / "map.tsv",
            "--samples", d / "samples.tsv", "--out-candidates", d / "c.faa",
            "--out-sources", d / "s.tsv", "--out-raw", d / "r.tsv"]


@suite.test
def parse_defensefinder_exact_matching():
    d = tmpdir("pdf_")
    base = parse_case(d)
    df_tables(d / "df/G1", "G1", [
        ("Prometheus", "Prometheus", [("G1_1_1", "Prometheus__ProA", 1, 1024.5),
                                      ("G1_1_2", "Prometheus__Accessory", 2, 50)]),
        ("NotPrometheusLike", "Prometheus_relative", [("G1_1_3", "X__Y", 3, 99)])])
    df_tables(d / "df/G2", "G2", [("Septu", "Septu", [("G2_1_1", "Septu__PtuA", 1, 77)])])
    df_tables(d / "df/STALE", "STALE", [("Prometheus", "Prometheus", [("STALE_1_1", "Prometheus__ProA", 1, 5)])])
    py("parse_defensefinder.py", *base, "--system", "Prometheus")
    assert utils.fasta_ids(str(d / "c.faa")) == ["G1__G1_1_1"]   # ProA gene only, no substring hits
    src = read_tsv(d / "s.tsv")
    assert src.columns[:5].tolist() == ["protein_id", "genome_id", "detect_defensefinder",
                                        "df_score", "df_subtype"]
    r0 = src.iloc[0]
    assert (r0.protein_id, r0.genome_id, r0.detect_defensefinder, float(r0.df_score),
            r0.df_subtype, r0.df_gene_name) == ("G1__G1_1_1", "G1", "1", 1024.5, "Prometheus",
                                                "Prometheus__ProA")
    assert read_tsv(d / "r.tsv")["protein_id"].tolist() == ["G1__G1_1_1"]
    py("parse_defensefinder.py", *base, "--system", "prometheus", "--gene", "")
    assert utils.fasta_ids(str(d / "c.faa")) == ["G1__G1_1_1", "G1__G1_1_2"]


@suite.test
def parse_defensefinder_failures_and_empty():
    d = tmpdir("pdf2_")
    base = parse_case(d) + ["--system", "Prometheus"]
    df_tables(d / "df/G1", "G1", [])
    df_tables(d / "df/G2", "G2", [])
    expect_fail(py("parse_defensefinder.py", *base, check=False), "ZERO systems")
    py("parse_defensefinder.py", *base, "--strict", "False")
    assert (d / "c.faa").read_text() == "" and read_tsv(d / "s.tsv").empty
    df_tables(d / "df/G1", "G1", [("Prometheus", "Prometheus", [("G1_4_4", "Prometheus__ProA", 1, 900)])])
    expect_fail(py("parse_defensefinder.py", *base, check=False), "absent from the pipeline proteome")
    write(d / "samples.tsv", "genome_id\nG1\nG2\nG9\n")
    df_tables(d / "df/G1", "G1", [("Septu", "Septu", [("G1_1_1", "Septu__PtuA", 1, 9)])])
    expect_fail(py("parse_defensefinder.py", *base, check=False), "missing DefenseFinder systems table")


# --------------------------------------------------------------------------- #
# 4. validate_prometheus
# --------------------------------------------------------------------------- #
@suite.test
def real_hmmscan_matches_fixture_and_motifs():
    domtbl, out = real_hmmscan()
    key = lambda p: sorted((f[3], f[0], f[17], f[18]) for f in
                           (l.split() for l in open(p) if not l.startswith("#")))
    assert key(domtbl) == key(FIXTURES / "proA_homologs.pfam_subset.domtbl"), "domtbl drift"
    hits, alns = utils.parse_domtbl(str(domtbl)), vp.parse_hmmscan_alignments(str(out))
    assert sum(len(v) for v in alns.values()) == sum(len(v) for v in hits.values()) == 62
    sites = vp.acidic_sites_by_model(alns)
    assert sites["MTES_1575"] == [28, 45]
    seqs = {k: vp.strip_stop(v) for k, v in fixture_seqs().items()}
    got = {}
    for q in hits:
        a = alns.get(q, [])
        acc = {h["domain"]: h["accession"] for h in hits[q]}
        for x in a:
            x["accession"] = acc.get(x["model"], "")
        rec = vp.architecture_record(hits[q], SPEC)
        got[acc_of(q)] = (rec, vp.catalytic_assessment(rec["present"], a, SPEC, "DNA2", "REase",
                                                       True, True, sites, seqs[q], hits[q]))
    assert got["A0A060RCK4"][0]["architecture"] == "DUF4011-DNA2-REase-DUF3320"
    h = got["A0A1H8VVA9"][0]
    assert (h["architecture"], h["missing_required"], h["canonical_order"], h["other_domains"]) == \
        ("DUF4011-DNA2-REase", [], True, ["ResIII", "Helicase_C"])
    assert got["A0AAV5N545"][0]["other_domains"] == ["Pkinase"]
    assert got["A0A060RCK4"][1][:3] == ("intact", "intact", "intact")
    assert "WalkerA:GPPGTGKS@305(AAA_11)" in got["A0A060RCK4"][1][3]
    assert got["A0A011V4U2"][1][:3] == ("partial", "intact", "disrupted")
    assert "D45:C@1642" in got["A0A011V4U2"][1][3]
    assert got["A0A0A6X990"][1][0] == "partial" and "D28:G@1461" in got["A0A0A6X990"][1][3]
    cons = got["Prometheus__ProA-consensus"][1]
    assert cons[0] == "intact" and "MotifII:DE@993(window)" in cons[3]
    assert sum(1 for _, st in got.values() if st[0] == "intact") == 9


@suite.test
def catalytic_status_edge_cases():
    assert vp.catalytic_assessment(set(), [], SPEC, "DNA2", "REase", True, False)[0] == "unknown"
    assert vp.catalytic_assessment({"DNA2", "REase"}, [], SPEC, "DNA2", "REase", False, True) == \
        ("intact", "present", "present", "domain_proxy")
    assert vp.catalytic_assessment({"DNA2"}, [], SPEC, "DNA2", "REase", False, True)[0] == "partial"
    aln = {"model": "AAA_11", "accession": "PF13086", "hmm_from": 20, "hmm_to": 33,
           "cons": "sliqGPPGtGKTkt", "tgt": "SLIQGPPGTGASKT", "pos": list(range(100, 114))}
    status, ev = vp._assess_helicase([aln])
    assert status == "disrupted" and "WalkerA:GPPGTGAS" in ev       # catalytic K -> A
    assert vp.catalytic_assessment({"DNA2"}, [aln], SPEC, "DNA2", "REase", True, True)[0] == "degenerate"


def tree_and_labels(d: Path, newick: str, labels: dict):
    write(d / "tree.nwk", newick)
    write(d / "labels.tsv", "".join(f"{k}\t{v}\n" for k, v in labels.items()))
    return d / "tree.nwk", vp.load_ref_labels(str(d / "labels.tsv"))


REFS5 = {"P1": "Prometheus_ref", "P2": "Prometheus_ref", "P3": "Prometheus_ref",
         "Q1": "Paralog_ref", "Q2": "Paralog_ref"}


@suite.test
def placement_outgroup_rooting_and_support():
    nwk = ("((Q1:0.1,(Q2:0.1,cB:0.05)99/100:0.1)100/100:0.5,"
           "(cS:0.3,((P1:0.1,cA:0.05)98/100:0.1,(P2:0.1,P3:0.1)95/99:0.1)100/100:0.2)90/97:0.3,"
           "cX:0.4);")
    t, lab = tree_and_labels(tmpdir("plc_"), nwk, REFS5)
    det = {}
    a = vp.assign_by_tree(str(t), lab, 80, 95, det)
    assert (a["cA"], det["cA"]["clade"], det["cA"]["support"]) == ("Prometheus", "crown", "98/100")
    assert (a["cB"], det["cB"]["clade"]) == ("Paralog", "outgroup")
    assert (a["cS"], det["cS"]["clade"], det["cS"]["support"]) == ("Prometheus", "stem", "90/97")
    assert (a["cX"], det["cX"]["clade"], det["cX"]["support"]) == ("Prometheus", "stem", "100/100")
    weak = re.sub(r"\d+/\d+", "50/60", nwk)                     # same topology, weak supports
    t2, _ = tree_and_labels(tmpdir("plc2_"), weak, REFS5)
    det2 = {}
    assert set(vp.assign_by_tree(str(t2), lab, 80, 95, det2).values()) == {"Unplaced"}
    assert det2["cA"]["note"].startswith("unsupported_ingroup")
    assert det2["cB"]["note"].startswith("unsupported_outgroup")
    t3, lab3 = tree_and_labels(tmpdir("plc3_"), "((P1,Q1)100/100,(P2,(Q2,cA)100/100)100/100,P3);", REFS5)
    det3 = {}
    assert vp.assign_by_tree(str(t3), lab3, 80, 95, det3) == {"cA": "Unplaced"}
    assert det3["cA"]["note"] == "references_not_separable"
    # IQ-TREE-renamed tips ('sp|A|PROA:1' -> 'sp|A|PROA_1') still match sanitised labels
    t4, lab4 = tree_and_labels(tmpdir("plc4_"), "((sp|A|PROA_1:0.1,cA:0.1)99/99:0.1,"
                               "(tr|B|DNA2:0.1,tr|C|NURA:0.1)99/99:0.1,sp|D|PROA_2:0.1);",
                               {"sp|A|PROA:1": "Prometheus_ref", "sp|D|PROA_2": "Prometheus_ref",
                                "tr|B|DNA2": "Paralog_ref", "tr|C|NURA": "Paralog_ref"})
    assert vp.assign_by_tree(str(t4), lab4, 80, 95) == {"cA": "Prometheus"}
    d = tmpdir("plc5_")
    assert vp.assign_by_tree(str(write(d / "ph.nwk", "();\n")), lab) == {}
    try:
        vp.assign_by_tree(str(write(d / "bad.nwk", "((a,b),(c;\n")), lab)
        raise AssertionError("malformed tree accepted")
    except vp.InputError:
        pass


@suite.test
def tier_rules():
    T = vp.assign_tier
    full = dict(detect_df=1, missing_required=[], canonical=True, tree_assignment="Prometheus",
                fragmented=False, catalytic="intact", defense_ctx=False, length_aa=2221,
                length_in_bounds=True, min_len=800, max_len=3000)
    assert T(**full) == ("A", "")
    # DUF3320 is optional (required: false): a canonical core without it can be tier A
    rec = vp.architecture_record([domtbl_hit("DUF4011", "PF13195.12", 50, 200),
                                  domtbl_hit("AAA_11", "PF13086.13", 300, 400),
                                  domtbl_hit("MTES_1575", "PF18741.8", 1300, 1390)], SPEC)
    assert (rec["missing_required"], rec["canonical_order"], rec["architecture"]) == \
        ([], True, "DUF4011-DNA2-REase")
    assert T(**{**full, "missing_required": rec["missing_required"]})[0] == "A"
    assert T(**{**full, "missing_required": ["DUF4011"]}) == ("B", "prometheus_placement;missing_required:DUF4011")
    assert T(**{**full, "canonical": False})[0] == "B"
    assert T(**{**full, "catalytic": "partial"}) == ("B", "prometheus_placement;catalytic:partial")
    assert T(**{**full, "fragmented": True})[0] == "B"
    assert T(**{**full, "length_aa": 3200, "length_in_bounds": False})[0] == "B"
    assert T(**{**full, "tree_assignment": "Paralog"}) == ("Rejected", "placed_with_cellular_paralog")
    assert T(**{**full, "tree_assignment": "Unplaced", "defense_ctx": True})[0] == "B"
    assert T(**{**full, "tree_assignment": "Unplaced"}, placement_note="no_tree") == ("C", "placement:no_tree")
    assert T(**{**full, "length_aa": 300, "length_in_bounds": False})[0] == "Rejected"
    assert T(**{**full, "length_aa": 300, "length_in_bounds": False, "fragmented": True})[0] == "B"
    assert T(**{**full, "length_aa": 7000, "length_in_bounds": False})[0] == "Rejected"


def domtbl_hit(name, acc, s, e):
    return {"domain": name, "accession": acc, "start": s, "end": e}


def build_validate_case(d: Path):
    """Five candidates in four genomes with known answers (validate_end_to_end)."""
    seqs = fixture_seqs()
    pick = {"A": "A0A060RCK4", "N": "A0A011V4U2", "H": "A0A1H8VVA9", "K": "A0AAV5N545",
            "S": "A0A0A7PKR5"}
    src = {k: [n for n in seqs if a in n][0] for k, a in pick.items()}
    cands = {"G1__G1_1_5": "A", "G2__G2_3_10": "N", "G3__G3_1_1": "H", "G3__G3_1_9": "K",
             "G4__G4_2_2": "S"}
    seq_of = lambda k: seqs[src[k]][:350] if k == "S" else seqs[src[k]]
    write(d / "cand.faa", "".join(f">{p}\n{seq_of(k)}*\n" for p, k in cands.items()))
    write(d / "sources.tsv", "protein_id\tgenome_id\tdetect_defensefinder\tdf_score\tdf_subtype\n"
          + "".join(f"{p}\t{utils.genome_of(p)}\t1\t1000\tPrometheus\n" for p in cands))
    write(d / "raw.tsv", "protein_id\tgenome_id\tmethod\tscore\tsubtype\n")
    domtbl, out = real_hmmscan()
    lines = []
    for p, k in cands.items():
        for ln in open(domtbl):
            f = ln.split()
            if ln.startswith("#") or f[3] != src[k]:
                continue
            if k == "S" and f[0] != "DUF4011":       # the 350-aa fragment keeps DUF4011 only
                continue
            lines.append(ln.replace(src[k], p, 1))
    write(d / "x.domtbl", "".join(lines))
    blocks = re.split(r"(?m)^(?=Query:)", open(out).read())
    keep = [b.replace(src[k], p) for p, k in cands.items() if k != "S"
            for b in blocks if b.startswith(f"Query:       {src[k]} ")]
    write(d / "x.hmmscan.txt", "".join(keep))
    # G1 supported with P1; G2 only weakly with Prometheus refs; G3_1_1 nested in paralogs;
    # G3_1_9 and G4 are absent from the tree
    write(d / "tree.nwk", "((Q1:0.1,(Q2:0.1,G3__G3_1_1:0.05)99/100:0.1)50/50:0.5,"
          "((P1:0.1,G1__G1_1_5:0.05)98/100:0.1,(P2:0.1,G2__G2_3_10:0.1)40/55:0.1)60/70:0.2,P3:0.3);")
    write(d / "labels.tsv", "".join(f"{k}\t{v}\n" for k, v in REFS5.items()))
    write(d / "meta.tsv", "genome_id\tcompleteness\n" + "".join(f"G{i}\t95\n" for i in range(1, 7)))
    write(d / "reps.txt", "G1\nG2\nG3\nG4\nG5\n")            # G6 is not a representative
    df_tables(d / "df/G1", "G1", [("Prometheus", "Prometheus", [("G1_1_5", "Prometheus__ProA", 5, 1000)]),
                                  ("Septu", "Septu", [("G1_1_9", "Septu__PtuA", 9, 100),
                                                      ("G1_1_10", "Septu__PtuB", 10, 90)])])
    df_tables(d / "df/G2", "G2", [("Prometheus", "Prometheus", [("G2_3_10", "Prometheus__ProA", 10, 1000)]),
                                  ("RM", "RM_Type_I", [("G2_3_21", "RM__X", 21, 100)]),
                                  ("Gabija", "Gabija", [("G2_4_10", "Gabija__GajA", 30, 100)])])
    df_tables(d / "df/G3", "G3", [("Prometheus", "Prometheus", [("G3_1_1", "Prometheus__ProA", 1, 1000)]),
                                  ("Prometheus", "Prometheus", [("G3_1_9", "Prometheus__ProA", 9, 900)]),
                                  ("CBASS", "CBASS_I", [("G3_1_12", "CBASS__Cyclase", 12, 80)])])
    df_tables(d / "df/G4", "G4", [("Prometheus", "Prometheus", [("G4_2_2", "Prometheus__ProA", 2, 700)])])
    L = lambda k: len(seq_of(k))
    genes_tsv(d / "genes", "G1", [("G1_1_5", 1, 5, "00", L("A")), ("G1_1_9", 1, 9, "00", 300),
                                  ("G1_1_10", 1, 10, "00", 300)])
    genes_tsv(d / "genes", "G2", [("G2_3_10", 3, 10, "00", L("N")), ("G2_3_21", 3, 21, "00", 300),
                                  ("G2_4_10", 4, 10, "00", 300)])
    genes_tsv(d / "genes", "G3", [("G3_1_1", 1, 1, "00", L("H")), ("G3_1_9", 1, 9, "00", L("K")),
                                  ("G3_1_12", 1, 12, "00", 300)])
    genes_tsv(d / "genes", "G4", [("G4_2_2", 2, 2, "10", 350)])


def validate_args(d: Path, **over):
    a = {"--sources": d / "sources.tsv", "--raw": d / "raw.tsv", "--domtbl": d / "x.domtbl",
         "--tree": d / "tree.nwk", "--ref-labels": d / "labels.tsv", "--candidates": d / "cand.faa",
         "--meta": d / "meta.tsv", "--reps": d / "reps.txt", "--df-dir": d / "df",
         "--context-window": "10", "--domains": SPEC_JSON, "--min-len": "800", "--max-len": "3000",
         "--out-catalog": d / "cat.tsv", "--out-proteins": d / "prot.faa",
         "--out-sensitivity": d / "sens.tsv", "--genes-dir": d / "genes",
         "--hmmscan-out": d / "x.hmmscan.txt"}
    a.update(over)
    return [x for k, v in a.items() if v is not None for x in (k, v)]


@suite.test
def validate_end_to_end_known_answers():
    d = tmpdir("val_")
    build_validate_case(d)
    py("validate_prometheus.py", *validate_args(d))
    cat = read_tsv(d / "cat.tsv").set_index("locus_id", drop=False)
    assert cat.columns[:13].tolist() == vp.BASE_COLUMNS            # SCHEMA columns, in order
    for c in ("has_DUF4011", "has_DNA2", "has_REase", "has_DUF3320", "other_domains",
              "canonical_order", "placement_support", "context_method", "partial"):
        assert c in cat.columns, c
    g1 = cat.loc["G1__G1_1_5"]
    assert (g1.confidence_tier, g1.tree_assignment, g1.catalytic_status, g1.defense_context,
            g1.fragmented, g1.architecture, g1.length_aa) == \
        ("A", "Prometheus", "intact", "True", "False", "DUF4011-DNA2-REase-DUF3320", "1828")
    assert (g1.context_n_systems, g1.context_systems, g1.context_method, g1.placement_support) == \
        ("1", "Septu:G1_Septu_2", "neighbourhood_+/-10_genes_tsv", "98/100")
    g2 = cat.loc["G2__G2_3_10"]       # weakly placed; RM 11 genes away; Gabija on contig 4
    assert (g2.tree_assignment, g2.defense_context, g2.catalytic_status, g2.confidence_tier) == \
        ("Unplaced", "False", "partial", "C")
    assert "placement:unsupported_ingroup" in g2.reason
    g3 = cat.loc["G3__G3_1_1"]
    assert (g3.confidence_tier, g3.tree_assignment, g3.other_domains, g3.has_DUF3320) == \
        ("Rejected", "Paralog", "ResIII,Helicase_C", "False")
    g3b = cat.loc["G3__G3_1_9"]       # absent from the tree; CBASS 3 genes away
    assert (g3b.tree_assignment, g3b.placement_note, g3b.defense_context, g3b.confidence_tier) == \
        ("Unplaced", "not_in_tree", "True", "B")
    assert g3b.context_systems == "CBASS:G3_CBASS_I_3"   # the other Prometheus copy is not context
    assert "core_architecture+defense_context" in g3b.reason
    g4 = cat.loc["G4__G4_2_2"]        # 350-aa ORF running off a contig end: fragmented, kept
    assert (g4.fragmented, g4.partial, g4.fragment_evidence, g4.confidence_tier) == \
        ("True", "10", "prodigal_partial=10", "C")
    assert "missing_required:DNA2,REase" in g4.reason
    assert set(utils.fasta_ids(str(d / "prot.faa"))) == set(cat.index) - {"G3__G3_1_1"}
    sens = read_tsv(d / "sens.tsv").set_index("genome_id")
    assert sens.index.tolist() == ["G1", "G2", "G3", "G4", "G5"]        # meta restricted to reps
    assert sens["prometheus_present"].tolist() == ["True", "True", "True", "True", "False"]
    assert sens.loc["G5", "detect_defensefinder"] == "0"
    py("validate_prometheus.py", *validate_args(d, **{"--context-window": "3",
                                                      "--out-catalog": d / "c3.tsv"}))
    assert read_tsv(d / "c3.tsv").set_index("locus_id").loc["G1__G1_1_5", "defense_context"] == "False"
    py("validate_prometheus.py", *validate_args(d, **{"--genes-dir": None, "--out-catalog": d / "cng.tsv"}))
    ng = read_tsv(d / "cng.tsv").set_index("locus_id")
    assert ng.loc["G1__G1_1_5", "context_method"] == "neighbourhood_+/-10_gene_id"
    assert ng.loc["G4__G4_2_2", "fragment_evidence"] == "length_proxy(length<min_len)"
    py("validate_prometheus.py", *validate_args(d, **{"--hmmscan-out": None, "--out-catalog": d / "cnp.tsv"}))
    npx = read_tsv(d / "cnp.tsv").set_index("locus_id")
    assert (npx.loc["G2__G2_3_10", "catalytic_status"], npx.loc["G2__G2_3_10", "catalytic_evidence"]) == \
        ("intact", "domain_proxy")
    _STATE["validate_case"] = d


@suite.test
def validate_empty_and_failures():
    d = tmpdir("val2_")
    build_validate_case(d)
    write(d / "esrc.tsv", "protein_id\tgenome_id\tdetect_defensefinder\tdf_score\tdf_subtype\n")
    py("validate_prometheus.py", *validate_args(d, **{"--sources": d / "esrc.tsv"}))
    cat = read_tsv(d / "cat.tsv")
    assert cat.empty and cat.columns.tolist() == vp.catalog_columns(SPEC)
    assert read_tsv(d / "sens.tsv")["prometheus_present"].tolist() == ["False"] * 5
    assert (d / "prot.faa").read_text() == ""
    write(d / "cmiss.faa", f">G1__G1_1_5\n{by_acc('A0A060RCK4')}\n")
    expect_fail(py("validate_prometheus.py", *validate_args(d, **{"--candidates": d / "cmiss.faa"}),
                   check=False), "have no sequence")
    write(d / "bad.nwk", "((P1,G1__G1_1_5),(Q1;")
    expect_fail(py("validate_prometheus.py", *validate_args(d, **{"--tree": d / "bad.nwk"}), check=False),
                "could not parse validation tree")
    expect_fail(py("validate_prometheus.py", *validate_args(d, **{"--domains": "{not json"}), check=False),
                "invalid --domains JSON")


# --------------------------------------------------------------------------- #
# 5. rules with real HMMER / MAFFT / ClipKIT / IQ-TREE; R figure
# --------------------------------------------------------------------------- #
@suite.test
def domain_scan_rule_fails_on_unpressed_pfam():
    d = tmpdir("dsr_")
    shutil.copy(FIXTURES / "pfam_proA_subset.hmm", d / "Pfam-A.hmm")
    write(d / "cand.faa", f">G1__G1_1_5\n{by_acc('A0A060RCK4')}*\n")
    ns = dict(input={"candidates": d / "cand.faa", "pfam": d / "Pfam-A.hmm"},
              output={"domtbl": d / "o/domains.domtbl", "hmmout": d / "o/domains.hmmscan.txt"},
              threads=THREADS, log=d / "log.txt")
    p = run_rule("prometheus", "prometheus_domain_scan", **ns)
    assert p.returncode != 0 and "run 'hmmpress" in (d / "log.txt").read_text()
    run(["hmmpress", "-f", str(d / "Pfam-A.hmm")])
    p = run_rule("prometheus", "prometheus_domain_scan", **ns)
    assert p.returncode == 0, p.stderr
    hits = utils.parse_domtbl(str(d / "o/domains.domtbl"))["G1__G1_1_5"]
    assert utils.architecture_tokens(hits, SPEC) == ["DUF4011", "DNA2", "REase", "DUF3320"]
    assert vp.parse_hmmscan_alignments(str(d / "o/domains.hmmscan.txt"))["G1__G1_1_5"]
    write(d / "cand.faa", "")                          # zero candidates: empty-input placeholder
    assert run_rule("prometheus", "prometheus_domain_scan", **ns).returncode == 0
    assert utils.parse_domtbl(str(d / "o/domains.domtbl")) == {}


@suite.test
def validation_tree_rule_real_tools():
    rnd = random.Random(7)

    def mut(s, f):
        s = list(s)
        for i in rnd.sample(range(len(s)), int(len(s) * f)):
            s[i] = rnd.choice("ACDEFGHIKLMNPQRSTVWY")
        return "".join(s)

    d = tmpdir("vtr_")
    prom = [("sp|A0A060RCK4|PROA_9BACT", "A0A060RCK4"), ("tr|A0A0A7PKR5|X:1", "A0A0A7PKR5"),
            ("tr|A0A073IP09|Y", "A0A073IP09")]
    para = [("tr|A0A1H8VVA9|HALO1", "A0A1H8VVA9"), ("tr|A0A2A2F1Z1|HALO2", "A0A2A2F1Z1")]
    # DUF4011 + N-terminal SF1 helicase region of each protein (homologous, short runtime)
    write(d / "prom.faa", "".join(f">{n} desc\n{by_acc(a)[:900]}\n" for n, a in prom))
    write(d / "para.faa", "".join(f">{n}\n{by_acc(a)[1150:2050]}\n" for n, a in para))
    write(d / "cand.faa", f">GA__GA_1_5\n{mut(by_acc('A0A060RCK4')[:900], 0.05)}*\n"
          f">GB__GB_2_9\n{mut(by_acc('A0A1H8VVA9')[1150:2050], 0.05)}*\n")
    out = d / "vt"
    ns = dict(input={"candidates": d / "cand.faa", "prom_ref": d / "prom.faa", "para_ref": d / "para.faa"},
              output={"tree": out / "validation.treefile", "labels": out / "reference_labels.tsv"},
              params={"outdir": out, "model": "LG+G4", "ufboot": 1000, "alrt": 1000, "seed": 12345},
              threads=THREADS, log=d / "log.txt")
    p = run_rule("prometheus", "prometheus_validation_tree", timeout=1500, **ns)
    assert p.returncode == 0, (d / "log.txt").read_text()[-2000:]
    labels = vp.load_ref_labels(str(out / "reference_labels.tsv"))
    tips = re.findall(r"[(,]([^():,;]+):", (out / "validation.treefile").read_text())
    assert set(labels) <= set(tips) and "tr_A0A0A7PKR5_X_1" in labels, (labels, tips)
    det = {}
    a = vp.assign_by_tree(str(out / "validation.treefile"), labels, 80, 95, det)
    assert a == {"GA__GA_1_5": "Prometheus", "GB__GB_2_9": "Paralog"}, (a, det)
    shutil.copy(out / "validation.treefile", d / "keep.treefile")
    shutil.copy(out / "reference_labels.tsv", d / "keep_labels.tsv")
    _STATE["vt"] = (d / "keep.treefile", d / "keep_labels.tsv")
    write(d / "nocand.faa", "")
    ns["input"]["candidates"] = d / "nocand.faa"
    p = run_rule("prometheus", "prometheus_validation_tree", **ns)
    assert p.returncode == 0 and (out / "validation.treefile").read_text().strip() == "();"
    write(d / "empty.faa", "")
    ns["input"]["candidates"], ns["input"]["para_ref"] = d / "cand.faa", d / "empty.faa"
    p = run_rule("prometheus", "prometheus_validation_tree", **ns)
    assert p.returncode != 0 and "are empty" in (d / "log.txt").read_text()


@suite.test
def plot_validation_tree_r():
    d = tmpdir("pvt_")
    tree, labels = _STATE.get("vt", (None, None))
    if tree is None:
        tree = write(d / "t.nwk", "((Q1:0.1,(Q2:0.1,GB__GB_2_9:0.05)99/100:0.1)100/100:0.5,"
                     "(P1:0.1,GA__GA_1_5:0.1)98/100:0.2,P2:0.3);")
        labels = write(d / "l.tsv", "P1\tPrometheus_ref\nP2\tPrometheus_ref\nQ1\tParalog_ref\nQ2\tParalog_ref\n")
    write(d / "cat.tsv", "locus_id\tprotein_id\tconfidence_tier\nGA__GA_1_5\tGA__GA_1_5\tA\n"
          "GB__GB_2_9\tGB__GB_2_9\tRejected\n")
    r = rscript("plot_validation_tree.R", "--tree", tree, "--labels", labels, "--catalog",
                d / "cat.tsv", "--plot-prefix", d / "vt")
    assert (d / "vt.png").stat().st_size > 1000 and (d / "vt.svg").exists()
    assert "rooted on the paralog outgroup" in r.stderr, r.stderr[-800:]
    rscript("plot_validation_tree.R", "--tree", write(d / "triv.nwk", "();\n"), "--labels", labels,
            "--catalog", d / "cat.tsv", "--plot-prefix", d / "ph")
    assert (d / "ph.png").exists()
    r = rscript("plot_validation_tree.R", "--tree", write(d / "bad.nwk", "((a,b"), "--labels", labels,
                "--catalog", d / "cat.tsv", "--plot-prefix", d / "bad", check=False)
    assert r.returncode != 0


# --------------------------------------------------------------------------- #
# 6. diversity: typology, clustering (real mmseqs), plots
# --------------------------------------------------------------------------- #
def typology_inputs(d: Path):
    domtbl, _ = real_hmmscan()
    hits = utils.parse_domtbl(str(domtbl))
    rows, lines = ["locus_id\tgenome_id\tarchitecture\tconfidence_tier"], []
    for i, (name, h) in enumerate(hits.items()):
        pid = f"G{i}__{utils.sanitize(name)}"
        rows.append(f"{pid}\tG{i}\t{'-'.join(utils.architecture_tokens(h, SPEC))}\t"
                    f"{'Rejected' if 'A0A0F0LN71' in name else 'A'}")
        lines += [ln.replace(name, pid, 1) for ln in open(domtbl)
                  if not ln.startswith("#") and ln.split()[3] == name]
    write(d / "cat.tsv", "\n".join(rows) + "\n")
    write(d / "x.domtbl", "".join(lines))


@suite.test
def architecture_typology_known_answers():
    d = tmpdir("typ_")
    typology_inputs(d)
    py("architecture_typology.py", "--domtbl", d / "x.domtbl", "--catalog", d / "cat.tsv",
       "--domains", SPEC_JSON, "--out", d / "typ.tsv")
    t = read_tsv(d / "typ.tsv").set_index("architecture")
    assert t.columns[:7].tolist() == ["n_loci", "example_locus", "has_DUF4011", "has_DUF3320",
                                      "has_DNA2", "has_nuclease", "has_REase"]
    assert t.index.tolist() == ["DUF4011-DNA2-REase-DUF3320", "DUF4011-DNA2-REase"]
    full, core = t.loc["DUF4011-DNA2-REase-DUF3320"], t.loc["DUF4011-DNA2-REase"]
    assert (full.n_loci, full.variant_class, full.has_DUF3320, full.has_nuclease) == \
        ("6", "canonical_full", "True", "False")
    assert (core.n_loci, core.variant_class, core.has_DUF3320, core.canonical_order) == \
        ("4", "canonical_core", "False", "True")
    assert (core.fusion_domains, core.n_with_fusion, core.n_tier_A, core.core_complete) == \
        ("Pkinase:2;ResIII+Helicase_C:2", "4", "4", "True")
    write(d / "r.domtbl", domline("MTES_1575", "PF18741.8", "G9__R", 10, 100)
          + domline("DUF4011", "PF13195.12", "G9__R", 200, 360)
          + domline("AAA_12", "PF13087.13", "G9__R", 500, 700)
          + domline("DUF4011", "PF13195.12", "G8__L", 20, 180))
    write(d / "r_cat.tsv", "locus_id\tgenome_id\tarchitecture\tconfidence_tier\n"
          "G9__R\tG9\tREase-DUF4011-DNA2\tB\nG8__L\tG8\tDUF4011\tC\n")
    py("architecture_typology.py", "--domtbl", d / "r.domtbl", "--catalog", d / "r_cat.tsv",
       "--domains", SPEC_JSON, "--out", d / "r_typ.tsv")
    r = read_tsv(d / "r_typ.tsv").set_index("architecture")
    assert (r.loc["REase-DUF4011-DNA2", "variant_class"], r.loc["REase-DUF4011-DNA2", "rearranged"]) == \
        ("rearranged", "True")
    assert (r.loc["DUF4011", "variant_class"], r.loc["DUF4011", "missing_required"]) == ("core_loss", "DNA2,REase")
    write(d / "bad_cat.tsv", "locus_id\tgenome_id\tarchitecture\tconfidence_tier\nG8__L\tG8\tDNA2\tC\n")
    expect_fail(py("architecture_typology.py", "--domtbl", d / "r.domtbl", "--catalog", d / "bad_cat.tsv",
                   "--domains", SPEC_JSON, "--out", d / "x.tsv", check=False), "disagree")
    write(d / "e_cat.tsv", "locus_id\tgenome_id\tarchitecture\tconfidence_tier\n")
    py("architecture_typology.py", "--domtbl", d / "r.domtbl", "--catalog", d / "e_cat.tsv",
       "--domains", SPEC_JSON, "--out", d / "e.tsv")
    assert read_tsv(d / "e.tsv").empty
    _STATE["typ"] = (d / "typ.tsv", d / "cat.tsv")


def renamed_fixture(d: Path) -> Path:
    utils.write_fasta([(f"G{i}__{utils.sanitize(n)}", s) for i, (n, s) in
                       enumerate(fixture_seqs().items())], str(d / "in.faa"))
    return d / "in.faa"


@suite.test
def mmseqs_cluster_real_and_rule():
    d = tmpdir("mm_")
    faa = renamed_fixture(d)
    ns = dict(input={"faa": faa}, output={"subfam": d / "sub.tsv", "summary": d / "summ.tsv"},
              params={"identities": "0.3,0.5,0.7,0.9", "coverage": 0.6, "cov_mode": 0,
                      "sensitivity": 7.5, "workdir": d / "work"},
              threads=THREADS, log=d / "log.txt")
    p = run_rule("diversity", "mmseqs_cluster", **ns)
    assert p.returncode == 0, (d / "log.txt").read_text()[-1500:]
    sub = read_tsv(d / "sub.tsv")
    assert sub.columns.tolist() == ["protein_id", "genome_id", "cluster_0.30", "cluster_0.50",
                                    "cluster_0.70", "cluster_0.90"]
    assert len(sub) == 11 and set(sub["protein_id"]) == set(utils.fasta_ids(str(faa)))
    c = sub.set_index("protein_id")
    halo = [x for x in c.index if "9EURY" in x]
    pk = [x for x in c.index if "9GAMM" in x or "9HYPH" in x]
    assert c.loc[halo[0], "cluster_0.50"] == c.loc[halo[1], "cluster_0.50"]   # 2.8-2.9 kaa pair
    assert c.loc[pk[0], "cluster_0.30"] == c.loc[pk[1], "cluster_0.30"]       # Pkinase fusions
    assert c["cluster_0.90"].nunique() == 11
    n30 = c["cluster_0.30"].nunique()
    assert n30 < c["cluster_0.50"].nunique() <= c["cluster_0.70"].nunique() <= 11
    summ = read_tsv(d / "summ.tsv").set_index("identity")
    assert (summ.loc["0.90", "n_clusters"], summ.loc["0.30", "n_clusters"]) == ("11", str(n30))
    _STATE["sub"] = d / "sub.tsv"


@suite.test
def mmseqs_failure_single_and_empty():
    d = tmpdir("mm2_")
    faa = renamed_fixture(d)
    fake = d / "fakebin"
    write(fake / "mmseqs", "#!/bin/bash\necho 'simulated crash' >&2\nexit 3\n").chmod(0o755)
    args = [PYTHON, str(SCRIPTS_PY / "cluster_diversity.py"), "--faa", str(faa), "--identities",
            "0.3", "--coverage", "0.6", "--workdir", str(d / "w"), "--out", str(d / "o.tsv")]
    e = env()
    e["PATH"] = str(fake) + os.pathsep + e["PATH"]
    expect_fail(subprocess.run(args, env=e, text=True, capture_output=True),
                "mmseqs easy-cluster failed", "simulated crash")
    e2 = env()
    e2["PATH"] = os.pathsep.join(x for x in e2["PATH"].split(os.pathsep)
                                 if not (Path(x) / "mmseqs").exists())
    expect_fail(subprocess.run(args, env=e2, text=True, capture_output=True), "mmseqs not on PATH")
    write(d / "one.faa", ">G1__G1_1_1\nMKV\n")
    py("cluster_diversity.py", "--faa", d / "one.faa", "--identities", "0.3,0.9", "--coverage",
       "0.6", "--workdir", d / "w1", "--out", d / "one.tsv")
    assert read_tsv(d / "one.tsv").values.tolist() == [["G1__G1_1_1", "G1", "G1__G1_1_1", "G1__G1_1_1"]]
    write(d / "none.faa", "")
    py("cluster_diversity.py", "--faa", d / "none.faa", "--identities", "0.30,0.50", "--coverage",
       "0.6", "--workdir", d / "w2", "--out", d / "none.tsv")
    t = read_tsv(d / "none.tsv")
    assert t.empty and t.columns.tolist() == ["protein_id", "genome_id", "cluster_0.30", "cluster_0.50"]
    expect_fail(py("cluster_diversity.py", "--faa", d / "one.faa", "--identities", "0.3,0.30",
                   "--coverage", "0.6", "--workdir", d / "w3", "--out", d / "x.tsv", check=False),
                "duplicates")


@suite.test
def rarefaction_exact_values():
    import numpy as np
    import plot_diversity as pdv
    # 4 genomes; subfamily a in {G1, G2}, b in {G3}; G4 carries nothing
    assert np.allclose(pdv.expected_richness([2, 1], 4, [0, 1, 2, 3, 4]), [0, 0.75, 4 / 3, 1.75, 2.0])
    assert abs(pdv.chao2([2, 1], 4) - 2.375) < 1e-9
    d = tmpdir("rar_")
    write(d / "sub.tsv", "protein_id\tgenome_id\tcluster_0.30\tcluster_0.90\n"
          "G1__p\tG1\ta\ta\nG2__q\tG2\ta\tq\nG3__r\tG3\tb\tr\n")
    write(d / "master.tsv", "genome_id\tgtdb_phylum\nG1\tp__A\nG2\tp__A\nG3\tp__B\nG4\tp__B\n")
    write(d / "typ.tsv", "architecture\tn_loci\n")
    write(d / "cat.tsv", "locus_id\tgenome_id\tarchitecture\tconfidence_tier\n"
          "G1__p\tG1\tDUF4011-DNA2-REase-DUF3320\tA\nG2__q\tG2\tDUF4011-DNA2-REase\tB\n"
          "G3__r\tG3\t\tC\nG4__s\tG4\tDUF4011\tRejected\n")
    py("plot_diversity.py", "--subfamilies", d / "sub.tsv", "--typology", d / "typ.tsv",
       "--master", d / "master.tsv", "--catalog", d / "cat.tsv", "--seed", "1",
       "--permutations", "50", "--curve", d / "curve.tsv", "--outdir", d / "plots")
    c = read_tsv(d / "curve.tsv")
    assert c["n_genomes"].tolist() == ["0", "1", "2", "3", "4"]
    assert np.allclose(c["expected_subfamilies"].astype(float), [0, 0.75, 4 / 3, 1.75, 2.0])
    assert float(c["perm_hi95"].iloc[-1]) == 2.0 and float(c["perm_lo95"].iloc[1]) <= 0.75
    for f in ("subfamily_rarefaction", "architecture_by_taxon"):
        assert (d / "plots" / f"{f}.png").stat().st_size > 1000
    write(d / "m2.tsv", "genome_id\tgtdb_phylum\nG1\tp__A\n")
    expect_fail(py("plot_diversity.py", "--subfamilies", d / "sub.tsv", "--typology", d / "typ.tsv",
                   "--master", d / "m2.tsv", "--outdir", d / "p2", check=False), "missing from the master")
    write(d / "esub.tsv", "protein_id\tgenome_id\tcluster_0.30\n")
    py("plot_diversity.py", "--subfamilies", d / "esub.tsv", "--typology", d / "typ.tsv",
       "--master", d / "master.tsv", "--outdir", d / "p3")
    assert (d / "p3" / "subfamily_rarefaction.svg").exists()


@suite.test
def diversity_extra_and_prometheus_plots():
    import plot_diversity_extra as pde
    if "typ" not in _STATE:
        architecture_typology_known_answers()
    typ, typ_cat = _STATE["typ"]
    assert pde.domain_categories(read_tsv(typ)) == ["DUF4011", "DNA2", "REase", "DUF3320"]
    d = tmpdir("pde_")
    sub = _STATE.get("sub") or write(d / "sub.tsv", "protein_id\tgenome_id\tcluster_0.30\nG1__a\tG1\tG1__a\n")
    write(d / "master.tsv", "genome_id\tgtdb_phylum\n" + "".join(f"G{i}\tp__X{i % 3}\n" for i in range(12)))
    py("plot_diversity_extra.py", "--typology", typ, "--catalog", typ_cat, "--master", d / "master.tsv",
       "--subfamilies", sub, "--outdir", d / "plots")
    for f in ("architecture_combinations", "architecture_by_phylum", "subfamily_identity_curve"):
        assert (d / "plots" / f"{f}.png").stat().st_size > 1000, f
    v = _STATE.get("validate_case")
    if v is None:
        v = tmpdir("ppv_")
        build_validate_case(v)
        py("validate_prometheus.py", *validate_args(v))
    py("plot_prometheus_validation.py", "--catalog", v / "cat.tsv", "--sensitivity", v / "sens.tsv",
       "--outdir", v / "plots")
    assert (v / "plots/prometheus_detection_summary.png").exists()
    assert (v / "plots/prometheus_architecture.svg").exists()
    write(v / "badtier.tsv", "locus_id\tconfidence_tier\tarchitecture\nx\tZ\tDUF4011\n")
    expect_fail(py("plot_prometheus_validation.py", "--catalog", v / "badtier.tsv", "--sensitivity",
                   v / "sens.tsv", "--outdir", v / "p2", check=False), "unexpected confidence_tier")


# --------------------------------------------------------------------------- #
# 7. real DefenseFinder 3.0.0 + models 3.1.0 / CasFinder 3.1.0
# --------------------------------------------------------------------------- #
@suite.test
def real_defensefinder_preflight_run_parse_context():
    if df_bin() is None:
        print("  skip  real DefenseFinder: not installed (set PROM_DF_BIN)")
        return
    models = df_models()
    if models is None:
        print("  skip  real DefenseFinder: pinned models unavailable (set PROM_DF_MODELS)")
        return
    d = tmpdir("dfreal_")
    e = df_env()
    p = subprocess.run([PYTHON, str(SCRIPTS_PY / "check_defensefinder.py"), "--run", "True",
                        "--software-version", "3.0.0", "--models-pin", "3.1.0",
                        "--casfinder-pin", "3.1.0", "--models-dir", str(models),
                        "--positive-control", "true", "--workdir", str(d / "pc"),
                        "--threads", str(THREADS),
                        "--out", str(d / "ok.tsv")], env=e, text=True, capture_output=True,
                       timeout=900)
    assert p.returncode == 0, p.stderr[-2000:]
    kv = dict(l.split("\t", 1) for l in (d / "ok.tsv").read_text().splitlines()[1:])
    assert kv["positive_control"].startswith("called Prometheus")
    assert kv["definition_grammar_versions"].count("'2.0'") == 2
    # two-contig genome: ProA homologue 2 genes from a real Septu system (DefenseFinder's
    # own test proteins) on contig 1; the ProA consensus alone on contig 2
    rnd = random.Random(3)
    with urllib.request.urlopen("https://raw.githubusercontent.com/mdmparis/defense-finder/"
                                "v3.0.0/tests/data/prot/df_test_prot.faa", timeout=120) as r:
        write(d / "df_test_prot.faa", r.read().decode())
    t = dict(utils.read_fasta(str(d / "df_test_prot.faa")))
    filler = lambda: "M" + "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(300))
    c1 = [filler(), filler(), t["CATH005.0523.00001.001C_00311"], t["CATH005.0523.00001.001C_00312"],
          filler(), by_acc("A0A060RCK4"), filler()]
    c2 = [filler() for _ in range(15)] + [by_acc("Prometheus__ProA-consensus")]
    utils.write_fasta([(f"GX_1_{i}", s) for i, s in enumerate(c1, 1)]
                      + [(f"GX_2_{i}", s) for i, s in enumerate(c2, 1)], str(d / "genes/GX/GX.faa"))
    out = d / "defense/GX"
    ns = dict(input={"faa": d / "genes/GX/GX.faa", "models_ok": d / "ok.tsv"},
              wildcards={"genome": "GX"},
              output={"systems": out / "GX_defense_finder_systems.tsv",
                      "genes": out / "GX_defense_finder_genes.tsv"},
              params={"run": "True", "provided": "", "outdir": out,
                      "models_opt": f"--models-dir {models}", "sys_header": "", "genes_header": ""},
              threads=THREADS, log=d / "df.log")
    p = run_rule("prometheus", "defensefinder", base_env=e, timeout=900, **ns)
    assert p.returncode == 0, (d / "df.log").read_text()[-2500:]
    assert sorted(read_tsv(out / "GX_defense_finder_systems.tsv")["type"]) == ["Prometheus", "Prometheus", "Septu"]
    assert set(read_tsv(out / "GX_defense_finder_genes.tsv")["replicon"]) == {"GX_1", "GX_2"}  # gembase
    write(d / "samples.tsv", "genome_id\nGX\n")
    py("combine_proteomes.py", "--samples", d / "samples.tsv", "--genes-dir", d / "genes",
       "--out-faa", d / "all.faa", "--out-map", d / "map.tsv")
    py("parse_defensefinder.py", "--df-dir", d / "defense", "--proteome", d / "all.faa", "--pmap",
       d / "map.tsv", "--samples", d / "samples.tsv", "--system", "Prometheus",
       "--out-candidates", d / "cand.faa", "--out-sources", d / "src.tsv", "--out-raw", d / "raw.tsv")
    assert utils.fasta_ids(str(d / "cand.faa")) == ["GX__GX_1_6", "GX__GX_2_16"]
    dfg = vp.load_df_genes(str(d / "defense"), "GX")
    own = lambda lid: set(dfg.loc[dfg["hit_id"] == lid, "sys_id"])
    c1ctx = vp.defense_context("GX", "GX_1_6", own("GX_1_6"), "Prometheus", 10, dfg, None)
    c2ctx = vp.defense_context("GX", "GX_2_16", own("GX_2_16"), "Prometheus", 10, dfg, None)
    assert c1ctx[0] is True and c1ctx[2] == 1 and c1ctx[3].startswith("Septu:")
    assert c2ctx[0] is False and c2ctx[2] == 0


if __name__ == "__main__":
    suite.main(TOKEN)
