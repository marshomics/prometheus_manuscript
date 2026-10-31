"""Synthetic test suite for Stage 5 (evolution): bash tests/run_leaf.sh evolution

Known-answer tests for every owned script, their empty-input and failure behaviour, and
real-tool runs of mafft, famsa, clipkit, seqkit, iqtree2, pal2nal, hyphy (pinned 2.5.62 when
the prom-evolution-hyphy env exists, plus the PATH version), R/phytools and GeneRax (bioconda
2.1.3; see test_generax_real_binary for what the sandboxed binary can and cannot do). ALE is
not installable here: its CLI/output contract is exercised with a mock built from the ALE
source. At most 2 threads per tool.
"""
import json
import os
import random
import re
import shutil
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _env import (Suite, py, rscript, run, tmpdir, FIXTURES, have, read_tsv, REPO,  # noqa: E402
                  SCRIPTS_PY, SCRIPTS_R, PYTHON, RSCRIPT, env, _scratch_base)

sys.path.insert(0, str(SCRIPTS_PY))
# MAFFT needs an ABSOLUTE TMPDIR; a sandbox may set a relative one (./.tmp). Use the same
# absolute scratch base as tmpdir() instead (never a directory relative to the cwd).
if not os.path.isabs(os.environ.get("TMPDIR", "")):
    os.environ["TMPDIR"] = _scratch_base()
suite = Suite("evolution")
T = "2"                                                      # threads per tool
DOMAINS = {"DUF4011": {"names": ["DUF4011"], "accessions": ["PF13195"], "required": True},
           "DNA2": {"names": ["AAA_11", "AAA_12", "AAA_19", "AAA_30"],
                    "accessions": ["PF13086", "PF13087", "PF13245", "PF13604"], "required": True},
           "REase": {"names": ["MTES_1575"], "accessions": ["PF18741"], "required": True},
           "DUF3320": {"names": ["DUF3320"], "accessions": ["PF11784"], "required": False}}
DJSON = json.dumps(DOMAINS)
_TB = os.environ.get("PROM_TOOLS_BIN", "")
ENVS = Path(_TB).resolve().parents[1] if _TB else Path("/nonexistent")
HYPHY_PINNED = ENVS / "prom-evolution-hyphy" / "bin"
GENERAX_DIR = ENVS / "prom-evolution-generax" / "bin"
OWNED_PY = ["domain_trees.py", "codon_align.py", "run_reconciliation.py", "run_selection.py",
            "plot_selection.py", "plot_evolution_extra.py"]
OWNED_R = ["ancestral_states.R", "plot_tanglegram.R"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def write(path, text):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(text)
    return str(path)


def fasta(path, recs):
    return write(path, "".join(f">{n}\n{s}\n" for n, s in recs))


def read_fa(path):
    import utils
    return list(utils.read_fasta(str(path)))


def must_fail(cmd, pattern=None, env_over=None, cwd=None):
    e = env()
    if env_over:
        e.update(env_over)
    p = subprocess.run(cmd, cwd=cwd or REPO, env=e, text=True, capture_output=True,
                       shell=isinstance(cmd, str), executable="/bin/bash" if isinstance(cmd, str) else None)
    assert p.returncode != 0, f"expected failure, got exit 0: {cmd}\n{p.stdout[-800:]}\n{p.stderr[-800:]}"
    if pattern:
        assert re.search(pattern, p.stdout + p.stderr), f"no /{pattern}/ in:\n{(p.stdout + p.stderr)[-1500:]}"
    return p


def path_without(*tools):
    """PATH with every directory that provides one of `tools` removed."""
    return os.pathsep.join(d for d in env()["PATH"].split(os.pathsep)
                           if not any(os.path.exists(os.path.join(d, t)) for t in tools))


def domtbl_line(target, acc, query, qlen, start, end):
    """One hmmscan --domtblout line (22 fixed columns + description)."""
    f = [target, acc, "100", query, "-", str(qlen), "1e-30", "100.0", "0.1", "1", "1",
         "1e-30", "1e-30", "100.0", "0.1", "1", "100", str(start), str(end), str(start),
         str(end), "0.95", "desc"]
    return " ".join(f) + "\n"


STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in STOPS]


def sim_codons(tree_newick, n_codons, seed, rate=1.5):
    import dendropy
    rng = random.Random(seed)
    t = dendropy.Tree.get(data=tree_newick, schema="newick", preserve_underscores=True)
    root = [rng.choice(CODONS) for _ in range(n_codons)]
    out = {}

    def ev(nd, seq):
        if nd.parent_node is not None:
            s = list(seq)
            for i in range(n_codons):
                if rng.random() < (nd.edge.length or 0) * rate:
                    c = list(s[i]); c[rng.randrange(3)] = rng.choice("ACGT"); c = "".join(c)
                    if c not in STOPS:
                        s[i] = c
            seq = s
        if nd.is_leaf():
            out[nd.taxon.label] = "".join(seq)
        for ch in nd.child_nodes():
            ev(ch, seq)
    ev(t.seed_node, root)
    return out


def fixture_family(d):
    """Fixture ProA proteins renamed to contract ids + the matching domtbl; returns paths."""
    import utils
    recs = list(utils.read_fasta(str(FIXTURES / "proA_homologs.faa")))
    ids = {n: f"G{i}__G{i}_1_{i}" for i, (n, _) in enumerate(recs, 1)}
    prot = fasta(Path(d) / "prot.faa", [(ids[n], s.rstrip("*")) for n, s in recs])
    lines = []
    for line in open(FIXTURES / "proA_homologs.pfam_subset.domtbl"):
        if not line.startswith("#"):
            f = line.split()
            if len(f) > 3 and f[3] in ids:
                line = line.replace(f[3], ids[f[3]], 1)
        lines.append(line)
    dom = write(Path(d) / "dom.tbl", "".join(lines))
    return prot, dom, ids


# --------------------------------------------------------------------------- #
# static checks
# --------------------------------------------------------------------------- #
@suite.test
def test_py_compile_owned():
    for f in OWNED_PY:
        run([PYTHON, "-m", "py_compile", str(SCRIPTS_PY / f)])
    run([PYTHON, "-m", "py_compile", __file__])


@suite.test
def test_r_parse_owned():
    for f in OWNED_R:
        run([RSCRIPT, "-e", f"invisible(parse(file='{SCRIPTS_R / f}'))"])


@suite.test
def test_rules_hygiene():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"), str(REPO / "workflow/rules/evolution.smk")])
    assert "RULES HYGIENE OK" in p.stdout


# --------------------------------------------------------------------------- #
# family_align / family_tree rule shell (rendered from evolution.smk)
# --------------------------------------------------------------------------- #
def rule_shell(rule):
    txt = (REPO / "workflow/rules/evolution.smk").read_text()
    m = re.search(rf"^rule {rule}:\n(.*?)(?=^rule |\Z)", txt, re.S | re.M)
    sh = re.search(r'shell:\n\s+r"""\n(.*?)"""', m.group(1), re.S).group(1)
    return textwrap.dedent(sh)


def render(sh, values):
    def sub(m):
        return str(values[m.group(1)])
    out = re.sub(r"(?<!\{)\{([a-z_.]+)\}(?!\})", sub, sh)
    return out.replace("{{", "{").replace("}}", "}")


def rule_paths(rule, root):
    """{'output.tree': Path, 'params.workdir': Path, 'log': Path, ...} for every f-string path the
    rule declares, with RESULTS/LOGS under `root`, so that path RELATIONSHIPS (e.g. an output that
    is also a tool's own <prefix> output) are tested exactly as written in evolution.smk."""
    txt = (REPO / "workflow/rules/evolution.smk").read_text()
    body = re.search(rf"^rule {rule}:\n(.*?)(?=^rule |\Z)", txt, re.S | re.M).group(1)
    out, section = {}, None
    for line in body.splitlines():
        m = re.match(r"^    (input|output|params|log):\s*$", line)
        if m:
            section = m.group(1)
            continue
        if re.match(r"^    \w+:", line):
            section = None
            continue
        if section == "log":
            m = re.match(r'^\s+f"([^"]+)",?\s*$', line)
            if m:
                out["log"] = m.group(1)
        elif section:
            m = re.match(r'^\s+(\w+)=f"([^"]+)",?\s*$', line)
            if m:
                out[f"{section}.{m.group(1)}"] = m.group(2)
    return {k: Path(v.replace("{RESULTS}", str(root / "results")).replace("{LOGS}", str(root / "logs")))
            for k, v in out.items()}


def align_values(d, faa, tool="mafft", linsi_max=200, trim="clipkit"):
    return {"input.faa": faa, "output.aln": d / "a.aln", "output.trim": d / "a.trim.aln",
            "output.method": d / "method.txt", "params.tool": tool, "params.linsi_max": linsi_max,
            "params.trim": trim, "params.clipkit_mode": "kpic-smart-gap",
            "params.clean": d / "clean.faa", "threads": T, "log": d / "log.txt"}


@suite.test
def test_family_align_shell_linsi_and_stop_stripping():
    d = tmpdir("evo_align_")
    recs = read_fa(FIXTURES / "proA_homologs.faa")[:6]
    faa = fasta(d / "in.faa", [(f"G{i}__G{i}_1_{i}", s.rstrip("*") + "*") for i, (n, s) in enumerate(recs, 1)])
    run(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_align"), align_values(d, faa))])
    aln = read_fa(d / "a.aln")
    assert len(aln) == 6 and len({len(s) for _, s in aln}) == 1, "not an alignment"
    assert all("*" not in s for _, s in aln)
    assert [s.replace("-", "").upper() for _, s in aln] == [s.rstrip("*").upper() for _, s in recs], \
        "aligned residues differ from the input"
    meth = (d / "method.txt").read_text()
    assert "L-INS-i" in meth and "n=6" in meth and "ClipKIT" in meth, meth
    trim = read_fa(d / "a.trim.aln")
    assert len(trim) == 6 and len(trim[0][1]) <= len(aln[0][1])


@suite.test
def test_family_align_shell_famsa_above_threshold_empty_and_bad_tool():
    d = tmpdir("evo_align2_")
    recs = read_fa(FIXTURES / "proA_homologs.faa")[:5]
    faa = fasta(d / "in.faa", [(f"G{i}__G{i}_1_{i}", s[:300]) for i, (n, s) in enumerate(recs, 1)])
    v = align_values(d, faa, linsi_max=3, trim="none")
    run(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_align"), v)])
    assert (d / "method.txt").read_text().startswith("FAMSA"), "n > linsi_max must use FAMSA"
    assert len({len(s) for _, s in read_fa(d / "a.aln")}) == 1
    assert (d / "a.aln").read_text() == (d / "a.trim.aln").read_text(), "trim_tool none must copy"
    v["input.faa"] = write(d / "empty.faa", "")
    run(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_align"), v)])
    assert (d / "a.aln").stat().st_size == 0 and "none" in (d / "method.txt").read_text()
    v["input.faa"], v["params.tool"] = faa, "bogus"
    must_fail(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_align"), v)],
              "unknown evolution.align_tool")


def tree_values(d, trim, model="LG"):
    v = rule_paths("family_tree", d)
    assert {"output.tree", "output.ufboot", "params.workdir", "log"} <= set(v), v
    v.update({"input.trim": trim, "params.model": model, "params.ufboot": 1000, "params.alrt": 1000,
              "params.seed": 1, "threads": T})
    return v


@suite.test
def test_family_tree_shell_small_input_placeholder_and_failure():
    d = tmpdir("evo_tree_")
    v = tree_values(d, fasta(d / "t.aln", [("a", "MKV"), ("b", "MKV"), ("c", "MRV")]))
    run(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_tree"), v)])
    assert v["output.tree"].read_text().strip() == "();"
    v["input.trim"] = fasta(d / "bad.aln", [("a", "MKV"), ("b", "MKV"), ("c", "MRV"), ("d", "MRVVVVV")])
    must_fail(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_tree"), v)])


@suite.test
def test_family_tree_real_iqtree_with_declared_paths():
    """Regression (cross-domain review): output.tree/.ufboot ARE IQ-TREE's own <prefix> outputs.
    The pristine `cp prefix.treefile output.tree || echo "();" > output.tree` copied a file onto
    itself (cp exits 1), so the ML tree was replaced by "();" on every run. The rule shell runs
    here with the paths exactly as declared, and the tree must have one tip per sequence,
    including identical sequences (IQ-TREE ignores them in the search and re-adds them)."""
    import dendropy
    d = tmpdir("evo_tree_real_")
    rng = random.Random(3)
    A = "ACDEFGHIKLMNPQRSTVWY"
    base = "".join(rng.choice(A) for _ in range(200))
    mut = lambda x, p: "".join(c if rng.random() > p else rng.choice(A) for c in x)  # noqa: E731
    c1, c2 = mut(base, 0.1), mut(base, 0.1)
    seqs = [mut(c1, 0.05), mut(c1, 0.05), mut(c2, 0.05), mut(c2, 0.05), mut(base, 0.2)]
    seqs += [seqs[0]] * 3                                          # three copies identical to G1
    ids = [f"G{i}__G{i}_1_1" for i in range(1, len(seqs) + 1)]
    v = tree_values(d, fasta(d / "t.aln", list(zip(ids, seqs))))
    run(["/bin/bash", "-euo", "pipefail", "-c", render(rule_shell("family_tree"), v)])
    t = dendropy.Tree.get(data=v["output.tree"].read_text(), schema="newick", preserve_underscores=True)
    assert sorted(lf.taxon.label for lf in t.leaf_node_iter()) == sorted(ids), "tips != input sequences"
    uf = v["output.ufboot"].read_text().splitlines()
    assert len(uf) == 1000 and all(x.count(",") + 1 == len(ids) for x in uf), "UFBoot trees incomplete"
    # a tree with fewer tips than sequences must fail the rule (a stale/foreign treefile)
    write(d / "fake_iqtree" / "iqtree2", "#!/bin/bash\nwhile [ $# -gt 0 ]; do [ \"$1\" = --prefix ] && p=$2; "
          "shift; done\necho '(G1__G1_1_1,G2__G2_1_1,G3__G3_1_1);' > \"$p.treefile\"; echo '();' > \"$p.ufboot\"\n")
    os.chmod(d / "fake_iqtree" / "iqtree2", 0o755)
    sh = f'export PATH="{d / "fake_iqtree"}:$PATH"\n' + render(rule_shell("family_tree"), v)
    must_fail(["/bin/bash", "-euo", "pipefail", "-c", sh], "has 3 tips but the alignment has 8 sequences")


# --------------------------------------------------------------------------- #
# codon_align.py
# --------------------------------------------------------------------------- #
PEP = [("a", "MKV-LE"), ("b", "MKVALE"), ("c", "MRVALE")]
CDS = [("a", "ATGAAAGTTCTGGAATAA"), ("b", "GTGAAAGTTGCGCTGGAATAG"), ("c", "TTGCGTGTTGCGCTGGAATGA")]


def codon_case(d, pep=PEP, cds=CDS, prot=None):
    a = fasta(d / "pep.aln", pep)
    c = fasta(d / "cds.fna", cds)
    p = fasta(d / "prot.faa", prot if prot is not None else [(n, s.replace("-", "") + "*") for n, s in pep])
    return a, c, p


def codon_cmd(a, c, p, out, *extra):
    return [PYTHON, str(SCRIPTS_PY / "codon_align.py"), "--prot-aln", a, "--cds", c, "--proteins", p,
            "--out", str(out), *extra]


@suite.test
def test_codon_align_known_answer_gtg_ttg_starts_terminal_stops():
    d = tmpdir("evo_codon_")
    a, c, p = codon_case(d)
    run(codon_cmd(a, c, p, d / "out.fa", "--report", str(d / "rep.tsv")))
    assert read_fa(d / "out.fa") == [("a", "ATGAAAGTT---CTGGAA"), ("b", "GTGAAAGTTGCGCTGGAA"),
                                     ("c", "TTGCGTGTTGCGCTGGAA")]
    assert list(read_tsv(d / "rep.tsv")["status"]) == ["kept"] * 3


@suite.test
def test_codon_align_table1_rejects_gtg_initiator():
    d = tmpdir("evo_codon1_")
    a, c, p = codon_case(d, PEP + [("e", "MKVALE")], CDS + [("e", "ATGAAAGTTGCGCTGGAATAA")])
    run(codon_cmd(a, c, p, d / "out.fa", "--report", str(d / "rep.tsv"), "--codon-table", "1",
                  "--max-drop-frac", "0.5"))
    rep = read_tsv(d / "rep.tsv")
    st = dict(zip(rep["id"], rep["status"]))
    assert st["b"] == "translation_mismatch" and st["c"] == "kept", st   # TTG is a table-1 initiator
    assert [n for n, _ in read_fa(d / "out.fa")] == ["a", "c", "e"]


@suite.test
def test_codon_align_failures_are_loud():
    d = tmpdir("evo_codonf_")
    a, c, p = codon_case(d, cds=CDS[:2])                       # expected locus without CDS
    must_fail(codon_cmd(a, c, p, d / "o.fa"), "could not be codon-aligned")
    a, c, p = codon_case(d, cds=[CDS[0], ("b", "GTGAAATAGGCGCTGGAATAG"), CDS[2]])
    must_fail(codon_cmd(a, c, p, d / "o.fa"), "internal_stop")
    a2 = fasta(d / "unal.aln", [("a", "MKVLE"), ("b", "MKVALE"), ("c", "MRVALE")])
    must_fail(codon_cmd(a2, c, p, d / "o.fa"), "not an alignment")
    a, c, p = codon_case(d, prot=[("a", "MKVLE"), ("b", "MKVALW"), ("c", "MRVALE")])
    must_fail(codon_cmd(a, c, p, d / "o.fa"), "protein_mismatch")
    a, c, p = codon_case(d)
    must_fail(codon_cmd(a, c, p, d / "o.fa"), "pal2nal", env_over={"PATH": path_without("pal2nal.pl", "pal2nal")})


@suite.test
def test_codon_align_unexpected_refs_dropped_and_multi_files():
    d = tmpdir("evo_codonm_")
    a = fasta(d / "pep.aln", PEP + [("ref_DNA2", "MKVALE")])
    c1, c2 = fasta(d / "c1.fna", CDS[:2]), fasta(d / "c2.fna", CDS[2:])
    p1 = fasta(d / "p1.faa", [("a", "MKVLE*"), ("b", "MKVALE*")])
    p2 = fasta(d / "p2.faa", [("c", "MRVALE*")])
    r = py("codon_align.py", "--prot-aln", a, "--cds", c1, c2, "--proteins", p1, "--proteins", p2,
           "--out", d / "out.fa", "--report", d / "rep.tsv")
    assert [n for n, _ in read_fa(d / "out.fa")] == ["a", "b", "c"]
    assert "WARNING: 1 aligned record(s) NOT listed in --proteins" in r.stderr
    rep = read_tsv(d / "rep.tsv").set_index("id")
    assert rep.loc["ref_DNA2", "status"] == "no_cds" and rep.loc["ref_DNA2", "in_proteins"] == "0"


@suite.test
def test_codon_align_empty_input():
    d = tmpdir("evo_codone_")
    py("codon_align.py", "--prot-aln", write(d / "pep.aln", ""), "--cds", fasta(d / "c.fna", CDS),
       "--proteins", write(d / "p.faa", ""), "--out", d / "out.fa")
    assert (d / "out.fa").exists() and (d / "out.fa").stat().st_size == 0


# --------------------------------------------------------------------------- #
# run_selection.py
# --------------------------------------------------------------------------- #
SEL_TREE = ("(((G1__G1_1_1:0.1,G2__G2_1_1:0.1):0.1,G3__G3_1_1:0.2):0.1,((G4__G4_1_1:0.1,"
            "G5__G5_1_1:0.1):0.1,(G6__G6_1_1:0.2,G7__G7_1_1:0.15):0.1):0.1);")


def selection_inputs(d, n_codons=30):
    seqs = sim_codons(SEL_TREE, n_codons, seed=11)
    names = list(seqs)
    gap = 25 if n_codons >= 25 else n_codons // 2         # G5 lacks this codon (1-based)
    s5 = seqs["G5__G5_1_1"]
    seqs["G5__G5_1_1"] = s5[:3 * (gap - 1)] + "---" + s5[3 * gap:]
    fasta(d / "codon.fasta", [(n, seqs[n]) for n in names])
    write(d / "gene.nwk", SEL_TREE + "\n")
    write(d / "catalog.tsv", "locus_id\tgenome_id\tcatalytic_status\n" + "".join(
        f"{n}\t{n.split('__')[0]}\t{'degenerate' if n.startswith('G7') else 'intact'}\n" for n in names))
    write(d / "master.tsv", "genome_id\tmethanogen\tis_anme\n" + "".join(
        f"G{i}\t{'True' if i <= 3 else 'False'}\t{'True' if i == 6 else 'False'}\n" for i in range(1, 8)))
    lines = []
    for n in names:
        L = n_codons - (1 if n.startswith("G5") else 0)
        lines.append(domtbl_line("DUF4011", "PF13195.9", n, L, 1, 15))
        lines.append(domtbl_line("MTES_1575", "PF18741.4", n, L, 16, L))
    write(d / "dom.tbl", "".join(lines))


def selection_cmd(d, methods="FEL,MEME,BUSTED,aBSREL,RELAX", extra=()):
    return [PYTHON, str(SCRIPTS_PY / "run_selection.py"), "--codon", str(d / "codon.fasta"),
            "--gene-tree", str(d / "gene.nwk"), "--catalog", str(d / "catalog.tsv"),
            "--master", str(d / "master.tsv"), "--domtbl", str(d / "dom.tbl"), "--domains", DJSON,
            "--methods", methods, "--threads", T, "--workdir", str(d / "wd"),
            "--out-summary", str(d / "sum.tsv"), "--out-sites", str(d / "sites.tsv"), *extra]


def hyphy_env():
    e = env()
    if (HYPHY_PINNED / "hyphy").exists():
        e["PATH"] = str(HYPHY_PINNED) + os.pathsep + e["PATH"]
    return e


# Long real-tool runs (HyPhy x5 methods, per-domain IQ-TREE) are started in the background when
# the suite runs as a script (prefetch(), before suite.main), so they overlap with the fast tests
# instead of adding to the wall time. Every tool still uses at most T threads. When the module is
# imported (single tests), background() simply runs the job synchronously.
PREFETCH = {}
BG_JOBS = {}


def bg_job(key):
    def deco(make):
        BG_JOBS[key] = make
        return make
    return deco


def background(key):
    """(workdir, completed-process-like) of job `key`."""
    from types import SimpleNamespace
    if key in PREFETCH:
        d, proc, out_p, err_p = PREFETCH.pop(key)
        rc = proc.wait()
        return d, SimpleNamespace(returncode=rc, stdout=Path(out_p).read_text(), stderr=Path(err_p).read_text())
    d, cmd, e = BG_JOBS[key]()
    return d, subprocess.run(cmd, cwd=REPO, env=e, text=True, capture_output=True)


def prefetch():
    for key, make in BG_JOBS.items():
        try:
            d, cmd, e = make()
        except Exception as exc:  # noqa: BLE001 -- the test itself re-runs the setup and reports it
            print(f"[prefetch] {key}: setup failed ({exc!r}); the test will run it synchronously")
            continue
        out_p, err_p = d / "bg.stdout", d / "bg.stderr"
        with open(out_p, "w") as fo, open(err_p, "w") as fe:
            PREFETCH[key] = (d, subprocess.Popen(cmd, cwd=REPO, env=e, stdout=fo, stderr=fe, text=True),
                             out_p, err_p)


@bg_job("selection_all")
def _job_selection_all():
    d = tmpdir("evo_sel_")
    selection_inputs(d)
    return d, selection_cmd(d), hyphy_env()


@suite.test
def test_selection_prunes_tips_without_codons_when_not_intact_only():
    """Driver/cross-domain request: the gene tree is pruned whenever its tips and the codon
    records differ, not only with --intact-only True (reference/paralog tips have no CDS)."""
    assert have("hyphy"), "hyphy not on PATH"
    import dendropy
    d = tmpdir("evo_selprune_")
    selection_inputs(d, n_codons=12)
    write(d / "gene.nwk", "(" + SEL_TREE.rstrip().rstrip(";") + ":0.1,(REF1__r_1_1:0.3,REF2__r_1_2:0.3):0.2);\n")
    e = env()
    if (HYPHY_PINNED / "hyphy").exists():
        e["PATH"] = str(HYPHY_PINNED) + os.pathsep + e["PATH"]
    p = subprocess.run(selection_cmd(d, methods="FEL", extra=("--intact-only", "False")), cwd=REPO, env=e,
                       text=True, capture_output=True)
    assert p.returncode == 0, p.stderr[-3000:]
    t = dendropy.Tree.get(path=str(d / "wd" / "gene.analysed.nwk"), schema="newick", preserve_underscores=True)
    codon_ids = sorted(n for n, _ in read_fa(d / "codon.fasta"))
    assert sorted(lf.taxon.label for lf in t.leaf_node_iter()) == codon_ids, "tree not pruned to the codon records"
    fel = read_tsv(d / "sum.tsv").set_index("statistic").loc["sites_positive"]
    assert "2 gene-tree tip(s) without a codon record pruned" in fel["note"], fel["note"]
    assert "intact-only" not in fel["note"] and fel["n_sites_or_branches"] == "12"


@suite.test
def test_selection_hyphy_seeded_busted_reproducible():
    """HyPhy seeds its RNG from time+pid; BUSTED's random initial grid then makes unseeded runs
    differ (seen with 2.5.62). run_selection.py passes ENV='SetParameter(RANDOM_SEED,s,0);',
    so two runs with the same --seed must give bit-identical BUSTED results."""
    assert have("hyphy"), "hyphy not on PATH"
    e, runs = hyphy_env(), []
    for k in (1, 2):                                 # the two runs execute concurrently
        d = tmpdir(f"evo_seed{k}_")
        selection_inputs(d, n_codons=12)
        runs.append((d, subprocess.Popen(selection_cmd(d, methods="BUSTED", extra=("--seed", "4242")), cwd=REPO,
                                         env=e, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)))
    res = []
    for d, proc in runs:
        _, err = proc.communicate()
        assert proc.returncode == 0, err[-3000:]
        assert "ENV=SetParameter(RANDOM_SEED,4242,0);" in err, err[-1500:]
        s = read_tsv(d / "sum.tsv").set_index("method")
        res.append((s.loc["BUSTED", "value"], s.loc["BUSTED", "p_value"]))
    assert res[0] == res[1], f"seeded BUSTED runs differ: {res}"


@suite.test
def test_selection_hyphy_path_version_fel_meme_parse():
    """The PATH HyPhy (2.5.103 here) inserts a MEME column; name-based parsing must still work."""
    d = tmpdir("evo_sel2_")
    selection_inputs(d, n_codons=20)
    write(d / "dom.tbl", "# none\n")
    run(selection_cmd(d, methods="FEL,MEME", extra=("--intact-only", "False")))
    sites = read_tsv(d / "sites.tsv")
    assert len(sites) == 20 and set(sites["domain"]) == {"NA"}
    assert sites["meme_p_value"].astype(float).between(0, 1).all()
    assert "n=7 loci" in read_tsv(d / "sum.tsv").set_index("method").loc["MEME", "note"]


@suite.test
def test_selection_parsers_known_answers():
    import run_selection as rs
    import dendropy
    import utils
    assert rs._norm("&beta;<sup>+</sup>") == "beta+" and rs._norm("p<sup>1</sup>") == "p1"
    assert rs._norm("FEL &alpha;") == "fel alpha" and rs._norm("alpha=beta") == "alpha=beta"
    hdr = [["alpha", ""], ["beta", ""], ["alpha=beta", ""], ["LRT", ""], ["p-value", ""], ["Total branch length", ""]]
    js = {"MLE": {"headers": hdr, "content": {"1": [[1, 3, 2, 5, 0.01, 1], [0, 2, 1, 3, 0.2, 1]],
                                              "0": [[2, 0.5, 1, 4, 0.05, 1], [0, 0, 0, 0, 1, 0]]}},
          "data partitions": {"0": {"coverage": [[0, 1]]}, "1": {"coverage": [[2, 3]]}}}
    fel = rs.parse_fel(js, 0.1)
    assert fel[0]["class"] == "negative" and fel[2]["class"] == "positive" and fel[3]["class"] == "neutral"
    assert fel[2]["ratio"] == 3 and fel[3]["ratio"] == float("inf") and fel[1]["ratio"] != fel[1]["ratio"]
    try:
        rs.parse_fel({"MLE": {"headers": [["beta", ""]], "content": {"0": [[1]]}}}, 0.1)
        raise AssertionError("a missing alpha column must raise")
    except SystemExit:
        pass
    t = dendropy.Tree.get(data="((A,B),(C,D),(E,F));", schema="newick", rooting="force-unrooted")
    st = {"A": "M", "B": "M", "C": "N", "D": "N", "E": "M", "F": "U"}
    key = lambda nd: nd.taxon.label if nd.is_leaf() else "|".join(sorted(x.taxon.label for x in nd.leaf_iter()))  # noqa: E731
    lab = {key(nd): v for nd, v in rs.relax_labels(t, st, "clade").items()}
    assert lab == {"A": "Test", "B": "Test", "A|B": "Test", "C": "Reference", "D": "Reference",
                   "C|D": "Reference", "E": "Test"}, lab
    lab_t = rs.relax_labels(t, st, "terminal")
    assert all(nd.is_leaf() for nd in lab_t) and len(lab_t) == 5
    recs = [("s1", "ATG---AAA"), ("s2", "ATGCCCAAA"), ("s3", "ATGCCC---")]
    h = lambda e: [{"domain": "DUF4011", "accession": "PF13195.1", "start": 1, "end": e, "env_from": 1, "env_to": e}]  # noqa: E731
    cd = rs.column_domains(recs, {"s1": h(2), "s2": h(1), "s3": h(2)}, utils.load_domain_spec(DJSON))
    assert cd[0] == ("DUF4011", 1.0, 3) and cd[1] == ("NA", 0.5, 2) and cd[2] == ("DUF4011", 0.5, 2), cd


@suite.test
def test_selection_disabled_small_and_failures():
    d = tmpdir("evo_self_")
    selection_inputs(d, n_codons=12)
    run(selection_cmd(d, extra=("--run", "False")))
    s = read_tsv(d / "sum.tsv")
    assert set(s["method"]) == {"FEL", "MEME", "BUSTED", "aBSREL", "RELAX", "GARD"}
    assert s["note"].str.contains("not run").all() and len(read_tsv(d / "sites.tsv")) == 0
    write(d / "catalog.tsv", "locus_id\tgenome_id\tcatalytic_status\nG1__G1_1_1\tG1\tintact\n")
    run(selection_cmd(d, methods="FEL"))
    assert read_tsv(d / "sum.tsv")["note"].str.contains("only 1 analysable").any()
    write(d / "catalog.tsv", "locus_id\tgenome_id\n")
    must_fail(selection_cmd(d, methods="FEL"), "lacks column 'catalytic_status'")
    selection_inputs(d, n_codons=12)
    must_fail(selection_cmd(d, methods="FEL"), "hyphy not found",
              env_over={"PATH": path_without("hyphy", "HYPHYMP")})
    recs = read_fa(d / "codon.fasta")
    recs[0] = (recs[0][0], "TAA" + recs[0][1][3:])
    fasta(d / "codon.fasta", recs)
    must_fail(selection_cmd(d, methods="FEL"), "stop codon")


# --------------------------------------------------------------------------- #
# domain_trees.py
# --------------------------------------------------------------------------- #
@bg_job("domain_trees_real")
def _job_domain_trees_real():
    d = tmpdir("evo_dom_")
    prot, dom, ids = fixture_family(d)
    run(f"mafft --amino --retree 2 --thread {T} {prot} > {d / 'fam.aln'}")
    run(["iqtree2", "-s", str(d / "fam.aln"), "-m", "LG", "-fast", "-T", T, "--seed", "3",
         "--prefix", str(d / "fam"), "-redo", "-quiet"])
    cmd = [PYTHON, str(SCRIPTS_PY / "domain_trees.py"), "--proteins", str(prot), "--domtbl", str(dom),
           "--family-tree", str(d / "fam.treefile"), "--domains", DJSON, "--threads", T, "--seed", "5",
           "--workdir", str(d / "wd"), "--out", str(d / "rf.tsv"), "--model", "LG+G4", "--ufboot", "0"]
    return d, cmd, env()


@suite.test
def test_domain_trees_rf_conflict_known_answer():
    import domain_trees as dt
    d = tmpdir("evo_rf_")
    fam = write(d / "fam.nwk", "((A,B)100,(C,D)100,(E,F)100);")
    same = write(d / "same.nwk", "((A,B)99,(C,D)99,(E,F)99);")
    nni = write(d / "nni.nwk", "((A,C)100,(B,D)100,(E,F)100);")
    weak = write(d / "weak.nwk", "((A,C)50,(B,D)50,(E,F)50);")
    sh = set("ABCDEF")
    assert dt.compare(same, fam, sh) == (6, 0, 0.0, 3, 0)
    n, rf, norm, nspl, conf = dt.compare(nni, fam, sh)
    assert (rf, round(norm, 4), conf) == (4, 0.6667, 2), (rf, norm, conf)
    assert dt.compare(weak, fam, sh)[4] == 0, "unsupported splits are not conflicts"
    assert dt.compare(nni, fam, set("ABEF"))[1] == 0, "on {A,B,E,F} both trees agree"


@suite.test
def test_domain_span_blocks():
    import domain_trees as dt
    import utils
    spec = utils.load_domain_spec(DJSON)
    h = lambda dom, acc, s, e: {"domain": dom, "accession": acc, "start": s, "end": e, "env_from": s, "env_to": e}  # noqa: E731
    split = [h("AAA_11", "PF13086.1", 100, 200), h("AAA_12", "PF13087.1", 230, 400)]
    assert dt.domain_span(split, spec, "DNA2") == (100, 400, False)
    rear = [h("AAA_11", "PF13086.1", 100, 200), h("DUF4011", "PF13195.1", 250, 300), h("AAA_12", "PF13087.1", 350, 600)]
    assert dt.domain_span(rear, spec, "DNA2") == (350, 600, True)
    assert dt.domain_span(rear, spec, "DUF3320") is None


@suite.test
def test_domain_trees_insufficient_and_tool_failure():
    d = tmpdir("evo_domf_")
    prot, dom, ids = fixture_family(d)
    write(d / "triv.nwk", "();\n")
    py("domain_trees.py", "--proteins", prot, "--domtbl", dom, "--family-tree", d / "triv.nwk",
       "--domains", DJSON, "--workdir", d / "wd", "--out", d / "rf.tsv")
    rf = read_tsv(d / "rf.tsv")
    assert len(rf) == 4 and rf["note"].str.contains("not computed").all() and (rf["rf_to_family"] == "").all()
    write(d / "fam.nwk", "(" + ",".join(ids.values()) + ");")
    py("domain_trees.py", "--proteins", prot, "--domtbl", dom, "--family-tree", d / "fam.nwk",
       "--domains", DJSON, "--workdir", d / "wd2", "--out", d / "rf2.tsv", "--min-taxa", "50")
    assert read_tsv(d / "rf2.tsv")["note"].str.contains(r"\(< 50\)").all()
    must_fail([PYTHON, str(SCRIPTS_PY / "domain_trees.py"), "--proteins", prot, "--domtbl", dom,
               "--family-tree", str(d / "fam.nwk"), "--domains", DJSON, "--workdir", str(d / "wd3"),
               "--out", str(d / "rf3.tsv")], "not on PATH", env_over={"PATH": path_without("iqtree2", "iqtree")})


# --------------------------------------------------------------------------- #
# run_reconciliation.py
# --------------------------------------------------------------------------- #
REC_SPECIES = "((((G1:0.1,G2:0.1)90/95:0.1,G3:0.2)100/100:0.1,(G4:0.15,G5:0.15)100/100:0.15):0.1,G6:0.4);"
REC_GENES = ("((((G1__G1_1_1:0.1,G2__G2_1_1:0.1):0.1,G3__G3_1_1:0.2):0.1,(G4__G4_1_1:0.15,"
             "G5__G5_1_1:0.15):0.15):0.1,(G6__G6_1_1:0.2,G2__G2_2_9:0.2):0.2,G9__G9_1_1:0.3);")

MOCK_GENERAX = r'''#!/usr/bin/env python3
"""Mock of the GeneRax 2.1.3 I/O contract (Scenario::saveEventsCounts / savePerSpeciesEventsCounts /
saveTransfers formats); validates the families file and Treerecs mapping like Families.cpp."""
import os, sys
a = sys.argv[1:]
arg = lambda k: a[a.index(k) + 1]
mode = os.environ.get("MOCK_GENERAX_MODE", "ok")
if mode == "crash":
    sys.exit(139)
fam = {}
for line in open(arg("--families")):
    line = "".join(line.split())
    if "=" in line:
        k, v = line.split("=", 1); fam[k] = v
sp = open(arg("--species-tree")).read()
for k in ("starting_gene_tree", "alignment", "mapping", "subst_model"):
    assert k in fam, k
genes = [l.split() for l in open(fam["mapping"]) if l.strip()]
assert all(len(g) == 2 for g in genes), "mapping must be 'gene species'"
assert "sp_n1" in sp and "/" not in sp, "internal species nodes must carry unique labels"
assert arg("--rec-model") == "UndatedDTL" and arg("--seed") == "777" and arg("--strategy") == "SPR"
p = arg("--prefix"); os.makedirs(os.path.join(p, "reconciliations"), exist_ok=True)
print("[00:00:00] Filtering invalid families...")
if mode == "discard":
    print("Error in family prometheus: Failed to map genes and species")
    print("WARNING!!! Found 1 invalid families (they will be discarded from the analysis)")
    sys.exit(0)
leaf = 0 if mode == "empty" else len(genes)
with open(os.path.join(p, "reconciliations", "prometheus_eventCounts.txt"), "w") as fh:
    fh.write(f"S:10\nSL:2\nD:1\nT:3\nTL:1\nL:4\nLeaf:{leaf}\n")
with open(os.path.join(p, "reconciliations", "prometheus_speciesEventCounts.txt"), "w") as fh:
    fh.write("species_label, speciations, duplications, losses, transfers, presence, origination\n")
    fh.write("sp_n2, 2, 0, 0, 1, 1, 1\nG2, 1, 0, 0, 0, 1, 0\n")
with open(os.path.join(p, "reconciliations", "prometheus_transfers.txt"), "w") as fh:
    fh.write("G6 G2 1\nG6 G2 1\nG4 sp_n3 1\nG1 sp_n5 1\n")
print("[00:00:00] Reconciliation likelihood: -12.3")
'''
MOCK_ALEOBS = r'''#!/usr/bin/env python3
import sys, shutil
shutil.copy(sys.argv[1], sys.argv[1] + ".ale")
'''
MOCK_ALEML = r'''#!/usr/bin/env python3
"""ALEml_undated output contract (src/ALEml_undated.cpp): '<species basename>_<ale basename>.uml_rec'
in the cwd; a gene maps to the species named by the text before the first separator character."""
import os, re, sys
sp, ale = sys.argv[1], sys.argv[2]
opts = dict(x.split("=", 1) for x in sys.argv[3:])
assert opts["separators"] == "." and opts["seed"] == "777" and opts["sample"] == "100"
species = set(re.findall(r"[(,]([^(),:;]+)", open(sp).read()))
for g in re.findall(r"[(,]([^(),:;]+)", open(ale).read()):
    assert g.split(".")[0] in species, g
out = os.path.basename(sp) + "_" + os.path.basename(ale) + ".uml_rec"
open(out, "w").write("#ALEml_undated\n# of\t Duplications\tTransfers\tLosses\tSpeciations\n"
                     "Total \t0.5\t2\t3\t6\n")
'''


def recon_inputs(d):
    import dendropy
    write(d / "species.nwk", REC_SPECIES + "\n")
    write(d / "gene.nwk", REC_GENES + "\n")
    t = dendropy.Tree.get(data=REC_GENES, schema="newick", preserve_underscores=True)
    fasta(d / "aln.faa", [(lf.taxon.label, "MKVLAEGHKRL") for lf in t.leaf_node_iter()])
    write(d / "ufboot.nwk", "\n".join([REC_GENES] * 3) + "\n")
    write(d / "master.tsv", "genome_id\tmethanogen\n" + "".join(
        f"G{i}\t{'True' if i in (1, 2, 3) else 'False'}\n" for i in range(1, 7)))


def mock_bin(d, name, code):
    b = Path(d) / "mockbin"
    b.mkdir(exist_ok=True)
    f = b / name
    f.write_text(code.replace("#!/usr/bin/env python3", f"#!{PYTHON}", 1))
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return str(b)


def recon_cmd(d, extra=()):
    return [PYTHON, str(SCRIPTS_PY / "run_reconciliation.py"), "--aln", str(d / "aln.faa"),
            "--gene-tree", str(d / "gene.nwk"), "--ufboot", str(d / "ufboot.nwk"),
            "--species-tree", str(d / "species.nwk"), "--master", str(d / "master.tsv"),
            "--seed", "777", "--workdir", str(d / "wd"), "--out", str(d / "events.tsv"),
            "--out-transfers", str(d / "transfers.tsv"), *extra]


def with_path(b, **kw):
    e = env()
    e["PATH"] = b + os.pathsep + e["PATH"]
    e.update(kw)
    return e


@suite.test
def test_reconciliation_mock_generax_known_answer():
    d = tmpdir("evo_rec_")
    recon_inputs(d)
    b = mock_bin(d, "generax", MOCK_GENERAX)
    p = subprocess.run(recon_cmd(d), cwd=REPO, env=with_path(b), text=True, capture_output=True)
    assert p.returncode == 0, p.stderr[-2000:]
    assert "1 gene(s) map to genomes absent from the species tree" in p.stderr      # G9
    ev = read_tsv(d / "events.tsv").set_index("event")
    assert (ev.loc["Speciation", "count"], ev.loc["Transfer", "count"]) == ("12", "4")
    assert (ev.loc["Loss", "count"], ev.loc["Duplication", "count"]) == ("7", "1")
    assert float(ev.loc["Speciation", "fraction"]) == 0.75 == float(ev.loc["Transfer", "fraction"])
    assert "transfer_fraction=T/(S+T)=0.25" in ev.loc["Transfer", "note"]
    assert "origination species node=sp_n2" in ev.loc["Origination", "note"]
    tr = read_tsv(d / "transfers.tsv")
    got = {(r["donor"], r["recipient"]): (r["count"], r["recipient_methanogen_clade"]) for _, r in tr.iterrows()}
    # preorder labels: sp_n3 = ((G1,G2),G3) (all methanogens), sp_n5 = (G4,G5) (none)
    assert got[("G6", "G2")] == ("2.0", "True") and got[("G4", "sp_n3")] == ("1.0", "True"), got
    assert got[("G1", "sp_n5")] == ("1.0", "False"), got
    gi = d / "wd" / "generax_input"
    mapping = (gi / "mapping.link").read_text().split("\n")
    assert "G9__G9_1_1 G9" not in mapping and "G2__G2_2_9 G2" in mapping
    assert "G9__G9_1_1" not in (gi / "gene_tree.nwk").read_text()
    assert "G9__G9_1_1" not in (gi / "gene_alignment.faa").read_text()


@suite.test
def test_reconciliation_failures_are_loud():
    d = tmpdir("evo_recf_")
    recon_inputs(d)
    b = mock_bin(d, "generax", MOCK_GENERAX)
    for mode, pat in (("empty", "inconsistent: Leaf=0"), ("crash", "exited 139"), ("discard", "discarded the family")):
        must_fail(recon_cmd(d), pat, env_over={"PATH": with_path(b)["PATH"], "MOCK_GENERAX_MODE": mode})
    must_fail(recon_cmd(d), "generax not found", env_over={"PATH": path_without("generax")})
    write(d / "species_unrooted.nwk", "(G1,G2,(G3,(G4,(G5,G6))));\n")
    cmd = recon_cmd(d)
    cmd[cmd.index(str(d / "species.nwk"))] = str(d / "species_unrooted.nwk")
    must_fail(cmd, "not rooted", env_over={"PATH": with_path(b)["PATH"]})
    must_fail(recon_cmd(d, ("--model", "Auto")), "invalid GeneRax", env_over={"PATH": with_path(b)["PATH"]})
    must_fail(recon_cmd(d, ("--run-generax", "False", "--run-ale", "True")), "ALEobserve")


@suite.test
def test_reconciliation_disabled_and_trivial():
    d = tmpdir("evo_recd_")
    recon_inputs(d)
    run(recon_cmd(d, ("--run-generax", "False")))
    ev = read_tsv(d / "events.tsv")
    assert len(ev) == 5 and ev["note"].str.contains("disabled").all() and (ev["count"] == "0").all()
    write(d / "gene.nwk", "();\n")
    run(recon_cmd(d))
    assert read_tsv(d / "events.tsv")["note"].str.contains("too few tips").all()


@suite.test
def test_reconciliation_parsers_and_species_labels():
    import run_reconciliation as rr
    d = tmpdir("evo_recp_")
    _, tips = rr.species_tree_for_tools(write(d / "s.nwk", REC_SPECIES), str(d / "o.nwk"))
    out = (d / "o.nwk").read_text()
    assert sorted(re.findall(r"\)([^:;,()]+)", out)) == sorted(f"sp_n{i}" for i in range(1, 6)), out
    assert "90/95" not in out and tips == {f"G{i}" for i in range(1, 7)}
    rows = rr.event_rows("t", D=1, T=4, L=7, S=12, note="n")
    assert [r[2] for r in rows] == [1, 4, 7, 12, 1] and rows[1][3] == rows[3][3] == 0.75
    uml = write(d / "x.uml_rec", "#ALEml_undated using ALE v1.0\n\nS:\t(a,b);\n\n>logl: -10\n"
                "rate of\t Duplications\tTransfers\tLosses\nML \t0.1\t0.2\t0.3\n2 reconciled G-s:\n\n(a);\n(a);\n"
                "# of\t Duplications\tTransfers\tLosses\tSpeciations\nTotal \t1.5\t3.25\t10\t44.5\n\n"
                "# of\t Duplications\tTransfers\tLosses\tOriginations\tcopies\tsingletons\textinction_prob\tpresence\tLL\n")
    assert rr.parse_uml_rec(uml) == (1.5, 3.25, 10.0, 44.5)


@suite.test
def test_reconciliation_ale_mock_contract():
    d = tmpdir("evo_ale_")
    recon_inputs(d)
    b = mock_bin(d, "ALEobserve", MOCK_ALEOBS)
    mock_bin(d, "ALEml_undated", MOCK_ALEML)
    p = subprocess.run(recon_cmd(d, ("--run-generax", "False", "--run-ale", "True")), cwd=REPO,
                       env=with_path(b), text=True, capture_output=True)
    assert p.returncode == 0, p.stderr[-2000:]
    ev = read_tsv(d / "events.tsv")
    ev = ev[ev["tool"] == "ale"].set_index("event")
    assert (ev.loc["Transfer", "count"], ev.loc["Speciation", "count"]) == ("2", "6")
    assert abs(float(ev.loc["Speciation", "fraction"]) - 0.75) < 1e-9
    m = read_tsv(d / "wd" / "ale" / "gene_name_map.tsv")
    assert "G9__G9_1_1" not in set(m["gene_id"]) and m["ale_name"].str.match(r"^S\d+\.\d+$").all()


@suite.test
def test_generax_real_binary():
    """Real GeneRax 2.1.3 (bioconda, MPI build). The sandbox forbids socket bind(), so the binary is
    started as an isolated MPI singleton; it then reads our families / Treerecs mapping / relabelled
    species tree, but crashes inside the reconciliation (it does so on GeneRax's own example data
    too). We assert that (a) GeneRax accepted every input we wrote and (b) run_reconciliation.py
    FAILS instead of writing zero counts. Set PROM_GENERAX_WORKS=1 with a working binary on PATH to
    check the full run instead."""
    gx = shutil.which("generax", path=str(GENERAX_DIR)) or shutil.which("generax", path=env()["PATH"])
    assert gx, "generax not installed (create env prom-evolution-generax with generax=2.1.3)"
    import dendropy
    d = tmpdir("evo_gx_")
    recon_inputs(d)
    t = dendropy.Tree.get(data=REC_GENES, schema="newick", preserve_underscores=True)
    rng = random.Random(5)
    base = "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(120))
    fasta(d / "aln.faa", [(lf.taxon.label, "".join(c if rng.random() > 0.15 else rng.choice("ACDEFGHIKLMNPQRSTVWY")
                                                   for c in base)) for lf in t.leaf_node_iter()])
    (d / "tmp").mkdir()
    e = with_path(str(Path(gx).parent), TMPDIR=str(d / "tmp"), OMPI_MCA_ess_singleton_isolated="1",
                  OMPI_MCA_btl="self")
    p = subprocess.run(recon_cmd(d, ("--strategy", "EVAL")), cwd=REPO, env=e, text=True, capture_output=True)
    logf = d / "wd" / "generax.stdout.log"
    log = logf.read_text() if logf.exists() else p.stderr
    assert "Number of gene families: 1" in log and "Total number of genes: 7" in log, log[-1500:]
    assert "Error in family" not in log and "invalid families (they will" not in log, log[-1500:]
    if os.environ.get("PROM_GENERAX_WORKS") == "1":
        assert p.returncode == 0, p.stderr[-2000:]
        ev = read_tsv(d / "events.tsv").set_index("event")
        assert int(float(ev.loc["Speciation", "count"])) > 0
    else:
        assert p.returncode != 0 and re.search(r"exited|inconsistent", p.stderr), p.stderr[-1500:]


# --------------------------------------------------------------------------- #
# ancestral_states.R
# --------------------------------------------------------------------------- #
ANC_SETUP = r'''
suppressMessages(library(ape))
tr <- stree(64, "balanced"); tr$tip.label <- sprintf("g%02d", 1:64); tr <- compute.brlen(tr, 1)
write.tree(tr, "species.nwk")
tier <- rep("", 64); meth <- rep("False", 64)
for (k in c(0, 16, 32, 48)) { tier[k + 1:4] <- "A"; meth[k + 1:4] <- "True" }
tier[9] <- "B"            # B-only genome: EXCLUDED under {A}, PRESENT under {A,B}
d <- data.frame(genome_id = tr$tip.label, tree_tip = tr$tip.label, prometheus_tier = tier,
                methanogen = meth, prometheus_present = ifelse(tier != "", "True", "False"))
write.table(d, "master.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
'''


@suite.test
def test_ancestral_known_answer_tier_rule_and_enrichment():
    d = tmpdir("evo_anc_")
    write(d / "setup.R", ANC_SETUP)
    run([RSCRIPT, "setup.R"], cwd=d)
    rscript("ancestral_states.R", "--tree", d / "species.nwk", "--master", d / "master.tsv",
            "--nsim", "100", "--nnull", "400", "--min-accept", "100", "--seed", "4", "--tiers", "A",
            "--tiers-sensitivity", "A,B", "--out", d / "gl.tsv", "--out-null", d / "null.tsv",
            "--plot-prefix", d / "plot")
    gl = read_tsv(d / "gl.tsv")
    v = {(r["tier_set"], r["metric"]): r["value"] for _, r in gl.iterrows()}
    # contract v2 item 6 regression: the B-only genome g09
    assert v[("A", "n_tips_excluded_tier")] == "1" and v[("A", "n_tips_present")] == "16"
    assert v[("A,B", "n_tips_excluded_tier")] == "0" and v[("A,B", "n_tips_present")] == "17"
    assert v[("A", "n_tips_analysed")] == "63" and v[("A,B", "n_tips_analysed")] == "64"
    acq = float(v[("A", "n_independent_acquisitions")])
    assert 3.5 <= acq <= 4.6, acq                       # four carrier clades = four acquisitions
    assert float(v[("A", "n_gains")]) >= acq - 1e-9
    assert float(v[("A", "gains_on_methanogen_branches")]) >= 3.0
    p = float(v[("A", "enrichment_p")])
    assert p < 0.05, p
    assert float(v[("A", "rate_gain")]) > 0, "rates must come from the simmap $Q"
    note = gl.loc[(gl["tier_set"] == "A") & (gl["metric"] == "enrichment_p"), "note"].iloc[0]
    assert "prevalence-matched" in note, note
    nl = read_tsv(d / "null.tsv")
    assert set(nl["kind"]) == {"sim", "observed"} and set(nl["tier_set"]) == {"A", "A,B"}
    sims = nl[(nl.tier_set == "A") & (nl.kind == "sim")]["value"].astype(float)
    k = (sims >= float(v[("A", "gains_on_methanogen_branches")]) - 1e-9).sum()
    assert abs(p - (k + 1) / (len(sims) + 1)) < 1e-6, "p must be (k+1)/(n+1)"
    assert (d / "plot.png").stat().st_size > 1000 and (d / "plot.svg").exists()


ANC_UNIT = r'''
source("ancestral_states.R", chdir = TRUE)
stopifnot(identical(tier_presence(c("A", "B", "", "C"), "A"), c(1L, NA, 0L, NA)))
stopifnot(identical(tier_presence(c("A", "B", "", "C"), c("A", "B")), c(1L, 1L, 0L, NA)))
r <- try(tier_presence("D", "A"), silent = TRUE); stopifnot(inherits(r, "try-error"))
tr <- read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);")
dv <- descendant_all(tr, matrix(c(TRUE, TRUE, FALSE, TRUE), ncol = 1))
stopifnot(identical(as.vector(dv[, 1]), c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)))
set.seed(1)
s0 <- sim_binary(tr, 0, 0, 1, 50); stopifnot(all(s0))                 # no change from root = 1
s1 <- sim_binary(compute.brlen(tr, 1000), 1, 3, 0, 4000)
stopifnot(abs(mean(s1) - 0.25) < 0.03)                                 # stationary q01/(q01+q10)
# handmade simmap: root 0; gain on the edge to (a,b); on the edge to b a loss and a regain
tr2 <- reorder(tr, "cladewise")
maps <- lapply(seq_len(nrow(tr2$edge)), function(i) c("0" = 1))
ab <- which(tr2$edge[, 2] == 6); ea <- which(tr2$edge[, 2] == 1); eb <- which(tr2$edge[, 2] == 2)
maps[[ab]] <- c("0" = 0.5, "1" = 0.5); maps[[ea]] <- c("1" = 1)
maps[[eb]] <- c("1" = 0.3, "0" = 0.3, "1" = 0.4)
tr2$maps <- maps
stopifnot(acquisitions_one(tr2) == 2)       # the (a,b) gain seeds a; the regain seeds b
g <- edge_gain_means(list(tr2)); stopifnot(sum(g) == 2, g[ab] == 1, g[eb] == 1)
cat("ANC UNIT OK\n")
'''


@suite.test
def test_ancestral_internal_functions():
    d = tmpdir("evo_ancu_")
    for f in ("ancestral_states.R", "plotting_theme.R"):
        shutil.copy(SCRIPTS_R / f, d / f)
    write(d / "unit.R", ANC_UNIT)
    p = run([RSCRIPT, "unit.R"], cwd=d)
    assert "ANC UNIT OK" in p.stdout, p.stdout + p.stderr


@suite.test
def test_ancestral_empty_and_failures():
    d = tmpdir("evo_ancf_")
    write(d / "setup.R", ANC_SETUP)
    run([RSCRIPT, "setup.R"], cwd=d)
    write(d / "triv.nwk", "();\n")
    rscript("ancestral_states.R", "--tree", d / "triv.nwk", "--master", d / "master.tsv",
            "--out", d / "gl.tsv", "--out-null", d / "null.tsv", "--plot-prefix", d / "p")
    gl = read_tsv(d / "gl.tsv")
    assert (gl["value"] == "NA").all() and gl["note"].str.contains("trivial").all()
    assert (d / "p.png").exists() and len(read_tsv(d / "null.tsv")) == 0
    write(d / "unrooted.nwk", "(g01:1,g02:1,(g03:1,(g04:1,g05:1):1):1);\n")
    must_fail([RSCRIPT, str(SCRIPTS_R / "ancestral_states.R"), "--tree", str(d / "unrooted.nwk"),
               "--master", str(d / "master.tsv"), "--out", str(d / "x.tsv"), "--plot-prefix", str(d / "x")],
              "not rooted")
    write(d / "bad_master.tsv", "genome_id\tmethanogen\ng01\tTrue\n")
    must_fail([RSCRIPT, str(SCRIPTS_R / "ancestral_states.R"), "--tree", str(d / "species.nwk"),
               "--master", str(d / "bad_master.tsv"), "--out", str(d / "x.tsv"), "--plot-prefix", str(d / "x")],
              "prometheus_tier")


# --------------------------------------------------------------------------- #
# plot_tanglegram.R and the Python plots
# --------------------------------------------------------------------------- #
@suite.test
def test_tanglegram_real_and_edge_cases():
    d = tmpdir("evo_tang_")
    write(d / "sp.nwk", REC_SPECIES + "\n")
    write(d / "gene.nwk", REC_GENES + "\n")
    write(d / "m.tsv", "genome_id\tmethanogen\tis_anme\n" + "".join(
        f"G{i}\t{'True' if i < 3 else 'False'}\t{'True' if i == 5 else 'False'}\n" for i in range(1, 7)))
    p = rscript("plot_tanglegram.R", "--gene-tree", d / "gene.nwk", "--species-tree", d / "sp.nwk",
                "--master", d / "m.tsv", "--plot-prefix", d / "tg")
    assert "7 links, 6 carrier genomes, 1 unmapped genes" in p.stderr, p.stderr[-800:]
    assert (d / "tg.png").stat().st_size > 5000 and (d / "tg.svg").stat().st_size > 1000
    write(d / "triv.nwk", "();\n")
    rscript("plot_tanglegram.R", "--gene-tree", d / "triv.nwk", "--species-tree", d / "sp.nwk",
            "--plot-prefix", d / "ph")
    assert (d / "ph.png").exists()
    write(d / "broken.nwk", "((G1__x,G2__y;\n")
    must_fail([RSCRIPT, str(SCRIPTS_R / "plot_tanglegram.R"), "--gene-tree", str(d / "broken.nwk"),
               "--species-tree", str(d / "sp.nwk"), "--plot-prefix", str(d / "bk")])


SITES = ("site\talignment_col\tdomain\tbeta_alpha_or_dnds\tp_value\tclassification\tfel_alpha\tfel_beta\t"
         "meme_p_value\tmeme_classification\tn_seqs_with_codon\tdomain_support\n"
         + "".join(f"{i}\t{i}\t{'DUF4011' if i <= 10 else ('NA' if i <= 12 else 'REase')}\t"
                   f"{'inf' if i == 3 else ('' if i == 4 else 0.5 + i / 10)}\t{0.01 if i in (2, 15, 16) else 0.5}\t"
                   f"{'positive' if i in (2, 15, 16) else 'neutral'}\t1\t1\t{0.02 if i == 15 else 0.6}\t"
                   f"{'episodic' if i == 15 else 'neutral'}\t5\t1.0\n" for i in range(1, 21)))


@suite.test
def test_plot_selection_known_counts_and_placeholder():
    import plot_selection as ps
    import pandas as pd
    d = tmpdir("evo_ps_")
    write(d / "sites.tsv", SITES)
    sites = pd.read_csv(d / "sites.tsv", sep="\t", keep_default_na=False, na_values=[""])
    tab = ps.per_domain(sites).set_index("domain")
    assert list(tab.index) == ["DUF4011", "REase", "NA"]
    assert (tab.loc["DUF4011", "n"], tab.loc["DUF4011", "k_fel"]) == (10, 1)
    assert (tab.loc["REase", "n"], tab.loc["REase", "k_fel"], tab.loc["REase", "k_meme"]) == (8, 2, 1)
    lo, hi = ps.wilson(5, 10)
    assert abs(lo - 0.2366) < 1e-3 and abs(hi - 0.7634) < 1e-3
    write(d / "sum.tsv", "method\tstatistic\tvalue\tp_value\tn_sites_or_branches\tnote\n"
          "FEL\tsites_positive\t3\t\t20\tp<=0.1 and beta>alpha; n=5\n")
    py("plot_selection.py", "--summary", d / "sum.tsv", "--sites", d / "sites.tsv", "--outdir", d / "o")
    assert (d / "o" / "selection_by_domain.png").stat().st_size > 1000
    write(d / "empty.tsv", SITES.split("\n")[0] + "\n")
    py("plot_selection.py", "--summary", d / "sum.tsv", "--sites", d / "empty.tsv", "--outdir", d / "o2")
    assert (d / "o2" / "selection_by_domain.svg").exists()


@suite.test
def test_plot_evolution_extra_inputs_and_tier_filter():
    import plot_evolution_extra as pe
    import pandas as pd
    d = tmpdir("evo_pe_")
    write(d / "sites.tsv", SITES)
    write(d / "rec.tsv", "tool\tevent\tcount\tfraction\tnote\n"
          "generax\tDuplication\t1\t\tx\ngenerax\tTransfer\t4\t0.75\tx\ngenerax\tLoss\t7\t\tx\n"
          "generax\tSpeciation\t12\t0.75\tx\ngenerax\tOrigination\t1\t\tx\n"
          "ale\tTransfer\t0\t\ttoo few genes on the species tree (2)\n")
    write(d / "gl.tsv", "metric\tvalue\tnote\ttier_set\n"
          "gains_on_methanogen_branches\t4\tx\tA\nexpected_gains_methanogen\t1.2\tx\tA\n"
          "enrichment_p\t0.01\tprevalence-matched\tA\ngains_on_methanogen_branches\t9\tx\tA,B\n")
    write(d / "null.tsv", "kind\tvalue\ttier_set\n" + "".join(f"sim\t{i % 3}\tA\n" for i in range(50))
          + "observed\t4\tA\nsim\t99\tA,B\n")
    gl = pe._primary(pd.read_csv(d / "gl.tsv", sep="\t"))
    assert set(gl["tier_set"]) == {"A"} and len(gl) == 3
    nl = pe._primary(pe._load_null(str(d / "null.tsv")))
    assert 99 not in set(nl["value"].astype(float))
    py("plot_evolution_extra.py", "--sites", d / "sites.tsv", "--reconciliation", d / "rec.tsv",
       "--gains", d / "gl.tsv", "--gains-null", d / "null.tsv", "--outdir", d / "o")
    for f in ("dnds_along_protein", "reconciliation_events", "gain_enrichment_null"):
        assert (d / "o" / f"{f}.png").stat().st_size > 1000 and (d / "o" / f"{f}.svg").exists(), f
    for f in ("e1", "e2", "e3", "e4"):
        write(d / f"{f}.tsv", "")
    py("plot_evolution_extra.py", "--sites", d / "e1.tsv", "--reconciliation", d / "e2.tsv",
       "--gains", d / "e3.tsv", "--gains-null", d / "e4.tsv", "--outdir", d / "o2")
    assert len(list((d / "o2").glob("*.png"))) == 3


# The two tests below wait for the background jobs started by prefetch(); they are defined
# (and therefore run) last so that the jobs overlap with all the other tests.
@suite.test
def test_selection_real_hyphy_pinned_all_methods():
    assert have("hyphy"), "hyphy not on PATH"
    d, p = background("selection_all")
    e = hyphy_env()
    assert p.returncode == 0, p.stderr[-3000:]
    ver = subprocess.run(["hyphy", "--version"], env=e, text=True, capture_output=True).stdout.strip()
    s = read_tsv(d / "sum.tsv")
    rows = {(m, st): r for m, st, r in zip(s["method"], s["statistic"], s.to_dict("records"))}
    for key in [("FEL", "sites_positive"), ("FEL", "sites_negative"), ("MEME", "sites_episodic"),
                ("BUSTED", "episodic_diversifying_p"), ("aBSREL", "branches_selected"), ("RELAX", "K"),
                ("GARD", "breakpoints")]:
        assert key in rows, f"{key} missing ({ver})"
    fel = rows[("FEL", "sites_positive")]
    assert fel["n_sites_or_branches"] == "30" and "intact-only: 6/7" in fel["note"], fel
    assert 0 <= float(rows[("BUSTED", "episodic_diversifying_p")]["p_value"]) <= 1
    assert rows[("aBSREL", "branches_selected")]["n_sites_or_branches"] == "9", "2n-3 branches, n=6"
    relax = rows[("RELAX", "K")]
    assert 0 <= float(relax["p_value"]) <= 1 and float(relax["value"]) >= 0
    # clade labelling on the unrooted 6-tip tree (G7 pruned as non-intact): Test = G1, G2, G3,
    # (G1,G2) and the {G1,G2,G3}|{G4,G5,G6} branch; Reference = G4, G5, (G4,G5); G6 has an ANME
    # host -> unclassified. 5 + 3 + 1 = 2n - 3 branches.
    assert "Test=5 methanogen branches (3 terminal), Reference=3, Unclassified=1" in relax["note"], relax["note"]
    lab = (d / "wd" / "relax.labelled.nwk").read_text()
    assert "G6__G6_1_1:" in lab and "G6__G6_1_1{" not in lab and "G7" not in lab
    sites = read_tsv(d / "sites.tsv")
    assert list(sites["alignment_col"].astype(int)) == list(range(1, 31))
    dom = dict(zip(sites["alignment_col"].astype(int), sites["domain"]))
    assert all(dom[c] == "DUF4011" for c in range(1, 16)) and all(dom[c] == "REase" for c in range(16, 31)), dom
    cov = dict(zip(sites["alignment_col"].astype(int), sites["n_seqs_with_codon"].astype(int)))
    assert cov[25] == 5 and cov[24] == 6, "codon 25 is a gap in G5; G7 excluded as non-intact"
    assert set(sites["classification"]) <= {"positive", "negative", "neutral"}
    assert set(sites["meme_classification"]) <= {"episodic", "neutral"}


@suite.test
def test_domain_trees_real_tools_on_fixture():
    import utils
    d, p = background("domain_trees_real")
    assert p.returncode == 0, p.stderr[-3000:]
    dom = d / "dom.tbl"                            # written by fixture_family() in the job
    rf = read_tsv(d / "rf.tsv")
    assert list(rf["domain"]) == list(DOMAINS), list(rf["domain"])
    spec, hits = utils.load_domain_spec(DJSON), utils.parse_domtbl(dom)
    for _, r in rf.iterrows():
        exp = sum(1 for hs in hits.values() if any(utils.hit_category(h, spec) == r["domain"] for h in hs))
        assert int(r["n_seqs"]) == exp, (r["domain"], r["n_seqs"], exp)
        assert r["note"].startswith("ok"), r["note"]
        n = int(r["n_tips"])
        assert 0 <= int(r["rf_to_family"]) <= 2 * (n - 3)
        assert abs(float(r["normalized_rf"]) - int(r["rf_to_family"]) / (2 * (n - 3))) < 1e-6
        assert 0 <= float(r["au_pvalue_family_topology"]) <= 1 and r["model"] == "LG+G4"
        assert r["n_supported_conflicts"] == "", "no UFBoot -> conflicts not computable"


if __name__ == "__main__":
    prefetch()
    suite.main("EVOLUTION TESTS PASSED")
