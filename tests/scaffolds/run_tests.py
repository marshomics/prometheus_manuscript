#!/usr/bin/env python3
"""Synthetic test suite for Stage 3 (scaffolds).

Covers root_tree.py, mcr_search.py, classify_methanogen.py, build_master_table.py and the
shell of the scaffolds.smk rules, on inputs with KNOWN expected answers, with real
IQ-TREE 2 / HMMER / MAFFT / R-ape / Snakemake runs. Run from the repo root:

    bash tests/run_leaf.sh scaffolds        # prints SCAFFOLDS TESTS PASSED on success

Every tool used here is required; a missing tool fails the suite (no silent skips).
At most 2 threads per tool.
"""
from __future__ import annotations

import contextlib
import gzip
import io
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True          # importing the owned scripts must not litter __pycache__
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, run, tmpdir, have, read_tsv, REPO, PYTHON,  # noqa: E402
                  SCRIPTS_PY, RSCRIPT, env, TOOLS_BIN)

sys.path.insert(0, str(SCRIPTS_PY))
FIX = HERE / "fixtures"
THREADS = "2"
OWNED_PY = ["root_tree.py", "mcr_search.py", "classify_methanogen.py", "build_master_table.py"]
suite = Suite("scaffolds")


# --------------------------------------------------------------------------- #
# environment: a TMPDIR that subprocesses (MAFFT's mktemp) can write
# --------------------------------------------------------------------------- #
def _writable_tmp():
    for cand in ("/tmp", os.environ.get("TMPDIR"), str(HERE / ".tmp")):
        if not cand:
            continue
        cand = os.path.abspath(cand)
        try:
            os.makedirs(cand, exist_ok=True)
            r = subprocess.run(["mktemp", "-d", os.path.join(cand, "promtest.XXXXXX")],
                               capture_output=True, text=True)
            if r.returncode == 0:
                shutil.rmtree(r.stdout.strip(), ignore_errors=True)
                return cand
        except OSError:
            continue
    return None


_TMP = _writable_tmp()
if _TMP:
    os.environ["TMPDIR"] = _TMP
    tempfile.tempdir = _TMP


def need(*tools):
    for t in tools:
        assert have(t), f"required tool {t!r} not on PATH (PROM_TOOLS_BIN={TOOLS_BIN!r})"


def write(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def tsv(path, header, rows):
    lines = ["\t".join(header)] + ["\t".join(str(x) for x in r) for r in rows]
    return write(path, "\n".join(lines) + "\n")


def fasta(path, recs):
    return write(path, "".join(f">{n}\n{s}\n" for n, s in recs))


def read_fa(path):
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0]
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


def fails(proc, *needles):
    """Assert a CompletedProcess failed and its output mentions every needle."""
    out = (proc.stdout or "") + (proc.stderr or "")
    assert proc.returncode != 0, f"expected failure, got exit 0:\n{out[-1500:]}"
    for n in needles:
        assert n in out, f"expected {n!r} in failure output:\n{out[-1500:]}"


def alisim(tree, length, seed, d, name="sim"):
    """Simulate a protein alignment with IQ-TREE AliSim on a known tree."""
    need("iqtree2")
    t = write(Path(d) / f"{name}.nwk", tree + "\n")
    run(["iqtree2", "--alisim", str(Path(d) / name), "-m", "LG+G", "-t", str(t),
         "--length", str(length), "--seed", str(seed), "-af", "fasta", "-redo"])
    return read_fa(Path(d) / f"{name}.fa")


def dtree(text_or_path):
    import dendropy
    s = str(text_or_path).strip()
    raw = s if s.startswith("(") else Path(s).read_text().strip()
    return dendropy.Tree.get(data=raw, schema="newick", preserve_underscores=True)


def leafset(nd):
    return frozenset(lf.taxon.label for lf in nd.leaf_iter())


def split_labels(tree):
    """{orientation-free split: label} for labelled internal non-root nodes."""
    tips = frozenset(lf.taxon.label for lf in tree.leaf_node_iter())
    anchor = min(tips)
    out = {}
    for nd in tree.preorder_internal_node_iter():
        if nd is tree.seed_node or not nd.label:
            continue
        s = leafset(nd)
        key = s if anchor not in s else tips - s
        out[key] = nd.label
    return out


def root_children(tree):
    return sorted((leafset(c) for c in tree.seed_node.child_nodes()), key=lambda s: sorted(s))


def patristic(tree):
    """{(a, b): distance} for every tip pair (independent of rooting)."""
    pdm = tree.phylogenetic_distance_matrix()
    taxa = sorted(tree.taxon_namespace, key=lambda t: t.label)
    return {(a.label, b.label): pdm.patristic_distance(a, b)
            for i, a in enumerate(taxa) for b in taxa[i + 1:]}


def root_tree(*args, check=False):
    return py("root_tree.py", *args, check=check)


# --------------------------------------------------------------------------- #
# static checks
# --------------------------------------------------------------------------- #
@suite.test
def test_py_compile_owned_scripts():
    import py_compile
    for s in OWNED_PY:
        py_compile.compile(str(SCRIPTS_PY / s), doraise=True)


@suite.test
def test_rules_hygiene():
    p = run([PYTHON, str(REPO / "tests" / "check_rules.py"),
             str(REPO / "workflow" / "rules" / "scaffolds.smk")])
    assert "RULES HYGIENE OK" in p.stdout, p.stdout


# --------------------------------------------------------------------------- #
# root_tree.py
# --------------------------------------------------------------------------- #
PMSF_TRUE = ("((((G1:0.1,G2:0.12):0.08,(G3:0.1,G4:0.15):0.07):0.1,(G5:0.2,G6:0.1):0.12):0.05,"
             "((OG1:0.3,OG2:0.25):0.3,XA:0.4):0.1);")


@suite.test
def test_root_outgroup_on_real_iqtree_pmsf_tree():
    """Real IQ-TREE 2 two-step PMSF (the rule's flags) -> outgroup rooting by EXACT GTDB
    taxon; root bifurcating; support labels stay on their bipartitions; ape reads a
    rooted binary tree with unquoted numeric supports."""
    need("iqtree2")
    d = tmpdir()
    seqs = alisim(PMSF_TRUE, 220, 7, d)
    msa = fasta(d / "msa.faa", sorted(seqs.items()))
    run(["iqtree2", "-s", str(msa), "-m", "LG+F+G", "-T", THREADS, "--seed", "12345",
         "--prefix", str(d / "guide"), "-redo"])
    p = run(["iqtree2", "-s", str(msa), "-m", "LG+C20+F+G", "-ft", str(d / "guide.treefile"),
             "-B", "1000", "-alrt", "1000", "-T", THREADS, "--seed", "12345",
             "--prefix", str(d / "species"), "-redo"])
    assert (d / "species.sitefreq").exists(), "PMSF did not write the site-frequency file"
    assert "COMPUTING SITE FREQUENCY MODEL" in p.stdout, p.stdout[-800:]
    tax = tsv(d / "tax.tsv", ["user_genome", "classification"],
              [["OG1", "d__Archaea;p__Altiarchaeota;c__X"], ["OG2", "d__Archaea;p__Altiarchaeota;c__X"],
               ["XA", "d__Archaea;p__Altiarchaeota_A;c__Z"]] +
              [[f"G{i}", "d__Archaea;p__Methanobacteriota;c__M"] for i in range(1, 7)])
    keep = write(d / "keep.txt", "genome_id\n" + "\n".join(sorted(seqs)) + "\n")
    out, log = d / "rooted.nwk", d / "rooting.txt"
    root_tree("--in-tree", d / "species.treefile", "--taxonomy", tax, "--method", "outgroup",
              "--outgroup-taxon", "p__Altiarchaeota", "--keep-tips", keep, "--msa", msa,
              "--out-membership", d / "mem.tsv", "--out-tree", out, "--log", log, check=True)
    text = out.read_text().strip()
    assert "[&R]" not in text and "'" not in text, text
    t = dtree(text)
    assert len(t.seed_node.child_nodes()) == 2
    assert root_children(t) == sorted([frozenset({"OG1", "OG2"}),
                                       frozenset(seqs) - {"OG1", "OG2"}], key=lambda s: sorted(s))
    before, after = split_labels(dtree(d / "species.treefile")), split_labels(t)
    assert before and all(after.get(k) == v for k, v in before.items()), (before, after)
    assert log.read_text().splitlines()[0].startswith("method=outgroup(2 tips;monophyletic)")
    assert "outgroup_tips=OG1,OG2" in log.read_text()          # XA (p__Altiarchaeota_A) excluded
    r = run([RSCRIPT, "-e",
             "suppressMessages(library(ape)); t <- read.tree('%s'); "
             "s <- suppressWarnings(as.numeric(sub('/.*', '', t$node.label[t$node.label != ''])));"
             "cat('ROOTED', is.rooted(t), 'BINARY', is.binary(t), 'NTIP', Ntip(t), "
             "'NUMERIC', all(!is.na(s)), '\\n')" % out])
    assert "ROOTED TRUE BINARY TRUE NTIP 9 NUMERIC TRUE" in r.stdout, r.stdout + r.stderr
    mem = read_tsv(d / "mem.tsv")
    assert set(mem["status"]) == {"in_tree"} and len(mem) == 9


@suite.test
def test_support_labels_follow_bipartitions_after_rerooting():
    """Regression: DendroPy 4.6.1 rerooting leaves node labels on nodes, which moves
    supports to the wrong splits; root_tree re-keys them by bipartition."""
    import dendropy
    nwk = "((A:1,B:1)s1:1,(C:1,D:1)s2:1,(E:1,(F:1,G:1)s3:1)s4:1);"
    d = tmpdir()
    write(d / "in.nwk", nwk)
    root_tree("--in-tree", d / "in.nwk", "--out-tree", d / "out.nwk", "--method", "outgroup",
              "--outgroup", "F,G", "--log", d / "log.txt", check=True)
    before, after = split_labels(dtree(nwk)), split_labels(dtree(d / "out.nwk"))
    for k, v in before.items():
        assert after.get(k) == v, f"split {sorted(k)} lost/moved its label {v}: {after}"
    naive = dendropy.Tree.get(data=nwk, schema="newick")
    naive.is_rooted = True
    naive.reroot_at_edge(naive.mrca(taxon_labels=["F", "G"]).edge, update_bipartitions=True)
    moved = [k for k, v in split_labels(naive).items() if before.get(k) != v]
    print(f"    (naive dendropy rerooting mislabels {len(moved)} split(s); root_tree: 0)")


def _random_tree(rng, n):
    nodes = [f"T{i}" for i in range(n)]
    while len(nodes) > 3:
        a, b = rng.sample(range(len(nodes)), 2)
        x, y = nodes[a], nodes[b]
        new = f"({x}:{rng.uniform(0.01, 1):.4f},{y}:{rng.uniform(0.01, 1):.4f})"
        nodes = [v for i, v in enumerate(nodes) if i not in (a, b)] + [new]
    return "(" + ",".join(f"{v}:{rng.uniform(0.01, 1):.4f}" for v in nodes) + ");"


@suite.test
def test_midpoint_rooting_is_exact_and_preserves_distances():
    """Midpoint root: bifurcating, every patristic distance unchanged, and the deepest
    tip on each side of the root at exactly half the diameter (incl. midpoint-on-node
    cases where DendroPy's own reroot_at_midpoint loses a tip or leaves a polytomy)."""
    rng = random.Random(42)
    trees = [_random_tree(rng, rng.randint(4, 30)) for _ in range(25)]
    trees += ["(A:1,B:1,C:1);", "((A:1,B:1):1,(C:1,D:1):1,E:0.5);"]
    d = tmpdir()
    import root_tree as rt
    for i, nwk in enumerate(trees):
        write(d / "in.nwk", nwk)
        args = ["--in-tree", str(d / "in.nwk"), "--out-tree", str(d / "out.nwk"),
                "--method", "midpoint"]
        if i < 2 or i >= len(trees) - 2:          # CLI for a few, in-process for the rest
            root_tree(*args, check=True)
        else:
            with contextlib.redirect_stderr(io.StringIO()):
                assert rt.main(args) == 0
        t_in, t_out = dtree(nwk), dtree(d / "out.nwk")
        assert len(t_out.seed_node.child_nodes()) == 2, (i, nwk)
        p_in, p_out = patristic(t_in), patristic(t_out)
        assert p_in.keys() == p_out.keys(), (i, "tip set changed")
        assert all(abs(p_in[k] - p_out[k]) < 1e-6 for k in p_in), (i, nwk)
        diam = max(p_in.values())
        depth = {}
        for lf in t_out.leaf_node_iter():
            top = lf
            while top.parent_node is not t_out.seed_node:
                top = top.parent_node
            depth[id(top)] = max(depth.get(id(top), 0.0), lf.distance_from_root())
        assert all(abs(v - diam / 2) < 1e-6 for v in depth.values()), (i, nwk, depth, diam)
    assert (d / "out.nwk").read_text().count("C") == 1


@suite.test
def test_outgroup_failures_and_explicit_fallback():
    nwk = "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);"
    d = tmpdir()
    write(d / "in.nwk", nwk)
    base = ["--in-tree", d / "in.nwk", "--out-tree", d / "o.nwk", "--log", d / "l.txt"]
    p = root_tree(*base, "--outgroup", "A,C")
    fails(p, "outgroup_not_monophyletic")
    log = (d / "l.txt").read_text()
    assert log.startswith("method=failed") and "outgroup_intruders=" in log, log
    p = root_tree(*base, "--outgroup", "A,C", "--on-outgroup-failure", "midpoint", check=True)
    assert (d / "l.txt").read_text().startswith("method=midpoint(fallback:outgroup_not_monophyletic)")
    fails(root_tree(*base, "--outgroup", "Z"), "outgroup_not_in_tree")
    fails(root_tree(*base, "--outgroup", "A,B,C,D,E,F"), "outgroup_is_all_tips")
    fails(root_tree(*base, "--method", "bogus"), "unknown --method")
    root_tree(*base, check=True)
    assert (d / "l.txt").read_text().startswith("method=midpoint(no_outgroup_configured)")
    root_tree(*base, "--outgroup", "E,F", check=True)
    t = dtree(d / "o.nwk")
    assert frozenset({"E", "F"}) in root_children(t)
    # already-rooted input is unrooted first, then re-rooted on the requested edge
    write(d / "in.nwk", "((A:1,B:1):1,((C:1,D:1):1,(E:1,F:1):1):1);")
    root_tree(*base, "--outgroup", "C,D", check=True)
    assert frozenset({"C", "D"}) in root_children(dtree(d / "o.nwk"))


@suite.test
def test_outgroup_taxon_is_exact_rank_token_match():
    import root_tree as rt
    tax = {"a": "d__Archaea;p__Altiarchaeota;c__X", "b": "d__Archaea;p__Altiarchaeota_A;c__Y",
           "c": "d__Archaea;p__Thermoproteota;c__Altiarchaeota"}
    assert rt.resolve_taxa("p__Altiarchaeota", [tax])[0] == {"a"}
    assert rt.resolve_taxa("p__altiarchaeota", [tax])[0] == {"a"}
    assert rt.resolve_taxa("Altiarchaeota", [tax])[0] == {"a", "c"}        # bare name: any rank
    assert rt.resolve_taxa("c__Y, p__Altiarchaeota", [tax])[0] == {"a", "b"}
    assert rt.resolve_taxa("p__Altiarch", [tax])[0] == set()               # no substring match
    assert rt.iqtree_name("Gen+1:x'y|a/b-c.d_e") == "Gen_1_x_y|a/b-c.d_e"


@suite.test
def test_trivial_tree_and_membership_reasons():
    d = tmpdir()
    write(d / "in.nwk", "();\n")
    keep = write(d / "keep.txt", "genome_id\nA\nB\nC\nD\n")
    msa = fasta(d / "msa.faa", [("A", "MKV"), ("B", "MKI")])
    dropped = write(d / "dropped.tsv", "C\tar53 --min_perc_aa filter: Insufficient (4.2%)\n")
    root_tree("--in-tree", d / "in.nwk", "--out-tree", d / "o.nwk", "--keep-tips", keep,
              "--msa", msa, "--gtdbtk-dropped", dropped, "--out-membership", d / "m.tsv",
              "--log", d / "l.txt", check=True)
    assert (d / "o.nwk").read_text().strip() == "();"
    m = {r["genome_id"]: r for _, r in read_tsv(d / "m.tsv").iterrows()}
    assert m["A"]["status"].startswith("msa_too_small_for_tree") and m["A"]["in_msa"] == "True"
    assert m["C"]["status"] == "dropped_by_gtdbtk: ar53 --min_perc_aa filter: Insufficient (4.2%)"
    assert m["D"]["status"] == "not_in_gtdbtk_msa: no reason reported"
    assert (d / "l.txt").read_text().startswith("method=none(trivial_tree")


@suite.test
def test_iqtree_renamed_tips_are_restored_to_genome_ids():
    """IQ-TREE 2 rewrites '+' (and other characters) to '_' in sequence names; the
    canonical tree must carry the genome_ids again. An id that would need Newick
    quoting is refused loudly."""
    need("iqtree2")
    d = tmpdir()
    seqs = alisim("((G1:0.1,G2:0.12):0.1,(G3:0.1,G4:0.15):0.07,(G5:0.2,G6:0.1):0.12);", 300, 3, d)
    ids = {"G1": "Gen+1", "G2": "MAG.2|b", "G3": "x-3", "G4": "G4", "G5": "G5", "G6": "G6"}
    msa = fasta(d / "msa.faa", [(ids[k], v) for k, v in sorted(seqs.items())])
    run(["iqtree2", "-s", str(msa), "-m", "LG+G", "-fast", "-T", "1", "--seed", "1",
         "--prefix", str(d / "t"), "-redo"])
    assert "Gen_1" in (d / "t.treefile").read_text()          # IQ-TREE renamed it
    keep = write(d / "keep.txt", "\n".join(ids.values()) + "\n")
    root_tree("--in-tree", d / "t.treefile", "--out-tree", d / "o.nwk", "--keep-tips", keep,
              "--log", d / "l.txt", check=True)
    tips = {lf.taxon.label for lf in dtree(d / "o.nwk").leaf_node_iter()}
    assert tips == set(ids.values()), tips
    assert "relabelled_iqtree_tips=1" in (d / "l.txt").read_text()
    ids["G1"] = "MAG:1"
    msa = fasta(d / "msa2.faa", [(ids[k], v) for k, v in sorted(seqs.items())])
    run(["iqtree2", "-s", str(msa), "-m", "LG+G", "-fast", "-T", "1", "--seed", "1",
         "--prefix", str(d / "t2"), "-redo"])
    keep = write(d / "keep2.txt", "\n".join(ids.values()) + "\n")
    fails(root_tree("--in-tree", d / "t2.treefile", "--out-tree", d / "o2.nwk",
                    "--keep-tips", keep), "unquoted")


@suite.test
def test_domain_scope_reference_outgroup_then_prune():
    """tree.scope=domain: root with a GTDB-reference outgroup (resolved through the
    headerless reference taxonomy), then prune the references; the analysed-genome
    tree keeps the root position the outgroup implied."""
    nwk = ("(((U1:0.1,U2:0.1):0.2,(U3:0.1,RS_GCF_1.1:0.1):0.2):0.3,"
           "(RS_GCF_2.1:0.2,RS_GCF_3.1:0.2):0.4,U4:0.5);")
    d = tmpdir()
    write(d / "in.nwk", nwk)
    reft = write(d / "ref.tsv", "RS_GCF_1.1\td__Archaea;p__In;c__I\n"
                                "RS_GCF_2.1\td__Archaea;p__Out;c__O\nRS_GCF_3.1\td__Archaea;p__Out;c__O\n")
    keep = write(d / "keep.txt", "U1\nU2\nU3\nU4\n")
    msa = fasta(d / "msa.faa", [(u, "MK") for u in ("U1", "U2", "U3", "U4")])
    root_tree("--in-tree", d / "in.nwk", "--ref-taxonomy", reft, "--outgroup-taxon", "p__Out",
              "--keep-tips", keep, "--msa", msa, "--out-membership", d / "m.tsv",
              "--out-full-tree", d / "full.nwk", "--out-tree", d / "o.nwk", "--log", d / "l.txt",
              check=True)
    full, t = dtree(d / "full.nwk"), dtree(d / "o.nwk")
    assert frozenset({"RS_GCF_2.1", "RS_GCF_3.1"}) in root_children(full)
    assert {lf.taxon.label for lf in t.leaf_node_iter()} == {"U1", "U2", "U3", "U4"}
    assert root_children(t) == [frozenset({"U1", "U2", "U3"}), frozenset({"U4"})]
    assert "n_reference_tips_pruned=3" in (d / "l.txt").read_text()
    assert set(read_tsv(d / "m.tsv")["status"]) == {"in_tree"}


@suite.test
def test_unknown_tip_and_missing_dendropy_are_errors():
    d = tmpdir()
    write(d / "in.nwk", "((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);")
    keep = write(d / "keep.txt", "A\nB\nC\nD\nE\n")
    fails(root_tree("--in-tree", d / "in.nwk", "--out-tree", d / "o.nwk", "--keep-tips", keep),
          "neither an analysed genome")
    # python -S hides site-packages: dendropy cannot import -> must fail, not pass through
    p = run([PYTHON, "-S", str(SCRIPTS_PY / "root_tree.py"), "--in-tree", str(d / "in.nwk"),
             "--out-tree", str(d / "o3.nwk"), "--method", "midpoint"], check=False)
    fails(p, "dendropy")


# --------------------------------------------------------------------------- #
# mcr_search.py (real HMMER 3.4)
# --------------------------------------------------------------------------- #
def run_env(cmd, extra_path=None, check=False):
    """Run with the toolchain PATH, optionally with `extra_path` FIRST (for fake tools)."""
    e = env()
    if extra_path:
        e["PATH"] = str(extra_path) + os.pathsep + e["PATH"]
    p = subprocess.run([str(c) for c in cmd], cwd=REPO, env=e, text=True, capture_output=True)
    if check and p.returncode != 0:
        raise AssertionError(f"command failed ({p.returncode}): {cmd}\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
    return p


def fake_tool(d, name, code=3, msg="simulated crash"):
    b = Path(d) / "fakebin"
    b.mkdir(parents=True, exist_ok=True)
    write(b / name, f"#!/bin/bash\necho '{msg}' >&2\nexit {code}\n").chmod(0o755)
    return b


def build_hmm(fixture, name, d, ga=None):
    need("mafft", "hmmbuild")
    d = Path(d)
    aln = d / f"{name}.afa"
    with open(aln, "w") as fh:
        p = subprocess.run(["mafft", "--quiet", "--thread", THREADS, str(fixture)], stdout=fh,
                           stderr=subprocess.PIPE, text=True, env=env())
    assert p.returncode == 0, p.stderr
    hmm = d / f"{name}.hmm"
    run(["hmmbuild", "-n", name, "--amino", str(hmm), str(aln)])
    if ga:   # HMMER3 header order: GA follows CKSUM (elsewhere hmmsearch 3.4 segfaults)
        out = []
        for ln in hmm.read_text().splitlines(keepends=True):
            out.append(ln)
            if ln.startswith("CKSUM"):
                out.append(f"GA    {ga[0]:.2f} {ga[1]:.2f};\n")
        hmm.write_text("".join(out))
    return hmm


def mcr_world(d):
    """Proteome / pmap / genes.tsv for six genomes with known mcr content."""
    A = list(read_fa(FIX / "uniprot_mcrA.faa").values())
    B = list(read_fa(FIX / "uniprot_mcrB.faa").values())
    G = list(read_fa(FIX / "uniprot_mcrG.faa").values())
    rng = random.Random(3)
    shuf = list(A[0])
    rng.shuffle(shuf)
    frag = A[2][:int(0.6 * len(A[2]))]
    genomes = {
        "G1": [(A[0], "00"), (B[0], "00"), (G[0], "00")],      # complete Mcr
        "G2": [(frag, "10"), (B[1], "00")],                     # McrA truncated at a contig end
        "G3": [(frag, "00")],                                   # same length, complete gene
        "G4": [("".join(shuf), "00"), ("MSTNPKPQRKTKRNTNRRPQDVKFPGG" * 10, "00")],  # decoys
        "G5": [(A[3], "00"), (A[7], "00")],                     # two McrA copies
        "G6": [(A[0][:250] + "GSGSGSGS" + A[0][:250], "00")],   # tandem: full >> best domain
    }
    recs, pmap = [], []
    gdir = Path(d) / "genes"
    for g, prots in genomes.items():
        rows = []
        for i, (seq, flag) in enumerate(prots, 1):
            local = f"{g}_1_{i}"
            recs.append((f"{g}__{local}", seq))
            pmap.append([f"{g}__{local}", g, local])
            rows.append([local, "c1", 1, i, 100 * i, 100 * i + 3 * len(seq), "+", flag, 90000, len(seq)])
        tsv(gdir / g / f"{g}.genes.tsv", ["protein_id", "contig", "contig_index", "gene_index",
                                          "start", "end", "strand", "partial", "contig_length",
                                          "length_aa"], rows)
    faa = fasta(Path(d) / "proteome.faa", recs)
    pm = tsv(Path(d) / "pmap.tsv", ["protein_id", "genome_id", "orig_id"], pmap)
    return faa, pm, gdir


def independent_scores(hmm, faa, d):
    """hmmsearch -T 0 tblout parsed by Biopython SearchIO: {target: (full, best-domain)}."""
    from Bio import SearchIO
    tbl = Path(d) / f"{Path(hmm).stem}.indep.tbl"
    run(["hmmsearch", "--cpu", THREADS, "-T", "0", "-o", os.devnull, "--tblout", str(tbl),
         str(hmm), str(faa)])
    out = {}
    for q in SearchIO.parse(str(tbl), "hmmer3-tab"):
        for h in q.hits:
            out[h.id] = (h.bitscore, h.hsps[0].bitscore)
    return out


def ko_list(path, rows):
    hdr = ["knum", "threshold", "score_type", "profile_type", "F-measure", "nseq", "nseq_used",
           "alen", "mlen", "eff_nseq", "re/pos", "definition"]
    return tsv(path, hdr, [[k, t, st, "all", "1.0", 1, 1, 1, 1, 1, 1, "x"] for k, t, st in rows])


def mcr(*args, check=False, extra_path=None):
    return run_env([PYTHON, SCRIPTS_PY / "mcr_search.py", *args], extra_path, check)


@suite.test
def test_mcr_search_real_hmmsearch_thresholds_fragments_copies():
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    ha = build_hmm(FIX / "uniprot_mcrA.faa", "K00399", hd)
    hb = build_hmm(FIX / "uniprot_mcrB.faa", "K00401", hd)
    hg = build_hmm(FIX / "uniprot_mcrG.faa", "K00402", hd)
    faa, pm, gdir = mcr_world(d)
    sa, sb, sg = (independent_scores(h, faa, d) for h in (ha, hb, hg))
    t_a = round(0.8 * sa["G1__G1_1_1"][0], 2)
    t_b = round(0.8 * min(sb["G1__G1_1_2"][0], sb["G2__G2_1_2"][0]), 2)
    t_g = round(0.8 * sg["G1__G1_1_3"][0], 2)
    frag = sa["G2__G2_1_1"][0]
    assert 0.5 * t_a <= frag < t_a, ("fixture premise: fragment between 0.5T and T", frag, t_a)
    assert sa["G6__G6_1_1"][0] >= t_a > sa["G6__G6_1_1"][1], ("tandem premise", sa["G6__G6_1_1"], t_a)
    kl = ko_list(d / "ko_list", [("K00399", t_a, "full"), ("K00401", t_b, "full"),
                                 ("K00402", t_g, "full")])
    out, allh = d / "mcr_hits.tsv", d / "mcr_all.tsv"
    mcr("--proteome", faa, "--pmap", pm, "--hmmdir", hd, "--score-fraction", "1.0",
        "--partial-score-fraction", "0.5", "--genes-dir", gdir, "--ko-list", kl,
        "--kos", "K00399,K00401,K00402", "--threads", THREADS, "--out", out, "--out-all", allh,
        check=True)
    h = read_tsv(out)
    assert list(h.columns[:5]) == ["genome_id", "gene", "protein_id", "score", "evalue"]
    calls = {(r.genome_id, r.gene): r for r in h.itertuples()}
    assert set(calls) == {("G1", "mcrA"), ("G1", "mcrB"), ("G1", "mcrG"), ("G2", "mcrA"),
                          ("G2", "mcrB"), ("G5", "mcrA"), ("G6", "mcrA")}, sorted(calls)
    assert calls[("G1", "mcrA")].call == "full_threshold"
    g2 = calls[("G2", "mcrA")]
    assert (g2.call, g2.partial) == ("partial_rescue", "10")            # contig-edge fragment
    g5 = calls[("G5", "mcrA")]
    best5 = max(("G5__G5_1_1", "G5__G5_1_2"), key=lambda p: sa[p][0])
    assert g5.n_copies == "2" and g5.protein_id == best5
    for r in h.itertuples():                               # scores == independent parse
        s = {"mcrA": sa, "mcrB": sb, "mcrG": sg}[r.gene][r.protein_id]
        assert abs(float(r.score) - s[0]) < 0.051 and abs(float(r.best_domain_score) - s[1]) < 0.051
        assert r.threshold_source == "ko_list" and r.score_type == "full"
    dec = {(r.protein_id, r.gene): r.decision for r in read_tsv(allh).itertuples()}
    assert dec[("G3__G3_1_1", "mcrA")] == "below_threshold"           # complete gene: no rescue
    assert ("G4__G4_1_1", "mcrA") not in dec or dec[("G4__G4_1_1", "mcrA")] == "below_threshold"


@suite.test
def test_mcr_search_kofam_domain_score_type_uses_best_domain():
    """KofamScan: score_type 'domain' compares the best-1-domain score (tblout col 9)."""
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    ha = build_hmm(FIX / "uniprot_mcrA.faa", "K00399", hd)
    faa, pm, gdir = mcr_world(d)
    sa = independent_scores(ha, faa, d)
    t = round(0.8 * sa["G1__G1_1_1"][0], 2)
    for st, expect_tandem in (("full", True), ("domain", False)):
        kl = ko_list(d / f"kl_{st}", [("K00399", t, st)])
        out = d / f"hits_{st}.tsv"
        mcr("--proteome", faa, "--pmap", pm, "--hmmdir", hd, "--ko-list", kl, "--kos", "K00399",
            "--threads", THREADS, "--out", out, check=True)
        genomes = set(read_tsv(out)["genome_id"])
        assert ("G6" in genomes) == expect_tandem, (st, genomes)
        assert "G1" in genomes


@suite.test
def test_mcr_search_ga_mode_matches_hmmer_cut_ga():
    """GA-threshold profiles (TIGRFAM naming): at fraction 1 the accepted set equals the
    sequences hmmsearch --cut_ga reports WITH a domain (domtblout). HMMER 3.4 lists a
    sequence passing GA1 with no domain >= GA2 in --tblout (N=0): the tandem protein."""
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    probe = build_hmm(FIX / "uniprot_mcrA.faa", "probe", d)
    faa, pm, gdir = mcr_world(d)
    sa = independent_scores(probe, faa, d)
    ga = round(0.8 * sa["G1__G1_1_1"][0], 2)
    h = build_hmm(FIX / "uniprot_mcrA.faa", "TIGR03256", hd, ga=(ga, ga))
    allh = d / "all.tsv"
    mcr("--proteome", faa, "--pmap", pm, "--hmmdir", hd, "--kos", "K00399", "--threads", THREADS,
        "--out", d / "hits.tsv", "--out-all", allh, check=True)
    mine = {r.protein_id for r in read_tsv(allh).itertuples() if r.decision == "accepted_full"}
    run(["hmmsearch", "--cut_ga", "--cpu", THREADS, "-o", os.devnull, "--tblout", str(d / "c.tbl"),
         "--domtblout", str(d / "c.dtbl"), str(h), str(faa)])
    with_domain = {ln.split()[0] for ln in (d / "c.dtbl").read_text().splitlines()
                   if ln.strip() and not ln.startswith("#")}
    tbl = {ln.split()[0] for ln in (d / "c.tbl").read_text().splitlines()
           if ln.strip() and not ln.startswith("#")}
    assert mine == with_domain, (mine, with_domain)
    assert "G6__G6_1_1" in tbl and "G6__G6_1_1" not in mine
    assert read_tsv(d / "hits.tsv")["threshold_source"].eq("GA").all()


@suite.test
def test_mcr_search_half_kofam_threshold_is_permissive():
    """Evidence for kofam_score_fraction 1.0: against the real KOfam K00399 threshold
    (741.47, ko_list fixture) a complete 250-aa McrA piece (partial '00') passes at the
    old 0.5 fraction but not at KofamScan's own 1.0."""
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    ha = build_hmm(FIX / "uniprot_mcrA.faa", "K00399", hd)
    a0 = next(iter(read_fa(FIX / "uniprot_mcrA.faa").values()))
    faa = fasta(d / "p.faa", [("H__H_1_1", a0[:250]), ("F__F_1_1", a0)])
    pm = tsv(d / "pm.tsv", ["protein_id", "genome_id", "orig_id"],
             [["H__H_1_1", "H", "H_1_1"], ["F__F_1_1", "F", "F_1_1"]])
    s = independent_scores(ha, faa, d)
    assert 0.5 * 741.47 <= s["H__H_1_1"][0] < 741.47 <= s["F__F_1_1"][0], s
    kl = FIX / "kofam_ko_list_mcr_subset.tsv"
    got = {}
    for f in ("0.5", "1.0"):
        out = d / f"h{f}.tsv"
        mcr("--proteome", faa, "--pmap", pm, "--hmmdir", hd, "--ko-list", kl, "--kos", "K00399",
            "--score-fraction", f, "--out", out, check=True)
        got[f] = set(read_tsv(out)["genome_id"])
    assert got == {"0.5": {"H", "F"}, "1.0": {"F"}}, got


@suite.test
def test_mcr_search_ko_list_reader_on_real_rows():
    import mcr_search as ms
    kl = ms.parse_ko_list(str(FIX / "kofam_ko_list_mcr_subset.tsv"))
    assert kl["K00399"] == {"threshold": 741.47, "score_type": "full"}
    assert kl["K00401"]["threshold"] == 506.80 and kl["K00402"]["threshold"] == 359.37
    assert ms.map_model("K00399", "", "x") == ("mcrA", "K00399")
    assert ms.map_model("TIGR03257", "TIGR03257.1", "y")[0] == "mcrB"
    assert ms.map_model("NAME", "", "K00402")[0] == "mcrG"
    assert ms.map_model("PF00001", "", "random")[0] is None


@suite.test
def test_mcr_search_failures_are_loud_and_empty_input_is_header_only():
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    faa, pm, gdir = mcr_world(d)
    base = ["--proteome", faa, "--pmap", pm, "--out", d / "o.tsv"]
    fails(mcr(*base, "--hmmdir", hd), "no *.hmm")                          # no profiles
    build_hmm(FIX / "uniprot_mcrA.faa", "K00399", hd)
    fails(mcr(*base, "--hmmdir", hd), "no curated threshold")              # KOfam HMM, no ko_list
    fails(mcr(*base, "--hmmdir", hd, "--ko-list", d / "missing"), "does not exist")
    kl = ko_list(d / "kl", [("K00399", 400.0, "full")])
    fails(mcr(*base, "--hmmdir", hd, "--ko-list", kl, "--kos", "K00399,K00401"),
          "configured subunit")                                            # mcrB HMM absent
    fails(mcr(*base, "--hmmdir", hd, "--ko-list", kl, "--kos", "K00399",
              extra_path=fake_tool(d, "hmmsearch")), "hmmsearch failed")  # crashed tool
    p = tsv(d / "pm_short.tsv", ["protein_id", "genome_id", "orig_id"], [["G1__G1_1_1", "G1", "x"]])
    fails(mcr("--proteome", faa, "--pmap", p, "--out", d / "o.tsv", "--hmmdir", hd,
              "--ko-list", kl), "not in the protein map")
    shutil.rmtree(gdir / "G2")
    fails(mcr(*base, "--hmmdir", hd, "--ko-list", kl, "--partial-score-fraction", "0.2",
              "--genes-dir", gdir), "genes table missing")
    empty = write(d / "empty.faa", "")
    mcr("--proteome", empty, "--pmap", pm, "--hmmdir", hd, "--ko-list", kl,
        "--out", d / "e.tsv", "--out-all", d / "e_all.tsv", check=True)
    import mcr_search as ms
    assert (d / "e.tsv").read_text() == ms.HEADER + "\n"


# --------------------------------------------------------------------------- #
# classify_methanogen.py (real MAFFT --add --keeplength + IQ-TREE 2)
# --------------------------------------------------------------------------- #
MCRA_TREE = ("((((M1:0.1,(M2:0.08,cM:0.03):0.04):0.05,(M3:0.1,M4:0.1):0.05):0.35,"
             "((A1:0.06,cA:0.03):0.06,(A2:0.1,A3:0.12):0.05):0.35):0.25,"
             "((O1:0.08,cO:0.03):0.05,O2:0.1):0.6);")


def classify(*args, check=False, extra_path=None):
    return run_env([PYTHON, SCRIPTS_PY / "classify_methanogen.py", *args], extra_path, check)


def mcra_world(d, labels=None, drop_anme=False):
    seqs = alisim(MCRA_TREE, 450, 21, d, "mcra")
    refs = [k for k in sorted(seqs) if not k.startswith("c")]
    if drop_anme:
        refs = [r for r in refs if not r.startswith("A")]
    ref_aln = fasta(Path(d) / "ref.afa", [(r, seqs[r]) for r in refs])
    lab = labels or {"M": "methanogen", "A": "ANME-1", "O": "other"}
    ref_lab = tsv(Path(d) / "ref_labels.tsv", ["ref_id", "group"], [[r, lab[r[0]]] for r in refs])
    rng = random.Random(5)
    cA = list(seqs["cA"])
    for _ in range(6):                     # small deletions: MAFFT must realign them
        del cA[rng.randrange(len(cA))]
    prots = {"GM__GM_1_1": seqs["cM"], "GA__GA_1_1": "".join(cA) + "*",
             "GO__GO_1_1": seqs["cO"][:int(0.7 * len(seqs["cO"]))]}
    faa = fasta(Path(d) / "prot.faa", sorted(prots.items()))
    pm = tsv(Path(d) / "pmap.tsv", ["protein_id", "genome_id", "orig_id"],
             [[p, p.split("__")[0], p.split("__")[1]] for p in prots] + [["GX__GX_1_1", "GX", "GX_1_1"]])
    hits = tsv(Path(d) / "hits.tsv", ["genome_id", "gene", "protein_id", "score", "evalue", "call",
                                      "n_copies"],
               [["GM", "mcrA", "GM__GM_1_1", 900, 0, "full_threshold", 1],
                ["GM", "mcrB", "GM__GM_1_2", 700, 0, "full_threshold", 1],
                ["GA", "mcrA", "GA__GA_1_1", 850, 0, "full_threshold", 2],
                ["GO", "mcrA", "GO__GO_1_1", 500, 0, "partial_rescue", 1]])
    return ref_aln, ref_lab, faa, pm, hits


@suite.test
def test_classify_placement_with_real_mafft_iqtree():
    need("mafft", "iqtree2")
    d = tmpdir()
    ref_aln, ref_lab, faa, pm, hits = mcra_world(d)
    out, plc = d / "calls.tsv", d / "placement.tsv"
    classify("--hits", hits, "--proteome", faa, "--pmap", pm, "--ref-aln", ref_aln,
             "--ref-labels", ref_lab, "--required", "K00399", "--anme-tree", "True",
             "--threads", THREADS, "--seed", "12345", "--workdir", d / "wd", "--out", out,
             "--out-placement", plc, check=True)
    c = {r.genome_id: r for r in read_tsv(out).itertuples()}
    assert set(c) == {"GM", "GA", "GO", "GX"}
    exp = {"GM": ("methanogen", "True", "False", "unresolved"),
           "GA": ("ANME", "False", "True", "NA"),
           "GO": ("other", "False", "False", "NA"),
           "GX": ("NA", "False", "False", "NA")}
    for g, (grp, meth, anme, sub) in exp.items():
        r = c[g]
        assert (r.mcr_group, r.methanogen, r.is_anme, r.methanogenesis_subtype) == (grp, meth, anme, sub), (g, r)
    assert c["GM"].n_mcr == "2" and c["GX"].evidence == "no_mcr"
    assert "mcrA=contig_edge_fragment" in c["GO"].evidence and c["GA"].n_mcrA_copies == "2"
    assert c["GM"].nearest_ref in ("M2",), c["GM"].nearest_ref
    p = read_tsv(plc)
    assert len(p) == 3 and set(p["mcr_group"]) == {"methanogen", "ANME", "other"}
    assert (d / "wd" / "mcrA.treefile").exists()
    hdr = read_tsv(out).columns.tolist()
    assert hdr[:7] == ["genome_id", "mcrA", "n_mcr", "methanogen", "is_anme",
                       "methanogenesis_subtype", "evidence"]


@suite.test
def test_nearest_labelled_clade_rule_units():
    import classify_methanogen as cm
    groups = {"M1": "methanogen", "M2": "methanogen", "A1": "ANME", "A2": "ANME", "A3": "ANME",
              "A4": "ANME", "O1": "other"}
    # candidate sister to a mixed clade -> ambiguous
    t = dtree("((c:0.1,(M1:0.1,A1:0.1):0.1):0.1,(M2:2.0,A2:2.1):0.1,O1:0.2);")
    r = cm.classify_tree(t, groups, {"c"})["c"]
    assert r["mcr_group"] == "ambiguous" and r["clade_ref_groups"] == "ANME:1,methanogen:1", r
    # long branches: clade says ANME although the patristically nearest reference is M1
    t = dtree("((c:0.05,(A3:1.0,A4:1.0):0.05):0.05,(M1:0.02,M2:0.02):0.4,O1:3.0);")
    r = cm.classify_tree(t, groups, {"c"})["c"]
    assert r["mcr_group"] == "ANME" and r["nearest_ref"] == "M1", r
    # candidates nested among candidates: walk up to the first node holding references
    t = dtree("(((c1:0.05,c2:0.05):0.05,A1:0.1):0.05,(M1:0.1,M2:0.1):0.3,O1:2.0);")
    r = cm.classify_tree(t, groups, {"c1", "c2"})
    assert r["c1"]["mcr_group"] == r["c2"]["mcr_group"] == "ANME", r
    assert cm.norm_group("ANME-2a") == "ANME" and cm.norm_group("Methanogen") == "methanogen"
    for bad in ("methanotroph", "ACR", ""):
        try:
            cm.norm_group(bad)
        except ValueError:
            continue
        raise AssertionError(f"label {bad!r} accepted")
    assert cm.parse_required("K00399,mcrb") == ["mcrA", "mcrB"]


@suite.test
def test_classify_without_placement_and_required_subunits():
    d = tmpdir()
    hits = tsv(d / "hits.tsv", ["genome_id", "gene", "protein_id", "score", "evalue", "call"],
               [["H1", "mcrA", "H1__a", 900, 0, "full_threshold"],
                ["H2", "mcrA", "H2__a", 900, 0, "full_threshold"],
                ["H2", "mcrB", "H2__b", 700, 0, "full_threshold"],
                ["H3", "mcrB", "H3__b", 700, 0, "full_threshold"],
                ["H4", "mcrA", "H4__a", 500, 0, "partial_rescue"]])
    faa = write(d / "p.faa", "")
    pm = tsv(d / "pm.tsv", ["protein_id", "genome_id", "orig_id"], [["H5__x", "H5", "x"]])
    for req, meth in (("K00399,K00401", {"H2"}), ("K00399", {"H1", "H2", "H4"})):
        out = d / "c.tsv"
        classify("--hits", hits, "--proteome", faa, "--pmap", pm, "--required", req,
                 "--anme-tree", "False", "--workdir", d / "wd", "--out", out, check=True)
        c = {r.genome_id: r for r in read_tsv(out).itertuples()}
        assert {g for g, r in c.items() if r.methanogen == "True"} == meth, (req, c)
        assert all(r.is_anme == "False" for r in c.values())
        assert c["H5"].evidence == "no_mcr" and c["H3"].mcrA == "False"
        if req == "K00399,K00401":
            assert "required_missing=mcrB" in c["H1"].evidence
        else:
            assert c["H1"].mcr_group == "not_placed" and "anme_check_disabled" in c["H1"].evidence
            assert "mcrA=contig_edge_fragment" in c["H4"].evidence
            assert {g for g, r in c.items() if r.methanogenesis_subtype == "unresolved"} == meth


@suite.test
def test_classify_failures_are_loud_and_empty_input_is_header_only():
    need("mafft", "iqtree2")
    d = tmpdir()
    ref_aln, ref_lab, faa, pm, hits = mcra_world(d)
    base = ["--hits", hits, "--proteome", faa, "--pmap", pm, "--required", "K00399",
            "--anme-tree", "True", "--workdir", d / "wd", "--out", d / "o.tsv"]
    fails(classify(*base, "--ref-aln", d / "missing.afa", "--ref-labels", ref_lab), "missing")
    bad = tsv(d / "bad_labels.tsv", ["ref_id", "group"],
              [[r, "methanotroph"] for r in read_fa(ref_aln)])
    fails(classify(*base, "--ref-aln", ref_aln, "--ref-labels", bad), "not one of")
    d2 = tmpdir()
    ra2, rl2, _, _, _ = mcra_world(d2, drop_anme=True)
    fails(classify(*base, "--ref-aln", ra2, "--ref-labels", rl2), "both ANME")
    faa2 = fasta(d / "p2.faa", [("GM__GM_1_1", "MKV")])
    fails(classify("--hits", hits, "--proteome", faa2, "--pmap", pm, "--required", "K00399",
                   "--anme-tree", "True", "--ref-aln", ref_aln, "--ref-labels", ref_lab,
                   "--workdir", d / "wd", "--out", d / "o.tsv"), "not in the proteome")
    fails(classify(*base, "--ref-aln", ref_aln, "--ref-labels", ref_lab,
                   extra_path=fake_tool(d, "mafft")), "mafft --add failed")
    fails(classify(*base, "--ref-aln", ref_aln, "--ref-labels", ref_lab, "--required", "K99999"),
          "neither an mcr KO")
    e = d / "empty"
    hh = tsv(e / "h.tsv", ["genome_id", "gene", "protein_id", "score", "evalue"], [])
    pp = tsv(e / "p.tsv", ["protein_id", "genome_id", "orig_id"], [])
    classify("--hits", hh, "--proteome", write(e / "x.faa", ""), "--pmap", pp,
             "--required", "K00399", "--anme-tree", "False", "--workdir", e / "wd",
             "--out", e / "o.tsv", check=True)
    import classify_methanogen as cm
    assert (e / "o.tsv").read_text() == cm.HEADER + "\n"


# --------------------------------------------------------------------------- #
# build_master_table.py
# --------------------------------------------------------------------------- #
SCHEMA_BACKBONE = [
    "genome_id", "representative", "tree_tip", "gtdb_domain", "gtdb_phylum", "gtdb_class",
    "gtdb_order", "gtdb_family", "gtdb_genus", "gtdb_species", "completeness", "contamination",
    "quality_tier", "pass_primary", "genome_size", "mcrA", "n_mcr", "methanogen", "is_anme",
    "methanogenesis_subtype", "prometheus_present", "prometheus_tier", "prometheus_n_loci",
    "detect_defensefinder", "n_defense_systems", "n_defense_systems_excl_prometheus"]
DF_HDR = ["sys_id", "type", "subtype", "activity", "sys_beg", "sys_end", "protein_in_syst",
          "genes_count", "name_of_profiles_in_sys"]
GIDS = ["G1", "G2", "G3", "G4", "G5", "G6"]


def master_world(d, catalog_rows=None, tree="((G1:1,G2:1):1,(G3:1,G4:1):1,G5:1);",
                 eco_rows=None, meta_ids=None):
    d = Path(d)
    ids = meta_ids or GIDS
    meta = tsv(d / "meta.tsv", ["genome_id", "source", "completeness", "contamination",
                                "genome_size", "n_contigs", "n50", "gc", "quality_tier",
                                "pass_primary"],
               [[g, "src", 95 + i * 0.5, 1.0, 2000000 + i, 10 + i, 50000, 0.4, "A", True]
                for i, g in enumerate(ids)])
    tax = tsv(d / "tax.tsv", ["user_genome", "classification"],
              [[g, f"d__Archaea;p__P{i % 2};c__C;o__O;f__F;g__Gen;s__Gen sp{i}"] for i, g in enumerate(ids)])
    write(d / "tree.nwk", tree + "\n")
    mem = tsv(d / "mem.tsv", ["genome_id", "in_msa", "in_tree", "status"],
              [[g, g != "G6", g != "G6",
                "in_tree" if g != "G6" else "dropped_by_gtdbtk: ar53 --min_perc_aa filter: x"]
               for g in ids])
    meth = tsv(d / "meth.tsv", ["genome_id", "mcrA", "n_mcr", "methanogen", "is_anme",
                                "methanogenesis_subtype", "evidence", "mcr_group"],
               [["G1", True, 3, True, False, "unresolved", "e", "methanogen"],
                ["G2", True, 1, False, True, "NA", "e", "ANME"]])
    eco = tsv(d / "species.tsv", ["genome_id", "classification", "niche_Gut", "group_Marine",
                                  "group_Human"],
              eco_rows or [["G1", "x", 1, 0, 1], ["G2", "x", 0, 1, ""], ["G3", "x", 2, 0, 0],
                           ["G4", "x", 0, 0, 0], ["G5", "x", 1, 1, 1], ["G6", "x", "", "", ""]])
    cat = tsv(d / "catalog.tsv", ["locus_id", "genome_id", "protein_id", "confidence_tier"],
              catalog_rows or [["l1", "G1", "p1", "A"], ["l2", "G1", "p2", "B"],
                               ["l3", "G2", "p3", "B"], ["l4", "G3", "p4", "C"],
                               ["l5", "G4", "p5", "Rejected"], ["l6", "G6", "p6", "A"],
                               ["l7", "G6", "p7", "Rejected"]])
    sens = tsv(d / "sens.tsv", ["genome_id", "detect_defensefinder", "prometheus_present"],
               [[g, 1 if g in ("G1", "G2", "G3", "G4", "G6") else 0, ""] for g in ids])
    reps = write(d / "reps.txt", "\n".join(ids) + "\n")
    dfd = d / "df"
    for g in ids:
        rows = []
        if g == "G1":
            rows = [["s1", "Prometheus", "Prometheus", "Defense", 1, 2, "p1", 1, "Prometheus__ProA"],
                    ["s2", "RM", "RM_Type_I", "Defense", 3, 5, "q", 3, "x"],
                    ["s3", "PrometheusLike", "x", "Defense", 6, 7, "q", 1, "x"],
                    ["s4", "Anti_X", "Anti_X", "Antidefense", 8, 9, "q", 1, "x"]]
        elif g == "G2":
            rows = [["s1", "prometheus", "Prometheus", "Defense", 1, 2, "p3", 1, "x"]]
        tsv(dfd / g / f"{g}_defense_finder_systems.tsv", DF_HDR, rows)
    return dict(meta=meta, taxonomy=tax, tree=d / "tree.nwk", methanogen=meth, ecology=eco,
                catalog=cat, sensitivity=sens, reps=reps, df_dir=dfd, membership=mem)


def master(w, d, *extra, check=False):
    return run_env([PYTHON, SCRIPTS_PY / "build_master_table.py",
                    "--meta", w["meta"], "--taxonomy", w["taxonomy"], "--tree", w["tree"],
                    "--methanogen", w["methanogen"], "--ecology", w["ecology"],
                    "--envo-map", "config/envo_map.tsv", "--catalog", w["catalog"],
                    "--sensitivity", w["sensitivity"], "--reps", w["reps"], "--df-dir", w["df_dir"],
                    "--tree-membership", w["membership"], "--out-master", Path(d) / "master.tsv",
                    "--out-ecology", Path(d) / "eco.tsv", "--out-missingness", Path(d) / "miss.tsv",
                    *extra], None, check)


@suite.test
def test_master_table_values_and_presence_rule():
    d = tmpdir()
    w = master_world(d)
    p = master(w, d, "--primary-tiers", "A", check=True)
    m = read_tsv(d / "master.tsv")
    assert m.columns.tolist()[:len(SCHEMA_BACKBONE)] == SCHEMA_BACKBONE
    eco_block = m.columns.tolist()[len(SCHEMA_BACKBONE):len(SCHEMA_BACKBONE) + 3]
    assert eco_block == ["group_Human", "group_Marine", "niche_Gut"], eco_block
    r = {x.genome_id: x for x in m.itertuples()}
    exp = {  # tier, present, n_loci, status A, status A+B
        "G1": ("A", "True", "2", "present", "present"),
        "G2": ("B", "True", "1", "excluded", "present"),     # B-only: never a confident absence
        "G3": ("C", "True", "1", "excluded", "excluded"),
        "G4": ("", "False", "0", "absent", "absent"),       # Rejected only
        "G5": ("", "False", "0", "absent", "absent"),
        "G6": ("A", "True", "1", "present", "present")}     # A + Rejected
    for g, e in exp.items():
        got = (r[g].prometheus_tier, r[g].prometheus_present, r[g].prometheus_n_loci,
               r[g].prometheus_status_primary, r[g].prometheus_status_sensitivity)
        assert got == e, (g, got, e)
    assert (r["G1"].n_defense_systems, r["G1"].n_defense_systems_excl_prometheus) == ("3", "2")
    assert (r["G2"].n_defense_systems, r["G2"].n_defense_systems_excl_prometheus) == ("1", "0")
    assert r["G5"].n_defense_systems == "0" and r["G5"].detect_defensefinder == "0"
    assert r["G1"].tree_tip == "G1" and r["G6"].tree_tip == ""
    assert r["G6"].tree_status.startswith("dropped_by_gtdbtk") and r["G1"].tree_status == "in_tree"
    assert (r["G1"].methanogen, r["G2"].is_anme, r["G3"].methanogen) == ("True", "True", "False")
    assert r["G3"].methanogenesis_subtype == "NA" and r["G2"].mcr_group == "ANME"
    assert r["G1"].gtdb_phylum == "P0" and r["G1"].gtdb_species == "Gen sp0"
    assert (r["G2"].group_Human, r["G3"].niche_Gut, r["G6"].niche_Gut) == ("", "1", "")
    assert r["G1"].n_contigs == "10" and r["G1"].representative == "True"
    miss = {x.column: x for x in read_tsv(d / "miss.tsv").itertuples()}
    assert (miss["group_Human"].n_1, miss["group_Human"].n_0, miss["group_Human"].n_NA) == ("2", "2", "2")
    e = read_tsv(d / "eco.tsv")
    assert e.columns.tolist() == ["genome_id", "group_Human", "group_Marine", "niche_Gut"]
    assert "1 non-binary ecology cell" in p.stderr and "envo-map" in p.stderr
    p = master(w, d, "--primary-tiers", "A,B", check=True)
    m2 = {x.genome_id: x for x in read_tsv(d / "master.tsv").itertuples()}
    assert m2["G2"].prometheus_status_primary == "present"


@suite.test
def test_master_table_inconsistent_inputs_fail():
    d = tmpdir()
    w = master_world(d, catalog_rows=[["l1", "G1", "p1", "D"]])
    fails(master(w, d), "not one of")                                   # invalid tier
    w = master_world(tmpdir(), catalog_rows=[["l1", "G9", "p1", "A"]])
    fails(master(w, d), "outside the genome universe")                  # foreign genome
    w = master_world(tmpdir(), tree="((G1:1,G2:1):1,(G3:1,G9:1):1,G5:1);")
    fails(master(w, d), "species tree")
    w = master_world(tmpdir(), meta_ids=GIDS + ["G1"])
    fails(master(w, d), "duplicated genome_id")
    w = master_world(tmpdir(), eco_rows=[["G1", "x", 1, 0, 1]])
    fails(master(w, d), "missing from the species table")
    w = master_world(tmpdir())
    os.remove(w["catalog"])
    fails(master(w, d), "prometheus catalog not found")
    w = master_world(tmpdir())
    shutil.rmtree(w["df_dir"] / "G3")
    fails(master(w, d), "DefenseFinder systems for G3")
    w = master_world(tmpdir())
    fails(master(w, d, "--primary-tiers", "A,Z"), "tier set")


@suite.test
def test_master_table_empty_universe_is_header_only():
    d = tmpdir()
    w = master_world(d)
    tsv(w["meta"], ["genome_id", "completeness"], [])
    tsv(w["catalog"], ["locus_id", "genome_id", "protein_id", "confidence_tier"], [])
    tsv(w["methanogen"], ["genome_id", "mcrA"], [])
    tsv(w["sensitivity"], ["genome_id", "detect_defensefinder"], [])
    tsv(w["taxonomy"], ["user_genome", "classification"], [])
    tsv(w["membership"], ["genome_id", "in_msa", "in_tree", "status"], [])
    write(w["tree"], "();\n")
    write(w["reps"], "genome_id\n")
    master(w, d, check=True)
    hdr = (d / "master.tsv").read_text().splitlines()
    assert len(hdr) == 1 and hdr[0].split("\t")[:len(SCHEMA_BACKBONE)] == SCHEMA_BACKBONE
    assert "niche_Gut" in hdr[0]


# --------------------------------------------------------------------------- #
# scaffolds.smk rule shells: real Snakemake + stub GTDB-Tk + real IQ-TREE PMSF
# --------------------------------------------------------------------------- #
STUB_GTDBTK = r"""#!/bin/bash
# Faithful to GTDB-Tk 2.4.0 where the rule depends on it: every output names genomes by
# the batchfile id (column 2), and identify ABORTS when a user id equals a GTDB reference
# id -- which include the bare NCBI accession of every GTDB genome (gtdbtk/tools.py
# get_reference_ids; main.py _genomes_to_process).
echo "$@" >> "$STUB_LOG"
sub=$1; shift
out=""; bf=""
while [ $# -gt 0 ]; do case "$1" in --out_dir) out=$2; shift;; --batchfile) bf=$2; shift;; esac; shift; done
if [ "$STUB_MODE" = fail_identify ] && [ "$sub" = identify ]; then echo "stub: identify crashed" >&2; exit 1; fi
if [ "$sub" = identify ]; then
  if cut -f2 "$bf" | grep -Eq '^(RS_|GB_)?GC[AF]_[0-9]{9}[.][0-9]+$'; then
    echo "You have genomes with the same id as GTDB-Tk reference genomes" >&2; exit 1
  fi
  mkdir -p "$out/identify"
  printf "user_GFAIL\tno genes called\n" > "$out/identify/gtdbtk.failed_genomes.tsv"
  exit 0
fi
mkdir -p "$out/align/intermediate_results"
if [ "$STUB_MODE" = nomsa ]; then exit 0; fi
printf "user_GFILT\tInsufficient number of amino acids in MSA (4.2%%)\n" > "$out/align/gtdbtk.ar53.filtered.tsv"
printf ">user_GBAC\nMKV\n" | gzip > "$out/align/gtdbtk.bac120.user_msa.fasta.gz"
printf "user_GNOM\tNo bacterial or archaeal marker\n" > "$out/align/intermediate_results/gtdbtk.align.failed.tsv"
awk '/^>/ {sub(/^>/, ">user_")} {print}' "$STUB_MSA" | gzip -c > "$out/align/gtdbtk.ar53.user_msa.fasta.gz"
"""


# Snakemake 8 hashes each rule's conda env for its job metadata even without --use-conda,
# which runs `conda info --json`, `conda --version` and `conda config --get
# channel_priority --json`; this stub answers only those.
STUB_CONDA = r"""#!/bin/bash
case "$1" in
  info) echo '{"conda_prefix": "@PREFIX@", "platform": "osx-arm64", "envs_dirs": [], "pkgs_dirs": []}';;
  --version) echo "conda 24.9.2";;
  config) echo '{"get": {"channel_priority": "strict"}}';;
  *) echo "stub conda: unsupported: $*" >&2; exit 1;;
esac
"""


def snakemake_bin():
    cands = [os.environ.get("PROM_SNAKEMAKE"), shutil.which("snakemake"),
             str(Path(PYTHON).parents[2] / "prom-smk" / "bin" / "snakemake")]
    for c in cands:
        if c and os.path.exists(c):
            return c
    raise AssertionError("snakemake not found (set PROM_SNAKEMAKE)")


def rule_world(d):
    d = Path(d)
    shutil.copytree(REPO / "workflow", d / "workflow", ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(REPO / "config", d / "config")
    seqs = alisim(PMSF_TRUE, 200, 13, d, "stubmsa")
    stub_msa = fasta(d / "stub_msa.fa", sorted(seqs.items()))
    analysed = sorted(seqs) + ["GBAC", "GFAIL", "GFILT", "GNOM"]
    pre = d / "R" / "preprocess"
    for g in analysed:
        write(pre / "genomes" / f"{g}.fna", ">c1\nACGT\n")
    write(pre / "genomes" / ".staged", "genome_id\n" + "\n".join(analysed) + "\n")
    cls = {g: "d__Archaea;p__Methanobacteriota;c__M" for g in analysed}
    cls.update(OG1="d__Archaea;p__Altiarchaeota;c__A", OG2="d__Archaea;p__Altiarchaeota;c__A",
               XA="d__Archaea;p__Altiarchaeota_A;c__B")
    tsv(pre / "species_metadata.tsv", ["genome_id", "classification", "group_Marine"],
        [[g, cls[g], i % 2] for i, g in enumerate(analysed)])
    tsv(pre / "taxonomy.tsv", ["user_genome", "classification"], [[g, cls[g]] for g in analysed])
    write(pre / "representatives.txt", "\n".join(analysed) + "\n")
    write(pre / "genome_metadata.tsv", "genome_id\n")
    (d / "gtdbdb").mkdir()
    write(d / "cfg.yaml", "\n".join([
        'samples: "R/preprocess/species_metadata.tsv"',
        'metadata: "R/preprocess/species_metadata.tsv"',
        'genomes_dir: "R/preprocess/genomes"', 'genome_extension: "fna"',
        'results_dir: "R"', 'logs_dir: "L"',
        'input: {metadata: "no_phase1_input.tsv", genomes_dir: "", genome_extension: "fasta"}',
        'db: {gtdbtk_data: "gtdbdb"}',
        'tree: {outgroup_taxon: "p__Altiarchaeota", ufboot: 1000, alrt: 1000, '
        'pmsf_model: "LG+C10+F+G"}',   # PMSF branch; C10 keeps the suite fast (unit test: C20)
        'threads: {small: 1, medium: 2, large: 2}']) + "\n")
    b = d / "bin"
    b.mkdir()
    write(b / "gtdbtk", STUB_GTDBTK).chmod(0o755)
    write(b / "conda", STUB_CONDA.replace("@PREFIX@", str(d / "condaprefix"))).chmod(0o755)
    os.symlink(PYTHON, b / "python")
    e = env()
    e.update(PATH=str(b) + os.pathsep + e["PATH"], HOME=str(d / "home"),
             XDG_CACHE_HOME=str(d / "home" / ".cache"), STUB_LOG=str(d / "stub.log"),
             STUB_MSA=str(stub_msa))
    (d / "home").mkdir()
    return seqs, analysed, e


def smk(d, e, *args, timeout=900):
    cmd = [snakemake_bin(), *args, "-s", "workflow/Snakefile", "--configfile", "cfg.yaml",
           "--allowed-rules", "gtdbtk_align", "species_tree", "--nolock"]
    return subprocess.run(cmd, cwd=d, env=e, text=True, capture_output=True, timeout=timeout)


@suite.test
def test_rules_gtdbtk_align_and_species_tree_end_to_end():
    need("iqtree2")
    d = tmpdir()
    seqs, analysed, e = rule_world(d)
    p = smk(d, e, "R/scaffolds/species_tree.treefile", "--cores", "2")
    assert p.returncode == 0, (p.stdout + p.stderr)[-3000:]
    calls = (d / "stub.log").read_text().splitlines()
    assert calls[0].startswith("identify --batchfile ") and "--cpus 2" in calls[0], calls
    assert calls[1].startswith("align ") and "--skip_gtdb_refs" in calls[1] \
        and "--min_perc_aa 10" in calls[1], calls
    batch = [ln.split("\t") for ln in (d / "R/scaffolds/gtdbtk_align/batchfile.tsv").read_text().splitlines()]
    gdir = os.path.realpath(d / "R" / "preprocess" / "genomes")
    # GTDB-Tk sees 'user_<genome_id>' (reference-id collision guard); outputs are stripped.
    assert sorted(b[1] for b in batch) == sorted("user_" + g for g in analysed)
    assert all(os.path.realpath(os.path.dirname(b[0])) == gdir and b[1].startswith("user_")
               and b[0].endswith(f"/{b[1][len('user_'):]}.fna") for b in batch), batch[:3]
    with gzip.open(d / "R/scaffolds/gtdbtk_align/user_msa.fasta.gz", "rt") as fh:
        msa_ids = [ln[1:].split()[0] for ln in fh if ln.startswith(">")]
    assert msa_ids and not any(i.startswith("user_") for i in msa_ids), msa_ids[:3]
    dropped = dict(ln.split("\t", 1) for ln in
                   (d / "R/scaffolds/gtdbtk_align/gtdbtk_dropped.tsv").read_text().splitlines())
    assert dropped == {"GFILT": "ar53 --min_perc_aa filter: Insufficient number of amino acids in MSA (4.2%)",
                       "GBAC": "assigned to the bac120 marker set by GTDB-Tk",
                       "GNOM": "No bacterial or archaeal marker",
                       "GFAIL": "GTDB-Tk identify failed: no genes called"}, dropped
    mem = {r.genome_id: r.status for r in read_tsv(d / "R/scaffolds/species_tree/species_tree_membership.tsv").itertuples()}
    assert {g for g, s in mem.items() if s == "in_tree"} == set(seqs)
    assert mem["GBAC"] == "dropped_by_gtdbtk: assigned to the bac120 marker set by GTDB-Tk"
    t = dtree(d / "R/scaffolds/species_tree.treefile")
    assert {lf.taxon.label for lf in t.leaf_node_iter()} == set(seqs)
    assert frozenset({"OG1", "OG2"}) in root_children(t)
    rooting = (d / "R/scaffolds/species_tree/rooting_method.txt").read_text()
    assert rooting.startswith("method=outgroup(2 tips;monophyletic)"), rooting
    assert (d / "R/scaffolds/species_tree/species.sitefreq").exists()      # PMSF really ran
    log = (d / "L/scaffolds/species_tree.log").read_text()
    assert "9 sequences in the marker MSA" in log and "LG+C10+F+G" in log


def _rendered_gtdbtk_shell(d, e):
    p = smk(d, e, "R/scaffolds/gtdbtk_align/user_msa.fasta.gz", "-n", "-p",
            "--forcerun", "gtdbtk_align", "--cores", "1")
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]
    txt = p.stdout
    i = txt.index("        out=")
    j = txt.index("genomes dropped by GTDB-Tk", i)
    return txt[i:txt.index("\n", j)]


@suite.test
def test_rule_gtdbtk_align_shell_fails_loudly():
    """The rendered rule shell (run as Snakemake runs it: bash, set -euo pipefail) fails
    when GTDB-Tk crashes, when no MSA is produced, and when the database is missing."""
    d = tmpdir()
    _, _, e = rule_world(d)
    script = _rendered_gtdbtk_shell(d, e)
    out_msa = d / "R/scaffolds/gtdbtk_align/user_msa.fasta.gz"
    (d / "L/gtdbtk").mkdir(parents=True, exist_ok=True)        # Snakemake creates log/output
    out_msa.parent.mkdir(parents=True, exist_ok=True)           # directories before a job runs
    for mode, needle in (("fail_identify", "stub: identify crashed"), ("nomsa", "GTDB-Tk wrote no")):
        e2 = dict(e, STUB_MODE=mode)
        p = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], cwd=d, env=e2,
                           text=True, capture_output=True)
        assert p.returncode != 0, (mode, p.stdout, p.stderr)
        assert not out_msa.exists(), mode
        assert needle in (d / "L/gtdbtk/align.log").read_text(), mode
    os.rename(d / "gtdbdb", d / "gtdbdb_moved")
    p = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], cwd=d, env=e,
                       text=True, capture_output=True)
    assert p.returncode != 0 and "not a GTDB-Tk reference data directory" in \
        (d / "L/gtdbtk/align.log").read_text()


def _rendered_species_tree_shell(d, e):
    p = smk(d, e, "R/scaffolds/species_tree.treefile", "-n", "-p", "--forcerun", "species_tree",
            "--cores", "1")
    assert p.returncode == 0, (p.stdout + p.stderr)[-2000:]
    txt = p.stdout
    i = txt.index("        wd=")
    j = txt.index("--log R/scaffolds/species_tree/rooting_method.txt", i)
    return txt[i:txt.index("\n", j)]


@suite.test
def test_rule_species_tree_shell_iqtree_crash_and_small_msa():
    """A crashing IQ-TREE fails the rule (no '();' fallback); an MSA with < 4 genomes is
    the one documented nonfatal case (trivial tree + membership says why)."""
    d = tmpdir()
    seqs, analysed, e = rule_world(d)
    ga = d / "R/scaffolds/gtdbtk_align"
    ga.mkdir(parents=True)
    with gzip.open(ga / "user_msa.fasta.gz", "wt") as fh:
        fh.write("".join(f">{k}\n{v}\n" for k, v in sorted(seqs.items())))
    write(ga / "gtdbtk_dropped.tsv", "GBAC\tassigned to the bac120 marker set by GTDB-Tk\n")
    script = _rendered_species_tree_shell(d, e)
    (d / "L/scaffolds").mkdir(parents=True, exist_ok=True)
    (d / "R/scaffolds/species_tree").mkdir(parents=True, exist_ok=True)
    e2 = dict(e)
    e2["PATH"] = str(fake_tool(d, "iqtree2", 2, "ERROR: simulated IQ-TREE crash")) + os.pathsep + e["PATH"]
    p = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], cwd=d, env=e2,
                       text=True, capture_output=True)
    assert p.returncode != 0, (p.stdout, p.stderr)
    assert not (d / "R/scaffolds/species_tree.treefile").exists()
    assert "simulated IQ-TREE crash" in (d / "L/scaffolds/species_tree.log").read_text()
    with gzip.open(ga / "user_msa.fasta.gz", "wt") as fh:         # 3 genomes only
        fh.write("".join(f">{k}\n{seqs[k]}\n" for k in ("G1", "G2", "G3")))
    p = subprocess.run(["bash", "-c", "set -euo pipefail\n" + script], cwd=d, env=e,
                       text=True, capture_output=True)
    assert p.returncode == 0, (p.stdout, p.stderr)
    assert (d / "R/scaffolds/species_tree.treefile").read_text().strip() == "();"
    mem = {r.genome_id: r.status for r in
           read_tsv(d / "R/scaffolds/species_tree/species_tree_membership.tsv").itertuples()}
    assert mem["G1"].startswith("msa_too_small_for_tree") and mem["GBAC"].startswith("dropped_by_gtdbtk")
    assert mem["OG1"] == "not_in_gtdbtk_msa: no reason reported"


@suite.test
def test_mcr_search_one_subunit_per_protein():
    """A protein accepted by two subunit models counts once, for the model it exceeds
    most (score / threshold); the other call is recorded as superseded."""
    need("hmmsearch")
    d = tmpdir()
    hd = d / "hmm"
    hd.mkdir()
    probe = build_hmm(FIX / "uniprot_mcrA.faa", "probe", d)
    faa, pm, gdir = mcr_world(d)
    full = independent_scores(probe, faa, d)["G1__G1_1_1"][0]
    build_hmm(FIX / "uniprot_mcrA.faa", "TIGR03256", hd, ga=(0.8 * full, 0.8 * full))   # mcrA
    build_hmm(FIX / "uniprot_mcrA.faa", "TIGR03257", hd, ga=(0.9 * full, 0.9 * full))   # 'mcrB'
    allh = d / "all.tsv"
    mcr("--proteome", faa, "--pmap", pm, "--hmmdir", hd, "--kos", "K00399,K00401",
        "--threads", THREADS, "--out", d / "hits.tsv", "--out-all", allh, check=True)
    h = read_tsv(d / "hits.tsv")
    g1 = h[h.genome_id == "G1"]
    assert set(g1.gene) == {"mcrA"}, g1[["gene", "protein_id"]]   # McrA not double-counted as McrB
    dec = {(r.protein_id, r.gene): r.decision for r in read_tsv(allh).itertuples()}
    assert dec[("G1__G1_1_1", "mcrA")] == "accepted_full"
    assert dec[("G1__G1_1_1", "mcrB")] == "superseded_by_mcrA"


if __name__ == "__main__":
    suite.main("SCAFFOLDS TESTS PASSED")
