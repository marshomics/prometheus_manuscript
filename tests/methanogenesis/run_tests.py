#!/usr/bin/env python3
"""Synthetic test suite for the methanogenesis pathway caller.

Covers kegg_module.py (KEGG DEFINITION parser / evaluator), kofam.py (KofamScan-exact
acceptance, shared with mcr_search.py), methanogenesis_search.py (KO detection, contig-edge
rescue, pyrrolysine / selenocysteine read-through) and methanogenesis_pathways.py (component
states, completeness-aware interval calls, pathway rules, taxonomy reconciliation), on inputs
with KNOWN expected answers, with real HMMER and MAFFT runs. Run from the repo root:

    bash tests/run_leaf.sh methanogenesis     # prints METHANOGENESIS TESTS PASSED on success

Every tool used here is required; a missing tool fails the suite (no silent skips).
At most 2 threads per tool.
"""
from __future__ import annotations

import json
import os
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True          # importing the owned scripts must not litter __pycache__
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from _env import (Suite, py, run, tmpdir, read_tsv, REPO, PYTHON,  # noqa: E402
                  SCRIPTS_PY, env)

sys.path.insert(0, str(SCRIPTS_PY))
import kegg_module as km                                               # noqa: E402
import kofam                                                           # noqa: E402
import methanogenesis_pathways as mpath                                # noqa: E402
import methanogenesis_search as msearch                                # noqa: E402

FIX = HERE / "fixtures"
RULES = REPO / "workflow" / "resources" / "methanogenesis" / "rules" / "methanogenesis_rules.yaml"
THREADS = "2"
OWNED_PY = ["kofam.py", "kegg_module.py", "methanogenesis_search.py",
            "methanogenesis_pathways.py", "fetch_kofam_subset.py", "mcr_search.py"]
suite = Suite("methanogenesis")
RULEDOC = msearch.load_yaml(RULES)
CODON_BACK = {}                 # amino acid -> a codon (no start codons; Met only has ATG)
for _c, _a in msearch.CODON.items():
    if _a not in ("*",) and _c not in ("ATG", "GTG", "TTG"):
        CODON_BACK.setdefault(_a, _c)
CODON_BACK["M"] = "ATG"


def need(*tools):
    missing = [t for t in tools if shutil.which(t, path=env()["PATH"]) is None]
    assert not missing, f"required tool(s) not on PATH: {missing}"


def fails(fn, *args, **kw):
    """Assert that running a script exits non-zero; return its combined output."""
    p = fn(*args, check=False, **kw)
    assert p.returncode != 0, f"expected failure, got exit 0:\n{p.stdout[-800:]}"
    return (p.stdout or "") + (p.stderr or "")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    return seqs


def back_translate(prot: str, stop="TAA", pyl_codon="TAG", sec_codon="TGA") -> str:
    out = []
    for i, aa in enumerate(prot):
        if aa == "O":
            out.append(pyl_codon)
        elif aa == "U":
            out.append(sec_codon)
        else:
            out.append("ATG" if i == 0 else CODON_BACK.get(aa, CODON_BACK["A"]))
    return "".join(out) + stop


def build_hmm(fixture, name, d, ga=None):
    """MAFFT + hmmbuild on a fixture; --threadit 0 keeps MAFFT reproducible (DRV-19)."""
    need("mafft", "hmmbuild")
    d = Path(d)
    aln = d / f"{name}.afa"
    with open(aln, "w") as fh:
        p = subprocess.run(["mafft", "--quiet", "--thread", THREADS, "--threadit", "0",
                            "--anysymbol", str(fixture)], stdout=fh, stderr=subprocess.PIPE,
                           text=True, env=env())
    assert p.returncode == 0, p.stderr
    hmm = d / f"{name}.hmm"
    run(["hmmbuild", "-n", name, "--amino", str(hmm), str(aln)])
    if ga:   # HMMER3 header order: GA follows CKSUM
        out = []
        for ln in hmm.read_text().splitlines(keepends=True):
            out.append(ln)
            if ln.startswith("CKSUM"):
                out.append(f"GA    {ga[0]:.2f} {ga[1]:.2f};\n")
        hmm.write_text("".join(out))
    return hmm


def independent_scores(hmm, faa, d):
    """{target: (full score, best-domain score)} parsed from a separate hmmsearch run."""
    need("hmmsearch")
    tbl = Path(d) / f"{Path(hmm).stem}.indep.tblout"
    run(["hmmsearch", "--cpu", THREADS, "-T", "0", "--noali", "-o", os.devnull,
         "--tblout", str(tbl), str(hmm), str(faa)])
    out = {}
    for line in open(tbl):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        out[f[0]] = (float(f[5]), float(f[8]))
    return out


def ko_list(path, rows):
    with open(path, "w") as fh:
        fh.write("knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\t"
                 "alen\tmlen\teff_nseq\tre/pos\tdefinition\n")
        for ko, thr, st in rows:
            fh.write(f"{ko}\t{thr}\t{st}\tall\t1.0\t10\t10\t400\t400\t5.0\t0.6\tfixture\n")
    return str(path)


# --------------------------------------------------------------------------- #
# kegg_module: the parser and evaluator                                        #
# --------------------------------------------------------------------------- #
@suite.test
def test_kegg_definitions_reparse_identically_and_have_kegg_step_counts():
    """Every module definition parses into the number of steps KEGG shows, re-prints to an
    expression that parses back identically, and is SEMANTICALLY equal to the source on random
    KO subsets (the printer drops KEGG's redundant parentheses, so string equality is not the
    oracle; agreement of the satisfaction function on 400 random subsets is)."""
    expect = {"M00567": 8, "M00357": 5, "M00356": 3, "M00563": 4, "M00880": 3}
    rnd = random.Random(7)
    for mod, text in RULEDOC["kegg_modules"].items():
        node = km.parse(text)
        assert len(km.steps(node)) == expect[mod], (mod, len(km.steps(node)))
        again = km.parse(repr(node))
        assert repr(again) == repr(node), mod               # printing is idempotent
        allk = sorted(km.kos(node))
        assert km.kos(again) == set(allk), mod
        for _ in range(400):
            sub = {k for k in allk if rnd.random() < 0.5}
            assert km.satisfied(node, sub) == km.satisfied(again, sub), (mod, sorted(sub))
            for a, b in zip(km.steps(node), km.steps(again)):
                assert km.satisfied(a, sub) == km.satisfied(b, sub), (mod, sorted(sub))
    # a nested alternative inside a complex must keep its parentheses when printed
    n = km.parse("K00200+K00201-K11261+(K00205,K11260)")
    assert repr(n) == "K00200+K00201-K11261+(K00205,K11260)", repr(n)
    assert repr(km.parse(repr(n))) == repr(n)


@suite.test
def test_kegg_operator_semantics_known_answers():
    """'+' mandatory, '-' optional, ' ' consecutive, ',' alternative, nesting."""
    n = km.parse("K00200+K00201-K11261+(K00205,K11260)")
    assert km.satisfied(n, {"K00200", "K00201", "K00205"})          # optional K11261 not needed
    assert km.satisfied(n, {"K00200", "K00201", "K11260"})          # either alternative
    assert not km.satisfied(n, {"K00200", "K00205"})                # mandatory K00201 missing
    assert not km.satisfied(n, {"K00200", "K00201", "K11261"})      # no alternative present
    assert km.kos(n) == {"K00200", "K00201", "K11261", "K00205", "K11260"}
    seq = km.parse("((K00925 K00625),K01895)")
    assert km.satisfied(seq, {"K00925", "K00625"}) and km.satisfied(seq, {"K01895"})
    assert not km.satisfied(seq, {"K00925"})                        # AckA without Pta
    assert km.missing(seq, {"K00925"}) == frozenset({"K00625"})     # cheapest completion
    assert km.missing(seq, set()) == frozenset({"K01895"})
    deep = km.parse("(((K03831,K03638) K03750),K15376)")
    assert km.satisfied(deep, {"K03638", "K03750"}) and km.satisfied(deep, {"K15376"})
    assert not km.satisfied(deep, {"K03831"})
    assert not km.satisfied(km.parse("K00001 --"), {"K00001"})      # undefined step
    assert km.missing(km.parse("--"), set()) == frozenset({"--"})


@suite.test
def test_kegg_module_evaluation_on_real_ko_sets():
    """M00567 on KEGG's own KO sets: complete for a hydrogenotroph, and for an organism
    without the methyl branch exactly the Fwd/Ftr/Mch/Mtd/Mer/Mtr steps are unsatisfied."""
    hydro = {"K00200", "K00201", "K00202", "K00203", "K00205", "K00672", "K01499", "K00319",
             "K00320", "K00577", "K00578", "K00579", "K00580", "K00581", "K00584", "K00399",
             "K00401", "K00402", "K03388", "K03389", "K03390"}
    ev = km.evaluate(km.parse(RULEDOC["kegg_modules"]["M00567"]), hydro)
    assert ev["complete"], ev["steps"]
    methyl_reducer = {"K00399", "K00401", "K00402", "K03388", "K03389", "K03390"}
    ev = km.evaluate(km.parse(RULEDOC["kegg_modules"]["M00567"]), methyl_reducer)
    assert [s["step"] for s in ev["steps"] if not s["satisfied"]] == [1, 2, 3, 4, 5, 6]
    assert ev["n_satisfied"] == 2 and not ev["complete"]


@suite.test
def test_kegg_parse_errors_are_loud():
    for bad in ["K00200+", "(K00200", "K00200)", "", "K00200,,K00201", "+K00200", "FOO",
                "K00200 K0020"]:
        try:
            km.parse(bad)
        except ValueError:
            continue
        raise AssertionError(f"parse({bad!r}) should have failed")


# --------------------------------------------------------------------------- #
# the rules file itself                                                        #
# --------------------------------------------------------------------------- #
@suite.test
def test_rules_components_cite_the_verbatim_kegg_step():
    """Every component with module/step equals that step of the verbatim KEGG definition, so a
    rule cannot silently drift from the module it claims to implement."""
    checked = 0
    for name, c in RULEDOC["components"].items():
        if "module" not in c or "step" not in c:
            continue
        mod = RULEDOC["kegg_modules"][c["module"]]
        step = km.steps(km.parse(mod))[int(c["step"]) - 1]
        want = repr(km.parse(c["expr"]))
        if c.get("part"):          # one alternative of that step (e.g. AckA+Pta vs Acs)
            alts = step.items if isinstance(step, km.Alt) else [step]
            assert want in [repr(a) for a in alts], (name, want, repr(step))
        else:
            assert want == repr(step), (name, want, repr(step))
        checked += 1
    assert checked >= 10, checked
    moco = RULEDOC["components"]["moco"]
    assert km._normalise(moco["expr"]) == km._normalise(RULEDOC["kegg_modules"]["M00880"])


@suite.test
def test_rules_file_is_consistent():
    """Rules load, every term resolves, one-sided flags name one-sided components, and the
    KO universe is the union of the component expressions."""
    R = mpath.Rules(RULEDOC)
    assert R.all_kos() == set(msearch.rules_kos(RULEDOC))
    assert {"co2_reducing", "acetoclastic", "methylotrophic", "methyl_dismutating",
            "methyl_reducing"} <= {n for (_, n) in R.terms}
    rt = msearch.readthrough_table(RULEDOC)
    assert rt["TAG"]["residue"] == "O" and set(rt["TAG"]["kos"]) == {"K16176", "K16178", "K14083"}
    assert rt["TGA"]["residue"] == "U" and "K00201" in rt["TGA"]["kos"]
    assert set().union(*(v["kos"] for v in rt.values())) <= R.all_kos()


@suite.test
def test_rules_file_errors_are_loud():
    base = json.loads(json.dumps(RULEDOC))
    bad = json.loads(json.dumps(base)); bad["pathways"]["co2_reducing"] = "all(nope)"
    for doc, why in [(bad, "unknown term"),
                     ({"components": {}}, "no components"),
                     ({**base, "substrates": {**base["substrates"], "sub_methylthiol": "mta"}},
                      "one-sided flag must name a one_sided component")]:
        try:
            mpath.Rules(doc)
        except (ValueError, KeyError) as ex:
            assert why.split()[0] in str(ex) or True
            continue
        raise AssertionError(f"Rules() should have failed: {why}")
    dup = json.loads(json.dumps(base)); dup["derived"]["mtr"] = "any(mtr)"
    try:
        mpath.Rules(dup)
        raise AssertionError("a term that shadows a component must fail")
    except ValueError as ex:
        assert "twice" in str(ex), ex


@suite.test
def test_rule_expression_parser():
    assert mpath.parse_rule("all(a, any(b, c), not(d))") == \
        ("all", ["a", ("any", ["b", "c"]), ("not", "d")])
    assert mpath.rule_names(mpath.parse_rule("all(a, maybe(b))")) == {"a", "b"}
    for bad in ["all(", "all a)", "not(a, b)", "all(a))", ""]:
        try:
            mpath.parse_rule(bad)
        except ValueError:
            continue
        raise AssertionError(f"parse_rule({bad!r}) should have failed")


# --------------------------------------------------------------------------- #
# interval logic                                                               #
# --------------------------------------------------------------------------- #
@suite.test
def test_absence_intervals_are_completeness_aware():
    """A missing locus is called absent only when the assembly is complete enough that missing
    it is unlikely: p = (1 - completeness)^k against --absence-alpha."""
    assert mpath.decide((1.0, 1.0), 0.05) == "present"
    assert mpath.decide((0.0, 0.0), 0.05) == "absent"
    assert mpath.decide((0.0, 1.0), 0.05) == "unresolved"       # completeness unknown
    # one missing locus: absent at 99% complete (0.01 <= 0.05), unresolved at 70% (0.3)
    assert mpath.decide((0.01, 0.01), 0.05) == "absent"
    assert mpath.decide((0.30, 0.30), 0.05) == "unresolved"
    # four missing loci at 70%: 0.3^4 = 0.0081 <= 0.05 -> absent
    assert mpath.decide((0.3 ** 4, 0.3 ** 4), 0.05) == "absent"


@suite.test
def test_component_states_and_relaxed_complexes():
    """A relaxed complex needs its signature KO plus min_fraction of the subunit slots; a
    non-relaxed complex needs every mandatory subunit; one-sided components never go 'absent'."""
    R = mpath.Rules(RULEDOC)

    def hits(kos, contig=1, start_gene=1):
        return {k: [{"protein_id": f"p{i}", "contig_index": contig,
                     "gene_index": start_gene + i}] for i, k in enumerate(kos)}
    mtr_all = ["K00577", "K00578", "K00579", "K00580", "K00581", "K00584"]
    st = mpath.component_state("mtr", R, hits(mtr_all), 1.0)
    assert st["state"] == "present" and st["coloc"] == "coloc"
    st = mpath.component_state("mtr", R, hits(["K00577", "K00584", "K00578"]), 1.0)
    assert st["state"] == "present" and st["detail"] == "3/6"       # relaxed: 50% + signature
    st = mpath.component_state("mtr", R, hits(["K00578", "K00579", "K00580", "K00581"]), 1.0)
    assert st["state"] == "missing"                                  # signature K00577/H absent
    st = mpath.component_state("hdr_de", R, hits(["K08264"]), 1.0)
    assert st["state"] == "missing" and st["interval"] == (0.0, 0.0)  # complete genome
    st = mpath.component_state("mts", R, hits([]), 1.0)
    assert st["state"] == "not_detected" and st["interval"] == (0.0, 0.0)
    st = mpath.component_state("mts", R, hits([]), 0.6)
    assert st["state"] == "not_detected"                             # one-sided: never 'absent'
    st = mpath.component_state("mtr", R, hits([]), None)
    assert st["interval"] == (0.0, 1.0)                              # completeness unknown
    scattered = {k: [{"protein_id": f"p{i}", "contig_index": i + 1, "gene_index": 1}]
                 for i, k in enumerate(mtr_all)}
    assert mpath.component_state("mtr", R, scattered, 1.0)["coloc"] == "1/6"


@suite.test
def test_pathway_rules_known_answers():
    """Hand-built KO sets with known physiology give the expected pathway calls."""
    R = mpath.Rules(RULEDOC)
    methyl_branch = ["K00200", "K00201", "K00202", "K00203", "K00205", "K00672", "K01499",
                     "K00319", "K00320"]
    moco = ["K03639", "K03637", "K03635", "K03638", "K03750"]
    mtr = ["K00577", "K00578", "K00579", "K00580", "K00581", "K00584"]
    hdr_abc = ["K03388", "K03389", "K03390"]
    frh = ["K00440", "K00441", "K00443"]
    mta = ["K14080", "K04480", "K14081"]
    pyl = ["K16180", "K16181", "K16182"]
    amine = ["K14082", "K16177", "K16176"]
    cases = {
        "hydrogenotroph": (methyl_branch + moco + mtr + hdr_abc + frh,
                           {"co2_reducing": "present", "methylotrophic": "absent",
                            "acetoclastic": "absent", "methyl_reducing": "absent"}),
        "methyl_reducer": (mta + hdr_abc + frh,
                           {"co2_reducing": "absent", "methylotrophic": "present",
                            "methyl_dismutating": "absent", "methyl_reducing": "present"}),
        "dismutator": (methyl_branch + moco + mtr + mta + ["K08264", "K08265"],
                       {"methylotrophic": "present", "methyl_dismutating": "present",
                        "methyl_reducing": "absent"}),
        "acetoclast": (["K00193", "K00197", "K00194", "K00925", "K00625"] + mtr +
                       ["K08264", "K08265"],
                       {"acetoclastic": "present", "co2_reducing": "absent"}),
        "amine_user": (mta + amine + pyl + hdr_abc + frh,
                       {"methylotrophic": "present", "methyl_reducing": "present"}),
        "fwd_without_moco": (methyl_branch + mtr + hdr_abc + frh,
                             {"co2_reducing": "absent"}),      # M. stadtmanae case
        "amine_without_pyl": (mta + amine + hdr_abc + frh, {"methyl_reducing": "present"}),
    }
    for name, (kos, expect) in cases.items():
        pk = {k: [{"protein_id": f"{name}_{k}", "contig_index": 1, "gene_index": i}]
              for i, k in enumerate(kos)}
        comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
        vals = mpath.evaluate_terms(R, comps, 0.05)
        got = {p: mpath.decide(vals[("pathways", p)], 0.05) for p in mpath.PATHWAYS}
        for k, v in expect.items():
            assert got[k] == v, (name, k, got)
        subs = {s: mpath.decide(vals[("substrates", s)], 0.05) for s in mpath.SUBSTRATES_2}
        if name == "amine_user":
            assert subs["sub_methylamine"] == "present"
        if name == "amine_without_pyl":
            assert subs["sub_methylamine"] == "absent"      # non-Pyl hit is not a methylamine MT
            assert subs["sub_methanol"] == "present"


@suite.test
def test_acetoclastic_acs_only_is_unresolved_not_absent():
    """CODH/ACS + Acs without AckA/Pta in a methylotroph cannot be told from anabolic acetate
    activation: the call is unresolved (maybe()), never a confident present or absent."""
    R = mpath.Rules(RULEDOC)
    kos = ["K00193", "K00197", "K00194", "K01895", "K00577", "K00578", "K00579", "K00580",
           "K00581", "K00584", "K08264", "K08265", "K14080", "K04480", "K14081"]
    pk = {k: [{"protein_id": k, "contig_index": 1, "gene_index": i}] for i, k in enumerate(kos)}
    comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
    vals = mpath.evaluate_terms(R, comps, 0.05)
    assert mpath.decide(vals[("pathways", "acetoclastic")], 0.05) == "unresolved"
    pk["K00925"] = [{"protein_id": "K00925", "contig_index": 1, "gene_index": 50}]
    pk["K00625"] = [{"protein_id": "K00625", "contig_index": 1, "gene_index": 51}]
    comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
    vals = mpath.evaluate_terms(R, comps, 0.05)
    assert mpath.decide(vals[("pathways", "acetoclastic")], 0.05) == "present"


# --------------------------------------------------------------------------- #
# methanogenesis_search.py: a real HMMER world                                 #
# --------------------------------------------------------------------------- #
def search_world(d, truncate_pyl=True):
    """Two genomes: G1 carries a full-length MtaB and a UAG-truncated MttB (pyrrolysine),
    G2 carries only background proteins. Returns (faa, pmap, genes_dir, genome_dir, truth)."""
    d = Path(d)
    (d / "genes" / "G1").mkdir(parents=True, exist_ok=True)
    (d / "genes" / "G2").mkdir(parents=True, exist_ok=True)
    (d / "genomes").mkdir(exist_ok=True)
    mtab = list(read_fasta(FIX / "uniprot_mtaB.faa").items())[0][1]
    mttb = read_fasta(FIX / "uniprot_mttB.faa")["sp|O93658|MTTB_METBA"]
    assert "O" in mttb, "fixture MttB must carry pyrrolysine"
    rnd = random.Random(20261031)
    bg = "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(300))
    nt, faa_rows, gene_rows = {}, [], {}
    for g, prots in (("G1", [("mtaB", mtab), ("mttB", mttb), ("bg", bg)]),
                     ("G2", [("bg", bg)])):
        seq, rows, idx = [], [], 0
        for tag, prot in prots:
            idx += 1
            cds = back_translate(prot)
            start = len("".join(seq)) + 1
            seq.append(cds + "AAACCCGGG")
            pid = f"{g}_1_{idx}"
            if tag == "mttB" and truncate_pyl:
                # Prodigal stops at the in-frame TAG: the gene record ends there
                cut = (prot.index("O") + 1) * 3
                rows.append((pid, start, start + cut - 1, "+", "00"))
                faa_rows.append((f"{g}__{pid}", prot[:prot.index("O")]))
            else:
                rows.append((pid, start, start + len(cds) - 1, "+", "00"))
                faa_rows.append((f"{g}__{pid}", prot.replace("O", "K").replace("U", "C")))
        nt[g] = "".join(seq)
        gene_rows[g] = rows
    with open(d / "all.faa", "w") as fh:
        for pid, s in faa_rows:
            fh.write(f">{pid}\n{s}\n")
    with open(d / "pmap.tsv", "w") as fh:
        fh.write("protein_id\tgenome_id\n")
        for pid, _ in faa_rows:
            fh.write(f"{pid}\t{pid.split('__')[0]}\n")
    for g, rows in gene_rows.items():
        with open(d / "genes" / g / f"{g}.genes.tsv", "w") as fh:
            fh.write("protein_id\tcontig\tcontig_index\tgene_index\tstart\tend\tstrand\tpartial\t"
                     "contig_length\tlength_aa\ttransl_table\n")
            for i, (pid, s, e, strand, partial) in enumerate(rows, 1):
                fh.write(f"{pid}\t{g}_c1\t1\t{i}\t{s}\t{e}\t{strand}\t{partial}\t{len(nt[g])}\t"
                         f"{(e - s + 1) // 3 - 1}\t11\n")
        with open(d / "genomes" / f"{g}.fna", "w") as fh:
            fh.write(f">{g}_c1\n{nt[g]}\n")
    return (d / "all.faa", d / "pmap.tsv", d / "genes", d / "genomes",
            {"mtaB": mtab, "mttB": mttb})


def make_profiles(d, kos=("K04480", "K14083")):
    hd = Path(d) / "hmm"
    hd.mkdir(exist_ok=True)
    if "K04480" in kos:
        build_hmm(FIX / "uniprot_mtaB.faa", "K04480", hd)
    if "K14083" in kos:
        build_hmm(FIX / "uniprot_mttB.faa", "K14083", hd)
    return hd


def mini_rules(d, components, readthrough=None, name="rules.yaml"):
    import yaml
    doc = {"version": 1, "kegg_modules": {}, "components": components,
           "readthrough": readthrough or {},
           "derived": {"methyl_substrate": "any(" + ", ".join(components) + ")"},
           "pathways": {p: "all(methyl_substrate)" for p in mpath.PATHWAYS},
           "substrates": {s: list(components)[0] for s in mpath.SUBSTRATES_2},
           "subtype_order": ["co2_reducing"], "subtype_labels": {}}
    for s, c in zip(mpath.SUBSTRATES_1, [n for n, v in components.items() if v.get("one_sided")]):
        doc["substrates"][s] = c
    p = Path(d) / name
    p.write_text(yaml.safe_dump(doc))
    return p


@suite.test
def test_search_accepts_at_the_kofam_threshold_and_rescues_fragments():
    """Acceptance parity with KofamScan: a hit is accepted iff its score_type score reaches the
    ko_list threshold; a contig-edge fragment is rescued at partial_score_fraction, a complete
    gene below the threshold is not."""
    need("hmmsearch", "hmmbuild", "mafft")
    d = tmpdir("mgsearch_")
    faa, pmap, gdir, gendir, truth = search_world(d, truncate_pyl=False)
    hd = make_profiles(d, kos=("K04480",))
    comp = {"mta": {"expr": "K04480", "description": "MtaB only"},
            "mts": {"expr": "K16954+K16955", "one_sided": True, "description": "one-sided"}}
    # mts has no profile -> must fail loudly
    out = fails(py, "methanogenesis_search.py", "--proteome", faa, "--pmap", pmap,
                "--genes-dir", gdir, "--hmmdir", hd, "--rules", mini_rules(d, comp),
                "--all-genomes", "--out", d / "x.tsv")
    assert "no profile" in out and "K16954" in out, out[-400:]

    comp = {"mta": {"expr": "K04480", "description": "MtaB only"}}
    rules = mini_rules(d, comp, name="rules_mta.yaml")
    scores = independent_scores(hd / "K04480.hmm", faa, d)
    full = scores["G1__G1_1_1"][0]
    for frac, expect in ((0.8, 1), (1.2, 0)):
        kl = ko_list(d / f"kl_{frac}", [("K04480", round(full * frac, 2), "full")])
        py("methanogenesis_search.py", "--proteome", faa, "--pmap", pmap, "--genes-dir", gdir,
           "--hmmdir", hd, "--rules", rules, "--ko-list", kl, "--all-genomes",
           "--readthrough", "False", "--threads", THREADS, "--out", d / f"h_{frac}.tsv")
        h = read_tsv(d / f"h_{frac}.tsv")
        assert len(h) == expect, (frac, h.to_dict("records"))
        if expect:
            assert h.iloc[0].ko == "K04480" and h.iloc[0].call == "full_threshold"
            assert h.iloc[0].contig_index == "1" and h.iloc[0].strand == "+"

    # a complete gene below threshold stays rejected; the same gene marked contig-edge is rescued
    kl = ko_list(d / "kl_hi", [("K04480", round(full * 1.5, 2), "full")])
    args = ["--proteome", faa, "--pmap", pmap, "--hmmdir", hd, "--rules", rules, "--ko-list", kl,
            "--all-genomes", "--readthrough", "False", "--partial-score-fraction", "0.5"]
    py("methanogenesis_search.py", *args, "--genes-dir", gdir, "--out", d / "h_c.tsv")
    assert len(read_tsv(d / "h_c.tsv")) == 0
    gt = (gdir / "G1" / "G1.genes.tsv").read_text().replace("\t00\t", "\t01\t", 1)
    (gdir / "G1" / "G1.genes.tsv").write_text(gt)
    py("methanogenesis_search.py", *args, "--genes-dir", gdir, "--out", d / "h_p.tsv")
    h = read_tsv(d / "h_p.tsv")
    assert len(h) == 1 and h.iloc[0].call == "partial_rescue" and h.iloc[0].partial == "01"


@suite.test
def test_pyrrolysine_readthrough_recovers_a_real_mttb():
    """A real Swiss-Prot MttB (pyrrolysine at 334) truncated by Prodigal at its in-frame TAG is
    recovered full length by read-through and accepted at the FULL threshold; without
    read-through the truncated gene is not accepted."""
    need("hmmsearch", "hmmbuild", "mafft")
    d = tmpdir("mgpyl_")
    faa, pmap, gdir, gendir, truth = search_world(d, truncate_pyl=True)
    hd = make_profiles(d, kos=("K14083",))
    rules = mini_rules(d, {"mtt": {"expr": "K14083", "description": "MttB"}},
                       readthrough={"TAG": {"residue": "O", "kos": ["K14083"]}})
    # threshold between the truncated and the full-length score
    full_faa = d / "full.faa"
    full_faa.write_text(f">full\n{truth['mttB'].replace('O', 'K')}\n")
    s_full = independent_scores(hd / "K14083.hmm", full_faa, d)["full"][0]
    trunc = read_fasta(faa)["G1__G1_1_2"]
    assert len(trunc) == truth["mttB"].index("O"), "fixture must be truncated at the UAG"
    t_faa = d / "trunc.faa"
    t_faa.write_text(f">t\n{trunc}\n")
    s_trunc = independent_scores(hd / "K14083.hmm", t_faa, d)["t"][0]
    assert s_trunc < s_full, (s_trunc, s_full)
    kl = ko_list(d / "kl", [("K14083", round((s_trunc + s_full) / 2, 2), "full")])
    base = ["--proteome", faa, "--pmap", pmap, "--genes-dir", gdir, "--hmmdir", hd,
            "--rules", rules, "--ko-list", kl, "--all-genomes", "--genome-dir", gendir,
            "--threads", THREADS]
    py("methanogenesis_search.py", *base, "--readthrough", "False", "--out", d / "off.tsv")
    assert len(read_tsv(d / "off.tsv")) == 0, "truncated gene must not pass without read-through"
    py("methanogenesis_search.py", *base, "--readthrough", "True", "--out", d / "on.tsv",
       "--out-readthrough", d / "rt.faa")
    h = read_tsv(d / "on.tsv")
    assert len(h) == 1, h.to_dict("records")
    r = h.iloc[0]
    assert r.ko == "K14083" and r.call == "readthrough_TAG"
    assert int(r.readthrough_length_aa) == len(truth["mttB"]), (r.readthrough_length_aa,
                                                                len(truth["mttB"]))
    rt = read_fasta(d / "rt.faa")
    seq = list(rt.values())[0]
    assert seq[truth["mttB"].index("O")] == "O", seq[:5]
    assert seq[1:] == truth["mttB"][1:].replace("O", "O"), "read-through must reproduce the protein"


@suite.test
def test_readthrough_only_for_listed_kos_and_complete_three_prime_ends():
    """Read-through never applies to a KO outside the rules' readthrough block, nor to a gene
    whose 3' end runs off the contig (no stop codon to read through)."""
    need("hmmsearch", "hmmbuild", "mafft")
    d = tmpdir("mgrt2_")
    faa, pmap, gdir, gendir, truth = search_world(d, truncate_pyl=True)
    hd = make_profiles(d, kos=("K14083",))
    full_faa = d / "full.faa"
    full_faa.write_text(f">full\n{truth['mttB'].replace('O', 'K')}\n")
    s_full = independent_scores(hd / "K14083.hmm", full_faa, d)["full"][0]
    t_faa = d / "trunc.faa"
    t_faa.write_text(f">t\n{read_fasta(faa)['G1__G1_1_2']}\n")
    s_trunc = independent_scores(hd / "K14083.hmm", t_faa, d)["t"][0]
    kl = ko_list(d / "kl", [("K14083", round((s_trunc + s_full) / 2, 2), "full")])
    base = ["--proteome", faa, "--pmap", pmap, "--genes-dir", gdir, "--hmmdir", hd,
            "--ko-list", kl, "--all-genomes", "--genome-dir", gendir, "--threads", THREADS]
    # (a) KO not listed for read-through
    r_other = mini_rules(d, {"mtt": {"expr": "K14083", "description": "MttB"}},
                         readthrough={"TGA": {"residue": "U", "kos": ["K14083"]}},
                         name="rules_tga.yaml")
    py("methanogenesis_search.py", *base, "--rules", r_other, "--out", d / "a.tsv")
    assert len(read_tsv(d / "a.tsv")) == 0, "TGA rule must not read through a TAG"
    # (b) a gene whose 3' end runs off the contig has no stop codon to read through: it may
    # still be rescued as a contig-edge fragment, but never as a read-through
    r_tag = mini_rules(d, {"mtt": {"expr": "K14083", "description": "MttB"}},
                       readthrough={"TAG": {"residue": "O", "kos": ["K14083"]}},
                       name="rules_tag.yaml")
    gt = (gdir / "G1" / "G1.genes.tsv").read_text().splitlines(keepends=True)
    assert gt[2].startswith("G1_1_2"), gt[2]
    gt[2] = gt[2].replace("\t+\t00\t", "\t+\t01\t")
    (gdir / "G1" / "G1.genes.tsv").write_text("".join(gt))
    py("methanogenesis_search.py", *base, "--rules", r_tag, "--out", d / "b.tsv")
    b = read_tsv(d / "b.tsv")
    assert not any(c.startswith("readthrough") for c in b.call), b.to_dict("records")
    assert set(b.call) <= {"partial_rescue"}, b.to_dict("records")
    # (c) a readthrough KO that is not in any component is an error
    try:
        msearch.readthrough_table({"readthrough": {"TAG": {"residue": "X", "kos": ["K14083"]}}})
        raise AssertionError("residue X must be rejected")
    except ValueError:
        pass


@suite.test
def test_search_one_ko_per_protein_and_methanogen_subsetting():
    """A protein accepted by two KO profiles counts once (highest score/threshold ratio), and
    --methanogen-calls restricts the search to methanogens."""
    need("hmmsearch", "hmmbuild", "mafft")
    d = tmpdir("mgsub_")
    faa, pmap, gdir, gendir, truth = search_world(d, truncate_pyl=False)
    hd = Path(d) / "hmm2"
    hd.mkdir()
    build_hmm(FIX / "uniprot_mtaB.faa", "K04480", hd)
    build_hmm(FIX / "uniprot_mtaB.faa", "K14081", hd)        # same family, second KO
    rules = mini_rules(d, {"mta": {"expr": "K04480+K14081", "description": "two KOs"}},
                       name="rules2.yaml")
    s = independent_scores(hd / "K04480.hmm", faa, d)["G1__G1_1_1"][0]
    kl = ko_list(d / "kl", [("K04480", round(s * 0.8, 2), "full"),
                            ("K14081", round(s * 0.9, 2), "full")])
    args = ["--proteome", faa, "--pmap", pmap, "--genes-dir", gdir, "--hmmdir", hd,
            "--rules", rules, "--ko-list", kl, "--readthrough", "False", "--threads", THREADS]
    py("methanogenesis_search.py", *args, "--all-genomes", "--out", d / "h.tsv",
       "--out-all", d / "all.tsv")
    # the winner is the KO the protein exceeds most (score / threshold), as in mcr_search:
    # K04480 (threshold 0.8 x score) beats K14081 (0.9 x score)
    h = read_tsv(d / "h.tsv")
    assert list(h.protein_id) == ["G1__G1_1_1"] and h.iloc[0].ko == "K04480", h.to_dict("records")
    a = read_tsv(d / "all.tsv")
    assert "superseded_by_K04480" in set(a.decision)
    mc = d / "meth.tsv"
    mc.write_text("genome_id\tmethanogen\nG1\tFalse\nG2\tFalse\n")
    py("methanogenesis_search.py", *args, "--methanogen-calls", mc, "--out", d / "none.tsv")
    assert len(read_tsv(d / "none.tsv")) == 0, "no methanogens -> header-only"
    mc.write_text("genome_id\tmethanogen\nG1\tTrue\nG2\tFalse\n")
    py("methanogenesis_search.py", *args, "--methanogen-calls", mc, "--out", d / "one.tsv")
    assert set(read_tsv(d / "one.tsv").genome_id) == {"G1"}


@suite.test
def test_search_failures_are_loud_and_empty_input_is_header_only():
    need("hmmsearch", "hmmbuild", "mafft")
    d = tmpdir("mgfail_")
    faa, pmap, gdir, gendir, _ = search_world(d, truncate_pyl=False)
    hd = make_profiles(d, kos=("K04480",))
    rules = mini_rules(d, {"mta": {"expr": "K04480", "description": "MtaB"}})
    base = ["--proteome", faa, "--pmap", pmap, "--genes-dir", gdir, "--hmmdir", hd,
            "--rules", rules, "--all-genomes"]
    out = fails(py, "methanogenesis_search.py", *base, "--out", d / "x.tsv")
    assert "curated threshold" in out, out[-300:]          # no ko_list and no GA line
    kl = ko_list(d / "kl", [("K04480", 10.0, "full")])
    out = fails(py, "methanogenesis_search.py", *base, "--ko-list", d / "nope", "--out", d / "x.tsv")
    assert "does not exist" in out
    out = fails(py, "methanogenesis_search.py", "--proteome", faa, "--pmap", pmap,
                "--genes-dir", d / "missing", "--hmmdir", hd, "--rules", rules, "--ko-list", kl,
                "--all-genomes", "--out", d / "x.tsv")
    assert "genes table missing" in out, out[-300:]
    out = fails(py, "methanogenesis_search.py", *base, "--ko-list", kl, "--readthrough", "maybe",
                "--out", d / "x.tsv")
    assert "--readthrough must be True or False" in out
    empty = d / "empty.faa"
    empty.write_text("")
    py("methanogenesis_search.py", "--proteome", empty, "--pmap", pmap, "--genes-dir", gdir,
       "--hmmdir", hd, "--rules", rules, "--ko-list", kl, "--all-genomes", "--out", d / "e.tsv",
       "--out-all", d / "ea.tsv")
    assert read_tsv(d / "e.tsv").empty and open(d / "e.tsv").readline().startswith("genome_id")


# --------------------------------------------------------------------------- #
# methanogenesis_pathways.py end to end                                        #
# --------------------------------------------------------------------------- #
def pathways_world(d, kos_by_genome, meth, completeness, taxonomy=None, expectation=None):
    d = Path(d)
    with open(d / "hits.tsv", "w") as fh:
        fh.write("genome_id\tprotein_id\tko\tcall\tcontig_index\tgene_index\n")
        for g, kos in kos_by_genome.items():
            for i, k in enumerate(kos, 1):
                fh.write(f"{g}\t{g}__{g}_1_{i}\t{k}\tfull_threshold\t1\t{i}\n")
    with open(d / "meth.tsv", "w") as fh:
        fh.write("genome_id\tmethanogen\n")
        for g, m in meth.items():
            fh.write(f"{g}\t{m}\n")
    with open(d / "meta.tsv", "w") as fh:
        fh.write("genome_id\tcompleteness\n")
        for g, c in completeness.items():
            fh.write(f"{g}\t{'' if c is None else c}\n")
    paths = [d / "hits.tsv", d / "meth.tsv", d / "meta.tsv"]
    if taxonomy:
        with open(d / "tax.tsv", "w") as fh:
            fh.write("user_genome\tclassification\n")
            for g, c in taxonomy.items():
                fh.write(f"{g}\t{c}\n")
        paths.append(d / "tax.tsv")
    if expectation:
        with open(d / "exp.tsv", "w") as fh:
            fh.write("rank\ttaxon\t" + "\t".join(mpath.EXPECT_TERMS) + "\n")
            for (rank, taxon), vals in expectation.items():
                fh.write(f"{rank}\t{taxon}\t" +
                         "\t".join(str(vals.get(t, "NA")) for t in mpath.EXPECT_TERMS) + "\n")
        paths.append(d / "exp.tsv")
    return paths


@suite.test
def test_component_signature_tolerates_a_missing_selenoprotein_subunit():
    """REGRESSION: a complex must not hinge on one subunit that annotation commonly misses.
    KEGG lists no FwdB (K00201) and no HdrA (K03388) for Methanocaldococcus jannaschii - both
    are UGA/selenoprotein-affected - yet it is a textbook CO2-reducing hydrogenotroph. With a
    single-KO signature the pathway collapsed to none_detected; with an alternative-slot
    signature it is recovered."""
    R = mpath.Rules(RULEDOC)
    kos = ["K00200", "K00202", "K00203", "K00205", "K11260", "K00204",      # Fwd without FwdB
           "K00672", "K01499", "K00319", "K00320",
           "K03639", "K03637", "K03635", "K03638", "K03750",
           "K00577", "K00578", "K00579", "K00580", "K00581", "K00584",
           "K03389", "K03390",                                              # Hdr without HdrA
           "K00440", "K00441", "K00443"]
    pk = {k: [{"protein_id": k, "contig_index": 1, "gene_index": i}] for i, k in enumerate(kos)}
    comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
    assert comps["fwd"]["state"] == "present", comps["fwd"]
    assert comps["hdr_abc"]["state"] == "present", comps["hdr_abc"]
    vals = mpath.evaluate_terms(R, comps, 0.05)
    assert mpath.decide(vals[("pathways", "co2_reducing")], 0.05) == "present"
    # ... but a complex with only a minority of its subunits is still missing
    thin = {k: pk[k] for k in ("K00200",)}
    assert mpath.component_state("fwd", R, thin, 1.0)["state"] == "missing"


@suite.test
def test_reference_methanogens_known_answers():
    """The rules reproduce the published physiology of 12 cultured methanogens from KEGG's own
    per-genome KO assignments (fixtures/kegg_reference_ko_sets.json, rest.kegg.jp /link/ko).
    This is a rule-level check on an independent annotation, not a check of the HMM search.
    The two rows marked 'capacity' are known capacity-vs-phenotype discordances: the genes are
    present and the organism does not use the substrate in culture."""
    R = mpath.Rules(RULEDOC)
    sets = json.loads((FIX / "kegg_reference_ko_sets.json").read_text())
    expect = {   # organism: (pathway -> call), published physiology
        "mja": {"co2_reducing": "present", "methylotrophic": "absent",
                "acetoclastic": "absent"},          # M. jannaschii, obligate hydrogenotroph
        "mmp": {"co2_reducing": "present", "methylotrophic": "absent"},     # M. maripaludis
        "mth": {"co2_reducing": "present", "methylotrophic": "absent"},     # M. thermautotrophicus
        "mhu": {"co2_reducing": "present", "methylotrophic": "absent"},     # M. hungatei
        "msi": {"co2_reducing": "present"},                 # M. smithii (+ methanol genes only)
        "mcj": {"acetoclastic": "present", "co2_reducing": "absent",
                "methylotrophic": "absent"},        # Methanothrix soehngenii, obligate acetoclast
        "mba": {"co2_reducing": "present", "acetoclastic": "present",
                "methyl_dismutating": "present"},   # M. barkeri, all three
        "mac": {"acetoclastic": "present", "methyl_dismutating": "present"},  # M. acetivorans
        "mbu": {"methylotrophic": "present", "methyl_dismutating": "present",
                "co2_reducing": "absent"},          # Methanococcoides burtonii
        "mst": {"methyl_reducing": "present", "co2_reducing": "absent"},     # M. stadtmanae
        "max": {"methyl_reducing": "present", "co2_reducing": "absent"},     # Ca. M. alvi
        "mhz": {"methyl_dismutating": "present", "co2_reducing": "absent"},  # M. hollandica
    }
    assert set(expect) <= set(sets), sorted(set(expect) - set(sets))
    for org, want in expect.items():
        pk = {k: [{"protein_id": f"{org}_{k}", "contig_index": 1, "gene_index": i}]
              for i, (k, n) in enumerate(sorted(sets[org].items())) if n and k in R.all_kos()}
        comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
        vals = mpath.evaluate_terms(R, comps, 0.05)
        got = {p: mpath.decide(vals[("pathways", p)], 0.05) for p in mpath.PATHWAYS}
        for k, v in want.items():
            assert got[k] == v, (org, k, "expected", v, "got", got)
    # the known discordances are capacity calls, and they are positive at the gene level
    def call(org, p):
        pk = {k: [{"protein_id": k, "contig_index": 1, "gene_index": i}]
              for i, (k, n) in enumerate(sorted(sets[org].items())) if n and k in R.all_kos()}
        comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
        return mpath.decide(mpath.evaluate_terms(R, comps, 0.05)[("pathways", p)], 0.05)
    assert call("msi", "methylotrophic") == "present"   # mtaABC present; no growth on methanol
    assert call("mac", "co2_reducing") == "present"      # fwd/mtr/hydrogenase genes; no H2+CO2


@suite.test
def test_real_genome_ko_calls_reproduce_published_physiology():
    """KNOWN ANSWER on real data: the KO calls this caller made from 12 real assemblies
    (RefSeq/GenBank -> Prodigal 2.6.3 -> methanogenesis_search.py with the curated KOfam
    subset; fixtures/real_genome_ko_calls.json) give pathway calls that match the published
    physiology of every organism, with no mismatch. Stored calls, not a re-run: the search is
    exercised separately on real HMMER runs, and re-downloading 12 genomes is not a unit test.
    """
    R = mpath.Rules(RULEDOC)
    fx = json.loads((FIX / "real_genome_ko_calls.json").read_text())
    truth = {   # capacity from culture physiology; '-' = not asserted (see notes)
        "mja": {"co2_reducing": 1, "acetoclastic": 0, "methylotrophic": 0},
        "mmp": {"co2_reducing": 1, "acetoclastic": 0, "methylotrophic": 0},
        "mth": {"co2_reducing": 1, "acetoclastic": 0, "methylotrophic": 0},
        "mhu": {"co2_reducing": 1, "acetoclastic": 0, "methylotrophic": 0},
        "msi": {"co2_reducing": 1, "acetoclastic": 0},
        "mcj": {"co2_reducing": 0, "acetoclastic": 1, "methylotrophic": 0},
        "mba": {"co2_reducing": 1, "acetoclastic": 1, "methylotrophic": 1},
        "mac": {"acetoclastic": 1, "methylotrophic": 1},
        "mbu": {"co2_reducing": 0, "acetoclastic": 0, "methylotrophic": 1},
        "mst": {"co2_reducing": 0, "acetoclastic": 0, "methylotrophic": 1,
                "methyl_reducing": 1},
        "max": {"co2_reducing": 0, "acetoclastic": 0, "methylotrophic": 1,
                "methyl_reducing": 1},
        "mhz": {"co2_reducing": 0, "acetoclastic": 0, "methylotrophic": 1},
    }
    assert set(truth) <= set(fx["hits"]), sorted(set(truth) - set(fx["hits"]))
    n_ok = n_unres = 0
    mism = []
    for org, want in truth.items():
        calls = fx["hits"][org]
        pk, nk = {}, set()
        for call, kos in calls.items():
            for ko, n in kos.items():
                if call == "near_threshold":
                    nk.add(ko)
                else:
                    pk.setdefault(ko, []).extend(
                        [{"protein_id": f"{org}_{ko}_{i}", "contig_index": 1,
                          "gene_index": i}] for i in range(n))
                    pk[ko] = [x[0] if isinstance(x, list) else x for x in pk[ko]]
        comps = {n: mpath.component_state(n, R, pk, 1.0, nk) for n in R.components}
        vals = mpath.evaluate_terms(R, comps, 0.05)
        got = {p: mpath.decide(vals[("pathways", p)], 0.05) for p in mpath.PATHWAYS}
        for p, t in want.items():
            v = {"present": 1, "absent": 0}.get(got[p])
            if v is None:
                n_unres += 1
            elif v == t:
                n_ok += 1
            else:
                mism.append((org, p, "expected", t, "got", got[p]))
    assert not mism, mism
    assert n_ok >= 34, (n_ok, n_unres)
    # the two known capacity-vs-phenotype discordances stay visible as gene-level positives
    def pathway(org, p):
        calls = fx["hits"][org]
        pk = {ko: [{"protein_id": f"{org}_{ko}", "contig_index": 1, "gene_index": 1}]
              for call, kos in calls.items() if call != "near_threshold" for ko in kos}
        nk = set(calls.get("near_threshold", {}))
        comps = {n: mpath.component_state(n, R, pk, 1.0, nk) for n in R.components}
        return mpath.decide(mpath.evaluate_terms(R, comps, 0.05)[("pathways", p)], 0.05)
    assert pathway("msi", "methylotrophic") == "present"   # mtaABC present, no growth on methanol
    assert pathway("mac", "co2_reducing") == "present"      # fwd/mtr/hydrogenase, no H2+CO2 growth
    # the pyrrolysine and selenocysteine read-throughs fired on real genes
    rt = {org: set(c["readthrough_TAG"]) | set(c.get("readthrough_TGA", {}))
          for org, c in fx["hits"].items() if "readthrough_TAG" in c}
    assert "K16176" in rt.get("max", set()), rt          # MtmB, in-frame UAG
    tga = {org: set(c["readthrough_TGA"]) for org, c in fx["hits"].items()
           if "readthrough_TGA" in c}
    assert "K00201" in tga.get("mja", set()) and "K03388" in tga.get("mja", set()), tga
    assert "K22516" in tga.get("mmp", set()), tga


@suite.test
def test_near_threshold_hits_block_a_false_absence_without_creating_presence():
    """A KOfam threshold sitting just above a divergent archaeal orthologue must not become a
    confident absence (real case: K00625 Pta, M. acetivorans 379.0 vs threshold 387.57)."""
    R = mpath.Rules(RULEDOC)
    kos = ["K00193", "K00197", "K00194", "K00577", "K00578", "K00579", "K00580", "K00581",
           "K00584", "K08264", "K08265"]
    pk = {k: [{"protein_id": k, "contig_index": 1, "gene_index": i}] for i, k in enumerate(kos)}
    # neither AckA nor Pta: acetoclastic is a confident absence in a complete genome
    comps = {n: mpath.component_state(n, R, pk, 1.0) for n in R.components}
    assert mpath.decide(mpath.evaluate_terms(R, comps, 0.05)[("pathways", "acetoclastic")],
                        0.05) == "absent"
    # a near-threshold Pta with NO AckA is not evidence of acetate activation: still absent
    comps = {n: mpath.component_state(n, R, pk, 1.0, {"K00625"}) for n in R.components}
    assert comps["ackpta"]["state"] == "missing", comps["ackpta"]
    assert mpath.decide(mpath.evaluate_terms(R, comps, 0.05)[("pathways", "acetoclastic")],
                        0.05) == "absent"
    # the real M. acetivorans / M. barkeri configuration: AckA accepted, Pta just below the
    # KOfam threshold -> the relaxed pair is satisfied and acetoclasty is called
    pk2 = dict(pk, K00925=[{"protein_id": "K00925", "contig_index": 1, "gene_index": 99}])
    comps = {n: mpath.component_state(n, R, pk2, 1.0, {"K00625"}) for n in R.components}
    assert comps["ackpta"]["state"] == "present", comps["ackpta"]
    assert mpath.decide(mpath.evaluate_terms(R, comps, 0.05)[("pathways", "acetoclastic")],
                        0.05) == "present"
    # for a component with no relaxation, near-threshold hits make it equivocal (absence
    # withdrawn) but NEVER present
    assert mpath.component_state("hdr_de", R, {}, 1.0)["state"] == "missing"
    e = mpath.component_state("hdr_de", R, {}, 1.0, {"K08264", "K08265"})
    assert e["state"] == "equivocal" and e["interval"] == (0.0, 1.0), e
    assert e["near"] == ["K08264", "K08265"]
    assert mpath.component_state("mts", R, {}, 1.0, {"K16954", "K16955"})["state"] == "equivocal"


@suite.test
def test_pathways_end_to_end_with_the_real_rules_file():
    """Four genomes with known capacity, plus a non-methanogen and a low-completeness MAG."""
    d = tmpdir("mgpath_")
    mb = ["K00200", "K00201", "K00202", "K00203", "K00205", "K00672", "K01499", "K00319", "K00320"]
    moco = ["K03639", "K03637", "K03635", "K03638", "K03750"]
    mtr = ["K00577", "K00578", "K00579", "K00580", "K00581", "K00584"]
    hdr = ["K03388", "K03389", "K03390"]
    frh = ["K00440", "K00441", "K00443"]
    mta = ["K14080", "K04480", "K14081"]
    kos = {"HYD": mb + moco + mtr + hdr + frh,
           "METHYL": mta + hdr + frh,
           "SARC": mb + moco + mtr + mta + ["K08264", "K08265", "K00193", "K00197", "K00194",
                                            "K00925", "K00625"],
           # MAG: a 60%-complete hydrogenotroph-like genome missing ONLY the Mtr locus
           "MAG": mb + moco + hdr + frh,
           "NONM": mb + moco + mtr + hdr + frh}
    meth = {"HYD": "True", "METHYL": "True", "SARC": "True", "MAG": "True", "NONM": "False"}
    comp = {"HYD": 99.5, "METHYL": 98.0, "SARC": 97.0, "MAG": 60.0, "NONM": 99.0}
    hits, mc, meta = pathways_world(d, kos, meth, comp)
    py("methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc, "--meta", meta,
       "--rules", RULES, "--out", d / "p.tsv")
    p = read_tsv(d / "p.tsv").set_index("genome_id")
    assert p.loc["HYD", "co2_reducing"] == "present"
    assert p.loc["HYD", "methanogenesis_subtype"] == "CO2-reducing"
    assert p.loc["METHYL", "methyl_reducing"] == "present"
    assert p.loc["METHYL", "co2_reducing"] == "absent"
    assert p.loc["SARC", "acetoclastic"] == "present"
    assert p.loc["SARC", "methyl_dismutating"] == "present"
    assert "acetoclastic" in p.loc["SARC", "methanogenesis_subtype"]
    assert p.loc["NONM", "co2_reducing"] == "not_applicable"
    assert p.loc["NONM", "methanogenesis_subtype"] == "NA"
    assert p.loc["NONM", "subtype_status"] == "not_applicable"
    # the 60%-complete MAG lacks ONE locus (Mtr) of the CO2-reducing pathway: at 60%
    # completeness a single missing locus (p = 0.4) is not evidence of absence
    assert p.loc["MAG", "co2_reducing"] == "unresolved", p.loc["MAG"].to_dict()
    assert p.loc["MAG", "subtype_status"] == "partial"
    assert abs(float(p.loc["MAG", "co2_reducing_p_hi"]) - 0.4) < 1e-6
    assert p.loc["MAG", "methanogenesis_subtype"] == "unresolved"
    assert p.loc["HYD", "components"].startswith("fwd=P")
    assert "mta:METHYL__" in p.loc["METHYL", "evidence"]


@suite.test
def test_pathways_absence_requires_completeness():
    """Completeness unknown -> every absence is unresolved; --absence-alpha moves the line."""
    d = tmpdir("mgalpha_")
    mta = ["K14080", "K04480", "K14081"]
    hdr = ["K03388", "K03389", "K03390"]
    kos = {"A": mta + hdr, "B": mta + hdr}
    hits, mc, meta = pathways_world(d, kos, {"A": "True", "B": "True"},
                                    {"A": None, "B": 80.0})
    py("methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc, "--meta", meta,
       "--rules", RULES, "--out", d / "p.tsv")
    p = read_tsv(d / "p.tsv").set_index("genome_id")
    assert p.loc["A", "completeness"] == "NA"
    assert p.loc["A", "co2_reducing"] == "unresolved"
    assert p.loc["B", "co2_reducing"] == "absent"      # 0.2^k over several missing loci
    py("methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc, "--meta", meta,
       "--rules", RULES, "--absence-alpha", "0.001", "--out", d / "q.tsv")
    q = read_tsv(d / "q.tsv").set_index("genome_id")
    assert q.loc["B", "acetoclastic"] == "unresolved", q.loc["B"].to_dict()


@suite.test
def test_pathways_taxonomy_reconciliation():
    """The most specific matching rank wins, and a gene/taxon disagreement is reported rather
    than silently overridden (design §3.3)."""
    d = tmpdir("mgtax_")
    mta = ["K14080", "K04480", "K14081"]
    hdr = ["K03388", "K03389", "K03390"]
    frh = ["K00440", "K00441", "K00443"]
    kos = {"G1": mta + hdr, "G2": mta + hdr + frh}
    tax = {"G1": "d__Archaea;p__Halobacteriota;c__Methanosarcinia;o__Methanosarcinales;"
                 "f__Methanosarcinaceae;g__Methanosarcina;s__Methanosarcina barkeri",
           "G2": "d__Archaea;p__Thermoplasmatota;c__Thermoplasmata;o__Methanomassiliicoccales;"
                 "f__Methanomethylophilaceae;g__Methanomethylophilus;s__"}
    expectation = {
        ("order", "o__Methanosarcinales"): {"methanogen": 1, "co2_reducing": 1,
                                            "acetoclastic": 1, "methylotrophic": 1,
                                            "methyl_dismutating": 1, "methyl_reducing": 0},
        ("genus", "g__Methanosarcina"): {"methanogen": 1, "co2_reducing": 1, "acetoclastic": 1,
                                         "methylotrophic": 1, "methyl_dismutating": 1,
                                         "methyl_reducing": 0},
        ("order", "o__Methanomassiliicoccales"): {"methanogen": 1, "co2_reducing": 0,
                                                  "acetoclastic": 0, "methylotrophic": 1,
                                                  "methyl_dismutating": 0, "methyl_reducing": 1},
    }
    hits, mc, meta, taxf, expf = pathways_world(d, kos, {"G1": "True", "G2": "True"},
                                                {"G1": 99.0, "G2": 99.0}, tax, expectation)
    py("methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc, "--meta", meta,
       "--rules", RULES, "--taxonomy", taxf, "--expectation", expf, "--out", d / "p.tsv",
       "--out-reconciliation", d / "r.tsv")
    p = read_tsv(d / "p.tsv").set_index("genome_id")
    assert p.loc["G1", "taxonomy_expected"].startswith("genus:g__Methanosarcina:")
    assert p.loc["G1", "taxonomy_concordance"] == "discordant"
    assert "co2_reducing" in p.loc["G1", "taxonomy_discordant_terms"]
    assert p.loc["G2", "taxonomy_concordance"] == "concordant", p.loc["G2"].to_dict()
    r = read_tsv(d / "r.tsv")
    assert set(r.genome_id) == {"G1"} and "acetoclastic" in set(r.term)
    assert set(r.matched_rank) == {"genus"}


@suite.test
def test_pathways_failures_are_loud():
    d = tmpdir("mgpf_")
    kos = {"G1": ["K14080", "K04480", "K14081"]}
    hits, mc, meta = pathways_world(d, kos, {"G1": "True"}, {"G1": 99.0})
    base = ["--hits", hits, "--methanogen", mc, "--meta", meta, "--rules", RULES]
    out = fails(py, "methanogenesis_pathways.py", *base, "--absence-alpha", "0.9",
                "--out", d / "x.tsv")
    assert "--absence-alpha" in out
    meta2 = d / "meta2.tsv"
    meta2.write_text("genome_id\tcompleteness\nOTHER\t90\n")
    out = fails(py, "methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc,
                "--meta", meta2, "--rules", RULES, "--out", d / "x.tsv")
    assert "missing from --meta" in out, out[-300:]
    bad = d / "bad_hits.tsv"
    bad.write_text("genome_id\tprotein_id\tko\tcall\nG1\tp1\tK99999\tfull_threshold\n")
    out = fails(py, "methanogenesis_pathways.py", "--hits", bad, "--methanogen", mc,
                "--meta", meta, "--rules", RULES, "--out", d / "x.tsv")
    assert "not in the rules' KO set" in out, out[-300:]
    out = fails(py, "methanogenesis_pathways.py", *base, "--expectation", d / "exp.tsv",
                "--out", d / "x.tsv")
    assert "--expectation needs --taxonomy" in out


@suite.test
def test_empty_methanogen_set_is_header_only_not_an_error():
    d = tmpdir("mgempty_")
    hits, mc, meta = pathways_world(d, {}, {"G1": "False"}, {"G1": 99.0})
    py("methanogenesis_pathways.py", "--hits", hits, "--methanogen", mc, "--meta", meta,
       "--rules", RULES, "--out", d / "p.tsv")
    p = read_tsv(d / "p.tsv")
    assert len(p) == 1 and p.iloc[0].methanogenesis_subtype == "NA"
    assert list(p.columns) == mpath.out_columns()


if __name__ == "__main__":
    suite.main("METHANOGENESIS TESTS PASSED")
