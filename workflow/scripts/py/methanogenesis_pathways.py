"""Call methanogenesis pathway capacity per genome (design §3.3) from the accepted KO hits.

Scope: genomes with methanogen == True in methanogen_calls.tsv (credible Mcr placed with the
methanogen references). ANME, divergent-Mcr and Mcr-less genomes are not_applicable: pathway
capacity is a property of methanogens only.

Components (rules file `components`) are KEGG expressions over the accepted KOs
(methanogenesis_ko_hits.tsv). A component is PRESENT when its expression is satisfied (KEGG rule:
every mandatory subunit) or, for a `relax` complex, when every signature KO is detected and at
least min_fraction of the mandatory subunit slots are. Otherwise it is MISSING with k missing
loci: 1 for a complex / single step, else the number of unsatisfied top-level steps.

Absence is an inference about an incomplete assembly, so each component carries an interval
[lo, hi] for the probability that it is truly encoded:
    present                       -> [1, 1]
    equivocal (a near-threshold   -> [0, 1]  (a curated threshold just above a divergent
      hit would complete it)                  archaeal orthologue is not an absence)
    missing, completeness c known -> [p, p], p = (1 - c)^k   (each missing locus independently
                                     absent from the assembly with probability 1 - c)
    missing, completeness unknown -> [0, 1]
    one_sided component missing   -> [0, 0]  (reported 'not_detected', never 'absent')
Rules combine intervals with all() = product (independent loci), any() = maximum (one
system provides the function; a missing function means its one locus was missed),
not() = [1 - hi, 1 - lo] and maybe() = [0, hi] (the genes allow the function but cannot show
it, e.g. an anabolic enzyme that may also be catabolic). A pathway / substrate is
    present    when lo >= 1 - alpha,
    absent     when hi <= alpha,
    unresolved otherwise                         (--absence-alpha, default 0.05).
So in a 100%-complete genome every call is resolved, while in a 70%-complete MAG one missing
locus (p = 0.3) leaves the call unresolved and four missing loci (p = 0.0081) resolve it as
absent. The bound assumes loci go missing independently and that detected genes are real;
CheckM2 completeness is an estimate (see RULES.md).

Outputs: --out methanogenesis_pathways.tsv (PLAN contract columns) and --out-reconciliation:
one row per (genome, term) where the gene-based call disagrees with the GTDB-taxon expectation
(--expectation, most specific matching rank: genus > family > order > class), for manual
review (design §3.3 "reconcile disagreements manually").
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils
import kegg_module

import argparse
import csv
import re

PATHWAYS = ["co2_reducing", "acetoclastic", "methylotrophic", "methyl_dismutating",
            "methyl_reducing"]
SUBSTRATES_2 = ["sub_h2", "sub_formate", "sub_methanol", "sub_methylamine"]
SUBSTRATES_1 = ["sub_methylthiol", "sub_quatamine", "sub_methoxyl"]
EXPECT_TERMS = ["methanogen"] + PATHWAYS
RANKS = [("g__", "genus"), ("f__", "family"), ("o__", "order"), ("c__", "class")]


def out_columns():
    cols = ["genome_id", "methanogen", "completeness"]
    for p in PATHWAYS:
        cols += [p, f"{p}_p_lo", f"{p}_p_hi"]
    cols += SUBSTRATES_2 + SUBSTRATES_1
    cols += ["methanogenesis_subtype", "subtype_status", "components", "evidence",
             "taxonomy_expected", "taxonomy_concordance", "taxonomy_discordant_terms"]
    return cols


RECON_COLS = ["genome_id", "classification", "matched_rank", "matched_taxon", "term",
              "expected", "called", "note"]


# ---- rule expressions -------------------------------------------------------------------
_TOK = re.compile(r"\s*(all|any|not|maybe|\(|\)|,|[A-Za-z_][A-Za-z0-9_]*)")


def parse_rule(text: str):
    """'all(a, any(b, c), not(d))' -> nested tuples ('all', [...]) / ('not', x) / name."""
    toks, pos, s = [], 0, str(text)
    while pos < len(s):
        m = _TOK.match(s, pos)
        if not m:
            if s[pos:].strip() == "":
                break
            raise ValueError(f"cannot parse rule {text!r} at {pos}")
        toks.append(m.group(1))
        pos = m.end()
    i = 0

    def expr():
        nonlocal i
        if i >= len(toks):
            raise ValueError(f"unexpected end of rule {text!r}")
        t = toks[i]
        i += 1
        if t in ("all", "any", "not", "maybe"):
            if i >= len(toks) or toks[i] != "(":
                raise ValueError(f"{t} needs '(' in {text!r}")
            i += 1
            args = [expr()]
            while i < len(toks) and toks[i] == ",":
                i += 1
                args.append(expr())
            if i >= len(toks) or toks[i] != ")":
                raise ValueError(f"missing ')' in {text!r}")
            i += 1
            if t in ("not", "maybe"):
                if len(args) != 1:
                    raise ValueError(f"{t}() takes one argument in {text!r}")
                return (t, args[0])
            return (t, args)
        if t in ("(", ")", ","):
            raise ValueError(f"unexpected {t!r} in {text!r}")
        return t

    node = expr()
    if i != len(toks):
        raise ValueError(f"trailing tokens in rule {text!r}")
    return node


def rule_names(node) -> set:
    if isinstance(node, str):
        return {node}
    if node[0] in ("not", "maybe"):
        return rule_names(node[1])
    return set().union(*(rule_names(a) for a in node[1]))


def signature_slots(relax) -> list:
    """`signature` entries: a KO string = that KO must be detected; a list of KOs = at least one
    of them must be. A signature is a diagnostic core of the complex, NOT a single point of
    failure: a lone subunit that is commonly missed (a selenoprotein read through a UGA, a
    lineage-specific paralogue) must never be the only signature.
    """
    out = []
    for item in (relax.get("signature") or []):
        out.append({item} if isinstance(item, str) else set(item))
    return out


def decide(iv, alpha):
    lo, hi = iv
    if lo >= 1 - alpha:
        return "present"
    if hi <= alpha:
        return "absent"
    return "unresolved"


class Rules:
    def __init__(self, rules: dict):
        self.raw = rules
        self.components = rules.get("components") or {}
        if not self.components:
            raise ValueError("rules file defines no components")
        self.parsed = {}
        for name, c in self.components.items():
            node = kegg_module.parse(c["expr"])
            relax = c.get("relax")
            if relax:
                mand, _ = kegg_module.complex_slots(node)
                slot_kos = set().union(*map(set, mand))
                for sig in signature_slots(relax):
                    if not sig <= slot_kos:
                        raise ValueError(f"component {name}: signature {sorted(sig - slot_kos)} "
                                         "is not a mandatory subunit")
                mf = float(relax.get("min_fraction", 1.0))
                if not (0 < mf <= 1):
                    raise ValueError(f"component {name}: min_fraction must be in (0, 1]")
            self.parsed[name] = node
        self.terms, seen = {}, set(self.components)
        for block in ("derived", "pathways", "substrates"):
            for name, text in (rules.get(block) or {}).items():
                if name in seen:
                    raise ValueError(f"term {name} is defined twice (components / derived / "
                                     "pathways / substrates share one namespace)")
                seen.add(name)
                self.terms[(block, name)] = parse_rule(text)
        for p in PATHWAYS:
            if ("pathways", p) not in self.terms:
                raise ValueError(f"rules file lacks pathway {p}")
        for s in SUBSTRATES_2 + SUBSTRATES_1:
            if ("substrates", s) not in self.terms:
                raise ValueError(f"rules file lacks substrate flag {s}")
        known = set(self.components) | {n for (_, n) in self.terms}
        for key, node in self.terms.items():
            unknown = rule_names(node) - known
            if unknown:
                raise ValueError(f"rule {key[1]} references unknown term(s) {sorted(unknown)}")
        for s in SUBSTRATES_1:
            node = self.terms[("substrates", s)]
            if not (isinstance(node, str) and self.components.get(node, {}).get("one_sided")):
                raise ValueError(f"one-sided flag {s} must name a one_sided component")
        self.order = list(rules.get("subtype_order") or ["co2_reducing", "acetoclastic",
                                                          "methyl_dismutating", "methyl_reducing"])
        self.labels = dict(rules.get("subtype_labels") or {})

    def all_kos(self) -> set:
        return set().union(*(kegg_module.kos(n) for n in self.parsed.values()))


def component_state(name, rules: Rules, present_kos: dict, c, near_kos=None):
    """State dict of one component for one genome (present_kos: KO -> [protein records]).

    near_kos = KOs with a near-threshold hit (methanogenesis_search call = near_threshold).
    They never make a component PRESENT, but a component that would be present with them is
    'equivocal' (interval [0, 1] -> unresolved) instead of absent: a curated threshold sitting
    just above a divergent archaeal orthologue must not become a confident absence.
    """
    spec = rules.components[name]
    node = rules.parsed[name]
    P = set(present_kos)
    near = set(near_kos or ()) - P
    one_sided = bool(spec.get("one_sided"))
    relax = spec.get("relax")
    if relax:
        mand, _ = kegg_module.complex_slots(node)
        ok_slots = [s for s in mand if P & set(s)]
        sigs = signature_slots(relax)
        frac = len(ok_slots) / len(mand)
        is_present = (all(s & P for s in sigs)
                      and frac + 1e-12 >= float(relax.get("min_fraction", 1.0))) \
            or kegg_module.satisfied(node, P)
        detail = f"{len(ok_slots)}/{len(mand)}"
        loci = 1
    else:
        ev = kegg_module.evaluate(node, P)
        is_present = ev["complete"]
        if len(ev["steps"]) > 1:
            detail = f"{ev['n_satisfied']}/{ev['n_steps']}"
            loci = ev["n_steps"] - ev["n_satisfied"]
        else:
            try:
                mand, _ = kegg_module.complex_slots(node)
            except ValueError:
                mand = []
            detail = f"{sum(1 for s in mand if P & set(s))}/{len(mand)}" if mand else ""
            loci = 1
    kos_used = sorted(kegg_module.kos(node) & P)
    proteins = sorted({r["protein_id"] for k in kos_used for r in present_kos[k]})
    equivocal = False
    if not is_present and near:
        PN = P | near
        if relax:
            mand2, _ = kegg_module.complex_slots(node)
            ok2 = [s for s in mand2 if PN & set(s)]
            equivocal = (all(sg & PN for sg in signature_slots(relax))
                         and len(ok2) / len(mand2) + 1e-12 >= float(relax.get("min_fraction", 1.0))
                         ) or kegg_module.satisfied(node, PN)
        else:
            equivocal = kegg_module.satisfied(node, PN)
    if is_present:
        iv, state = (1.0, 1.0), "present"
    elif equivocal:
        iv, state = (0.0, 1.0), "equivocal"
    elif one_sided:
        iv, state = (0.0, 0.0), "not_detected"
    elif c is None:
        iv, state = (0.0, 1.0), "missing"
    else:
        p = (1.0 - c) ** max(1, loci)
        iv, state = (p, p), "missing"
    coloc = ""
    if is_present and len(kos_used) >= 2:
        coloc = colocalisation(kos_used, present_kos)
    return {"state": state, "interval": iv, "detail": detail, "proteins": proteins,
            "coloc": coloc, "loci_missing": 0 if is_present else max(1, loci),
            "near": sorted(near & kegg_module.kos(node))}


def colocalisation(kos_used, present_kos, window=20) -> str:
    """'coloc' when one contig window of `window` genes holds a copy of every detected subunit,
    else 'k/n' = most subunits found together."""
    pos = {}
    for k in kos_used:
        for r in present_kos[k]:
            if r.get("contig_index") not in (None, "") and r.get("gene_index") not in (None, ""):
                pos.setdefault(k, []).append((int(r["contig_index"]), int(r["gene_index"])))
    if len(pos) < 2:
        return ""
    best = 0
    anchors = [p for v in pos.values() for p in v]
    for ci, gi in anchors:
        n = sum(1 for v in pos.values() if any(c2 == ci and abs(g2 - gi) <= window for c2, g2 in v))
        best = max(best, n)
    return "coloc" if best == len(pos) else f"{best}/{len(pos)}"


def evaluate_terms(rules: Rules, comps: dict, alpha: float) -> dict:
    """Intervals for every derived term, pathway and substrate flag (memoised, cycle-checked)."""
    cache, stack = {}, set()
    by_name = {}
    for (block, name), node in rules.terms.items():
        by_name.setdefault(name, []).append((block, node))

    def value(name):
        if name in comps:
            return comps[name]["interval"]
        if name in cache:
            return cache[name]
        if name in stack:
            raise ValueError(f"rule cycle through {name}")
        defs = by_name.get(name)
        if not defs:
            raise ValueError(f"unknown term {name}")
        stack.add(name)
        node = [n for b, n in defs if b != "substrates"] or [n for _, n in defs]
        iv = ev(node[0])
        stack.discard(name)
        cache[name] = iv
        return iv

    def ev(node):
        if isinstance(node, str):
            return value(node)
        op, arg = node
        if op == "not":
            lo, hi = ev(arg)
            return (1.0 - hi, 1.0 - lo)
        if op == "maybe":                      # possible but not inferable: [0, hi]
            return (0.0, ev(arg)[1])
        ivs = [ev(a) for a in arg]
        if op == "all":
            lo = hi = 1.0
            for a, b in ivs:
                lo *= a
                hi *= b
            return (lo, hi)
        return (max(a for a, _ in ivs), max(b for _, b in ivs))

    out = {}
    for (block, name), node in rules.terms.items():
        out[(block, name)] = ev(node) if block == "substrates" else value(name)
    return out


def subtype_label(calls: dict, rules: Rules) -> tuple:
    labels = rules.labels
    present = [p for p in rules.order if calls[p] == "present"]
    if not present and calls.get("methylotrophic") == "present":
        present = ["methylotrophic"]
    resolved = sum(calls[p] in ("present", "absent") for p in PATHWAYS)
    status = "resolved" if resolved == len(PATHWAYS) else ("partial" if resolved else "unresolved")
    if present:
        return "+".join(labels.get(p, p) for p in present), status
    if any(calls[p] == "unresolved" for p in PATHWAYS):
        return "unresolved", status
    return "none_detected", status


# ---- inputs -------------------------------------------------------------------------------
def read_tsv(path, required, what):
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None or not set(required) <= set(rd.fieldnames):
            raise ValueError(f"{what} {path} lacks column(s) "
                             f"{sorted(set(required) - set(rd.fieldnames or []))}")
        return list(rd)


def as_bool(x) -> bool:
    return str(x).strip().lower() in ("true", "1", "yes")


def completeness_of(meta_rows, col):
    out = {}
    for r in meta_rows:
        v = str(r.get(col, "")).strip()
        try:
            out[r["genome_id"].strip()] = None if v in ("", "NA", "nan") else float(v)
        except ValueError:
            raise ValueError(f"completeness '{v}' of {r['genome_id']} is not a number")
    return out


def load_expectation(path):
    if not path:
        return {}
    rows = read_tsv(path, ["rank", "taxon"] + EXPECT_TERMS, "taxonomy expectation")
    out = {}
    for r in rows:
        rank, taxon = r["rank"].strip(), r["taxon"].strip()
        if rank not in {x for _, x in RANKS}:
            raise ValueError(f"taxonomy expectation rank {rank!r} must be one of genus, family, "
                             "order, class")
        vals = {}
        for t in EXPECT_TERMS:
            v = str(r[t]).strip()
            if v not in ("1", "0", "NA", ""):
                raise ValueError(f"taxonomy expectation {taxon} {t}={v!r} (expected 1, 0 or NA)")
            vals[t] = None if v in ("NA", "") else int(v)
        if (rank, taxon) in out:
            raise ValueError(f"taxonomy expectation lists {rank} {taxon} twice")
        out[(rank, taxon)] = vals
    return out


def match_expectation(classification, table):
    if not table or not classification:
        return None, None, None
    toks = [t.strip() for t in classification.split(";")]
    for prefix, rank in RANKS:
        for t in toks:
            if t.startswith(prefix) and len(t) > len(prefix) and (rank, t) in table:
                return rank, t, table[(rank, t)]
    return None, None, None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="methanogenesis_ko_hits.tsv")
    ap.add_argument("--methanogen", required=True, help="methanogen_calls.tsv")
    ap.add_argument("--meta", required=True, help="genome metadata with completeness")
    ap.add_argument("--completeness-col", default="completeness")
    ap.add_argument("--rules", required=True)
    ap.add_argument("--taxonomy", default="", help="taxonomy.tsv (user_genome, classification)")
    ap.add_argument("--expectation", default="", help="taxonomy_expectation.tsv")
    ap.add_argument("--absence-alpha", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-reconciliation", default="")
    args = ap.parse_args(argv)
    if not (0 < args.absence_alpha < 0.5):
        raise ValueError("--absence-alpha must be in (0, 0.5)")
    if args.expectation and not args.taxonomy:
        raise ValueError("--expectation needs --taxonomy")
    for p in (args.out, args.out_reconciliation):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    import methanogenesis_search as msearch
    rules = Rules(msearch.load_yaml(args.rules))
    universe_kos = rules.all_kos()
    mrows = read_tsv(args.methanogen, ["genome_id", "methanogen"], "methanogen calls")
    meta = completeness_of(read_tsv(args.meta, ["genome_id", args.completeness_col],
                                    "genome metadata"), args.completeness_col)
    hits = read_tsv(args.hits, ["genome_id", "protein_id", "ko", "call"], "KO hits")
    by_genome, near_by_genome = {}, {}
    for h in hits:
        ko = h["ko"].strip()
        if ko not in universe_kos:
            raise ValueError(f"KO hit {ko} ({h['protein_id']}) is not in the rules' KO set: the "
                             "hits table was made with a different rules file")
        gid = h["genome_id"].strip()
        if h["call"].strip() == "near_threshold":
            near_by_genome.setdefault(gid, set()).add(ko)
        else:
            by_genome.setdefault(gid, {}).setdefault(ko, []).append(h)
    taxonomy = {}
    if args.taxonomy:
        for r in read_tsv(args.taxonomy, ["user_genome", "classification"], "taxonomy"):
            taxonomy[r["user_genome"].strip()] = r["classification"].strip()
    expectation = load_expectation(args.expectation)

    cols = out_columns()
    out_rows, recon = [], []
    for r in mrows:
        gid = r["genome_id"].strip()
        is_m = as_bool(r["methanogen"])
        comp = meta.get(gid)
        if gid not in meta:
            raise ValueError(f"genome {gid} of the methanogen calls is missing from --meta")
        row = {c: "" for c in cols}
        row.update(genome_id=gid, methanogen=str(is_m),
                   completeness="NA" if comp is None else f"{comp:g}")
        calls = {}
        if not is_m:
            for p in PATHWAYS:
                row[p], row[f"{p}_p_lo"], row[f"{p}_p_hi"] = "not_applicable", "NA", "NA"
                calls[p] = "not_applicable"
            for s in SUBSTRATES_2 + SUBSTRATES_1:
                row[s] = "not_applicable"
            row.update(methanogenesis_subtype="NA", subtype_status="not_applicable",
                       components="", evidence="")
        else:
            c = None if comp is None else min(max(comp / 100.0, 0.0), 1.0)
            pk = by_genome.get(gid, {})
            nk = near_by_genome.get(gid, set())
            comps = {n: component_state(n, rules, pk, c, nk) for n in rules.components}
            vals = evaluate_terms(rules, comps, args.absence_alpha)
            for p in PATHWAYS:
                iv = vals[("pathways", p)]
                calls[p] = decide(iv, args.absence_alpha)
                row[p] = calls[p]
                row[f"{p}_p_lo"], row[f"{p}_p_hi"] = f"{iv[0]:.4g}", f"{iv[1]:.4g}"
            for s in SUBSTRATES_2:
                row[s] = decide(vals[("substrates", s)], args.absence_alpha)
            for s in SUBSTRATES_1:
                node = rules.terms[("substrates", s)]
                row[s] = "present" if comps[node]["state"] == "present" else "not_detected"
            row["methanogenesis_subtype"], row["subtype_status"] = subtype_label(calls, rules)
            parts, ev = [], []
            for n in rules.components:
                st = comps[n]
                code = {"present": "P", "missing": "A", "not_detected": "N",
                        "equivocal": "E"}[st["state"]]
                extra = ",".join(x for x in (st["detail"], st["coloc"],
                                             "near:" + "+".join(st["near"]) if st["near"] else "")
                                 if x)
                parts.append(f"{n}={code}" + (f"[{extra}]" if extra else ""))
                if st["state"] == "present":
                    ev.append(f"{n}:{'|'.join(st['proteins'])}")
            row["components"] = ";".join(parts)
            row["evidence"] = ";".join(ev)
        # ---- taxonomy reconciliation ----------------------------------------------------
        cls = taxonomy.get(gid, "")
        rank, taxon, exp = match_expectation(cls, expectation)
        if not args.expectation:
            row.update(taxonomy_expected="", taxonomy_concordance="NA", taxonomy_discordant_terms="")
        elif not cls:
            row.update(taxonomy_expected="", taxonomy_concordance="NA", taxonomy_discordant_terms="")
        elif exp is None:
            row.update(taxonomy_expected="", taxonomy_concordance="no_expectation",
                       taxonomy_discordant_terms="")
        else:
            called = {"methanogen": 1 if is_m else 0}
            for p in PATHWAYS:
                called[p] = {"present": 1, "absent": 0}.get(calls[p])
            disc, comparable = [], 0
            for t in EXPECT_TERMS:
                e, v = exp[t], called.get(t)
                if e is None or v is None:
                    continue
                comparable += 1
                if e != v:
                    disc.append(t)
                    recon.append({"genome_id": gid, "classification": cls, "matched_rank": rank,
                                  "matched_taxon": taxon, "term": t, "expected": e, "called": v,
                                  "note": "gene-based call disagrees with the taxon expectation"})
            row["taxonomy_expected"] = f"{rank}:{taxon}:" + ",".join(
                f"{t}={'NA' if exp[t] is None else exp[t]}" for t in EXPECT_TERMS)
            row["taxonomy_concordance"] = ("discordant" if disc else
                                           "concordant" if comparable else "no_expectation")
            row["taxonomy_discordant_terms"] = ",".join(disc)
        out_rows.append(row)

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in out_rows:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")
    if args.out_reconciliation:
        with open(args.out_reconciliation, "w") as out:
            out.write("\t".join(RECON_COLS) + "\n")
            for r in recon:
                out.write("\t".join(str(r[c]) for c in RECON_COLS) + "\n")
    n_m = sum(1 for r in out_rows if r["methanogen"] == "True")
    summary = {p: sum(1 for r in out_rows if r[p] == "present") for p in PATHWAYS}
    utils.log(f"methanogenesis_pathways: {len(out_rows)} genomes, {n_m} methanogens; present "
              f"calls {summary}; {len(recon)} taxonomy discordance(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
