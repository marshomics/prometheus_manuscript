"""Integrate detection, domain architecture, phylogenetic placement, catalytic-motif
evidence, defense-island context, length and fragmentation into a confidence-tiered
Prometheus catalog (Stage 1.2-1.8 of the design).

Scientific framing
-------------------
A DefenseFinder Prometheus call is *not* accepted at face value. ProA is a single
~2,221-aa protein whose domains occur N->C as DUF4011 -> SF1 helicase (DNA2/Upf1-like,
Pfam AAA_11/AAA_12, often split by a long insertion) -> REase_MTES_1575 -> DUF3320, with
DUF3320 absent from ~54% of Prometheus (van den Berg et al. 2024, Cell Host Microbe
32:1427). Every line of evidence below is computed per locus and written to the catalog,
so each tier call can be audited.

* ARCHITECTURE -- only through the shared domain spec (utils.load_domain_spec /
  architecture_tokens / is_canonical_order / required_categories / other_domains).
  "Core complete" = every `required: true` category present; DUF3320 is recorded
  (has_DUF3320) but never required. other_domains lists Pfam hits outside the spec
  (fusions such as the N-terminal ResIII/Helicase_C extension of haloarchaeal ProA).
* CATALYTIC STATUS -- residue-level evidence when the hmmscan alignment output is given
  (--hmmscan-out): residues aligned to HMM match states are read from hmmscan's own
  alignments. Helicase module (helicase category, SF1): Walker A / motif I
  [AG]x4GK[ST] and Walker B / motif II DE, located on the HMM consensus (anchor residues
  upper-case = match emission >= 0.5, HMMER p7_hmm_SetConsensus) and required in the
  candidate. Nuclease module (nuclease category, REase_MTES_1575 = Vsr-REase fold,
  Pfam/InterPro PF18741): every upper-case (>= 0.5) D/E match state of the consensus
  must carry D or E in the candidate. PF18741 has no canonical (D/E)xK, so a
  PD-(D/E)XK regex would be wrong for this domain; these acidic residues are the
  putative metal-binding core, not experimentally mapped for ProA.
  intact = both modules intact; partial = one module compromised (absent / motif
  disrupted); degenerate = both compromised; unknown = no compromised module but at
  least one not assessable. Without --hmmscan-out the old domain-presence proxy is used
  and labelled `catalytic_evidence=domain_proxy`.
* PHYLOGENETIC PLACEMENT -- the validation tree is rooted on the paralog (DNA2/NurA/
  HerA) reference outgroup: the root edge is the edge whose one side holds every
  Paralog_ref and no Prometheus_ref, taken closest to the paralogs. A candidate is
  'Prometheus' when it lies on the ingroup side AND some clade containing it plus >= 1
  Prometheus_ref (and, by construction, no Paralog_ref) has SH-aLRT >= --min-alrt and
  UFBoot >= --min-ufboot (IQ-TREE guidance: 80 / 95); 'Paralog' symmetrically on the
  outgroup side. Unsupported, untestable (no tree, < 1 reference of either kind,
  references not separable) or absent tips are 'Unplaced'; the reason is in
  placement_note. placement_clade says crown (inside the MRCA of the Prometheus
  references) or stem (ingroup, outside that crown).
* DEFENSE CONTEXT -- a true +/-N-ORF neighbourhood test (--context-window, design 1.5):
  genes of OTHER DefenseFinder systems (not the candidate's own system and not other
  systems of the focal type) on the same contig within N genes. Gene order comes from
  genes/{g}/{g}.genes.tsv (--genes-dir), else the contig-aware gene id
  {genome}_{contig#}_{gene#}, else DefenseFinder's own replicon + hit_pos. Contigs are
  treated as linear. Context only ever raises confidence (tier B), never rejects.
* FRAGMENTATION -- Prodigal's partial code from genes.tsv ('00' complete; '10'/'01'/'11'
  run off a contig end) when --genes-dir is given; otherwise a labelled length proxy
  (length < --min-len). Fragmented loci are kept ("present-fragmented", design 1.7).

Tiers: Rejected (supported paralog placement, or a COMPLETE ORF < 0.5*min_len or >
2*max_len); A (DefenseFinder call + core complete in canonical order + supported
Prometheus placement + catalytic intact + not fragmented + length within
[min_len, max_len]); B (supported Prometheus placement, or core complete + defense
context); C (everything else retained). `reason` lists the basis and every shortfall.
"""
import os
import re
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import math
import pandas as pd

BASE_COLUMNS = [
    "locus_id", "genome_id", "protein_id", "length_aa", "detect_defensefinder",
    "architecture", "n_domains", "catalytic_status", "defense_context",
    "tree_assignment", "fragmented", "confidence_tier", "reason",
]
# Additive columns (SCHEMA Contract v2 item 5: never rename/remove, only add).
EXTRA_COLUMNS_HEAD = ["missing_required", "canonical_order", "other_domains"]
EXTRA_COLUMNS_TAIL = [
    "helicase_module", "nuclease_module", "catalytic_evidence",
    "placement_support", "placement_clade", "placement_note",
    "context_method", "context_n_systems", "context_systems",
    "partial", "fragment_evidence", "length_in_bounds",
]
# Backwards-compatible constant (13 base columns) for callers that import it.
CATALOG_HEADER = "\t".join(BASE_COLUMNS)
SENSITIVITY_HEADER = "genome_id\tdetect_defensefinder\tprometheus_present"


def catalog_columns(spec) -> list:
    """Full catalog column list: base + has_<category> per spec category + extras."""
    return (BASE_COLUMNS + EXTRA_COLUMNS_HEAD
            + [f"has_{c}" for c in utils.canonical_order(spec)] + EXTRA_COLUMNS_TAIL)


class InputError(RuntimeError):
    """Inconsistent or malformed input: the step must fail, not degrade silently."""


# --------------------------------------------------------------------------- #
# Small IO helpers
# --------------------------------------------------------------------------- #
def strip_stop(seq: str) -> str:
    """Prodigal writes a terminal '*' for complete genes; it is not a residue."""
    return seq[:-1] if seq.endswith("*") else seq


def read_id_list(path: str) -> list:
    """Genome ids from a TSV with a genome_id column, or a headerless one-per-line list
    (representatives.txt). Returns [] when the path is empty/absent."""
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return []
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
    if "\t" in first or first.strip() == "genome_id":
        df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        col = "genome_id" if "genome_id" in df.columns else df.columns[0]
        return [g.strip() for g in df[col].tolist() if g.strip()]
    with open(path) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


# --------------------------------------------------------------------------- #
# Domain architecture (exclusively through the shared spec helpers)
# --------------------------------------------------------------------------- #
def architecture_for(hits: list, spec: dict):
    """(architecture_string, n_distinct_categories, present_set) for one protein.

    Kept for callers (architecture_typology fallback). Tokens come from
    utils.architecture_tokens: N->C by alignment start, consecutive pieces of one
    category collapsed, genuine re-occurrences kept (rearrangements stay visible).
    """
    toks = utils.architecture_tokens(hits, spec)
    present = set(toks)
    return "-".join(toks), len(present), present


def architecture_record(hits: list, spec: dict) -> dict:
    toks = utils.architecture_tokens(hits, spec)
    present = set(toks)
    required = utils.required_categories(spec)
    missing = [c for c in required if c not in present]
    return {
        "tokens": toks,
        "architecture": "-".join(toks),
        "n_domains": len(present),
        "present": present,
        "missing_required": missing,
        # canonical = every present category once, in spec order (and core complete)
        "canonical_order": bool(toks) and utils.is_canonical_order(toks, spec),
        "other_domains": utils.other_domains(hits, spec),
    }


# --------------------------------------------------------------------------- #
# hmmscan main-output alignments -> residues per HMM match state
# --------------------------------------------------------------------------- #
_ALN_ROW = re.compile(r"^\s*(\S+)\s+(\d+)\s+(\S+)\s+(\d+)\s*$")


def parse_hmmscan_alignments(path: str) -> dict:
    """Parse `hmmscan -o` (best with --notextw) into
    {query: [{"model", "hmm_from", "hmm_to", "ali_from", "ali_to", "cons", "tgt", "pos"}]}

    For every reported domain, `cons` is the model consensus over MATCH states only
    (hmm_from..hmm_to, case preserved: upper = emission >= 0.5), `tgt` the candidate
    residue aligned to each match state ('-' for a deletion) and `pos` its 1-based
    sequence coordinate (None for deletions). Insert columns ('.' in the model row) are
    skipped. Wrapped (non --notextw) alignments are concatenated block by block.
    """
    out = {}
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return out
    query, model = None, None
    blocks = []          # rows of the current domain alignment: (name, start, seq, end)
    in_aln = False

    def flush():
        if query is None or model is None or not blocks:
            return
        mrows = [b for b in blocks if b[0] == model]
        trows = [b for b in blocks if b[0] != model]
        if not mrows or not trows or len(mrows) != len(trows):
            return
        mseq = "".join(b[2] for b in mrows)
        tseq = "".join(b[2] for b in trows)
        if len(mseq) != len(tseq):
            return
        hmm_from, hmm_to = mrows[0][1], mrows[-1][3]
        ali_from, ali_to = trows[0][1], trows[-1][3]
        cons, tgt, pos = [], [], []
        spos = ali_from - 1
        for m, t in zip(mseq, tseq):
            is_res = t not in "-."
            if is_res:
                spos += 1
            if m == ".":            # insert column (residue not aligned to a state)
                continue
            cons.append(m)
            tgt.append(t.upper() if is_res else "-")
            pos.append(spos if is_res else None)
        out.setdefault(query, []).append({
            "model": model, "hmm_from": hmm_from, "hmm_to": hmm_to,
            "ali_from": ali_from, "ali_to": ali_to,
            "cons": "".join(cons), "tgt": "".join(tgt), "pos": pos})

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("Query:"):
                flush(); blocks = []; in_aln = False; model = None
                query = line.split()[1]
                continue
            if line.startswith(">>"):
                flush(); blocks = []; in_aln = False
                model = line.split()[1]
                continue
            s = line.strip()
            if s.startswith("== domain"):
                flush(); blocks = []; in_aln = True
                continue
            if s == "//" or s.startswith("Internal pipeline statistics"):
                flush(); blocks = []; in_aln = False
                continue
            if not in_aln or not s:
                continue
            if s.endswith(" PP") or s.endswith(" CS") or s.endswith(" RF"):
                continue
            m = _ALN_ROW.match(line)
            if m and (m.group(1) == model or m.group(1) == query):
                blocks.append((m.group(1), int(m.group(2)), m.group(3), int(m.group(4))))
    flush()
    return out


WALKER_A_CONS = re.compile(r"[AaGg][A-Za-z]{4}GK[ST]")   # anchors G,K,[ST] high-information
MOTIF_II_CONS = re.compile(r"DE")                         # both residues high-information
# Sequence-level confirmation used ONLY where no HMM alignment covers a motif (Pfam's
# AAA_11 C-terminal piece is often below GA in ProA). Strict patterns, searched in an
# architecture-anchored window, can confirm presence but never call a motif disrupted:
#   motif I  (P-loop)          G-x2-G-x-G-K-[ST]           (ProA: GPPGTGKS)
#   motif II + motif III (SF1) h3-D-E-x(18-32)-G-D-x2-Q    (ProA: LVIFDEASQ...VVGDPKQ)
WALKER_A_SEQ = re.compile(r"G..G.GK[ST]")
MOTIF_II_III_SEQ = re.compile(r"[ACFGILMVWY]{3}(DE).{18,32}GD..Q")


def helicase_windows(seq_len: int, hits: list, spec: dict, helicase_cat: str):
    """(walker_window, motif2_window) as 1-based inclusive (lo, hi), or None.
    The helicase region runs from the end of the last hit of a category that precedes the
    helicase in canonical order to the start of the first hit of a category that follows
    it (DUF4011 ... REase in ProA). Motif I is searched across that region up to the end
    of the last helicase hit; motif II from 250 residues before the last helicase hit
    (the C-terminal RecA-like piece) to the region end."""
    order = utils.canonical_order(spec)
    if helicase_cat not in order:
        return None
    cat_of = [(h, utils.hit_category(h, spec)) for h in hits]
    hel = sorted([h for h, c in cat_of if c == helicase_cat], key=lambda h: h["start"])
    if not hel:
        return None
    k = order.index(helicase_cat)
    before, after = set(order[:k]), set(order[k + 1:])
    first, last = hel[0], hel[-1]
    lo = max([h["end"] for h, c in cat_of if c in before and h["end"] <= first["start"]] or [0]) + 1
    hi = min([h["start"] for h, c in cat_of if c in after and h["start"] >= last["end"]]
             or [seq_len + 1]) - 1
    return (lo, max(lo, last["end"])), (max(lo, last["start"] - 250), max(lo, hi))


def _assess_helicase(alns: list, seq: str = "", windows=None):
    """Walker A (motif I) + Walker B (motif II) of the SF1 helicase for one protein.
    HMM-anchored first; strict window patterns only confirm presence. Returns
    (status, evidence) with status in intact / disrupted / unassessed."""
    wa, m2 = [], []          # per assessed copy: (ok, text)
    for a in alns:
        cons, tgt, pos = a["cons"], a["tgt"], a["pos"]
        for mt in WALKER_A_CONS.finditer(cons):
            i = mt.start()
            seg = tgt[i:i + 8]
            if len(seg) < 8:
                continue
            ok = seg[0] in "AG" and seg[5] == "G" and seg[6] == "K" and seg[7] in "ST"
            p = pos[i + 6]
            wa.append((ok, f"WalkerA:{seg}@{p if p else 'del'}({a['model']})"))
        for mt in MOTIF_II_CONS.finditer(cons):
            i = mt.start()
            seg = tgt[i:i + 2]
            ok = seg == "DE"
            p = pos[i]
            m2.append((ok, f"MotifII:{seg}@{p if p else 'del'}({a['model']})"))
    if windows and seq:
        (wlo, whi), (mlo, mhi) = windows
        if not any(ok for ok, _ in wa):
            m = WALKER_A_SEQ.search(seq, max(0, wlo - 1), whi)
            if m:
                wa.append((True, f"WalkerA:{m.group()}@{m.start() + 7}(window)"))
        if not any(ok for ok, _ in m2):
            m = MOTIF_II_III_SEQ.search(seq, max(0, mlo - 1), mhi)
            if m:
                m2.append((True, f"MotifII:DE@{m.start(1) + 1}(window)"))
    parts, status = [], "intact"
    for name, found in (("WalkerA", wa), ("MotifII", m2)):
        good = [t for ok, t in found if ok]
        if good:
            parts.append(good[0])
            continue
        if found:            # aligned to the HMM and not conserved
            parts.append(found[0][1])
            status = "disrupted"
        else:
            parts.append(f"{name}:not_assessable")
            if status == "intact":
                status = "unassessed"
    return status, ",".join(parts)


def acidic_sites_by_model(alignments: dict) -> dict:
    """{model: sorted HMM match states whose consensus residue is an upper-case D/E}.
    The consensus is a property of the model, so the union over every candidate's
    alignment recovers the sites even where one candidate's alignment is truncated."""
    sites = {}
    for alns in alignments.values():
        for a in alns:
            st = sites.setdefault(a["model"], set())
            for i, c in enumerate(a["cons"]):
                if c in "DE":
                    st.add(a["hmm_from"] + i)
    return {m: sorted(v) for m, v in sites.items()}


def _assess_nuclease(alns: list, model_sites: dict = None):
    """Every high-information acidic match state of the nuclease-category model must
    keep D/E. Sites outside a candidate's alignment are 'not aligned'. The copy with the
    most assessable sites decides. Returns (status, evidence) with status in
    intact / disrupted / unassessed."""
    best = None
    for a in alns:
        cons, tgt, pos = a["cons"], a["tgt"], a["pos"]
        sites = (model_sites or {}).get(a["model"]) or [
            a["hmm_from"] + i for i, c in enumerate(cons) if c in "DE"]
        if not sites:
            continue
        res, n_missing = [], 0
        for state in sites:
            i = state - a["hmm_from"]
            if i < 0 or i >= len(cons):
                n_missing += 1
                res.append((None, f"{state}:not_aligned"))
                continue
            ok = tgt[i] in "DE"
            res.append((ok, f"{cons[i]}{state}:{tgt[i]}@{pos[i] if pos[i] else 'del'}"))
        n_ok = sum(1 for ok, _ in res if ok)
        n_bad = sum(1 for ok, _ in res if ok is False)
        cand = (len(res) - n_missing, n_ok, n_bad, res, a["model"])
        if best is None or (cand[0], cand[1]) > (best[0], best[1]):
            best = cand
    if best is None:
        return "unassessed", "acidic_core:not_aligned"
    n_assessed, n_ok, n_bad, res, model = best
    if n_bad:
        status = "disrupted"
    elif n_assessed < len(res):
        status = "unassessed"
    else:
        status = "intact"
    return status, f"acidic_core({model}):" + ",".join(t for _, t in res)


def catalytic_assessment(present: set, alns: list, spec: dict, helicase_cat: str,
                         nuclease_cat: str, have_alignments: bool, has_hits: bool,
                         model_sites: dict = None, seq: str = "", hits: list = None):
    """(catalytic_status, helicase_module, nuclease_module, evidence) -- see module doc."""
    if not has_hits:
        return "unknown", "absent", "absent", "no_domain_hits"
    if not have_alignments:
        # Labelled fallback: domain presence only (no residue evidence available).
        h = helicase_cat in present
        n = nuclease_cat in present
        status = "intact" if (h and n) else ("partial" if (h or n) else "degenerate")
        return (status, "present" if h else "absent", "present" if n else "absent",
                "domain_proxy")
    by_cat = {}
    for a in alns:
        cat = utils.hit_category({"domain": a["model"], "accession": a.get("accession", "")}, spec)
        if cat is not None:
            by_cat.setdefault(cat, []).append(a)
    mods = {}
    ev = []
    for label, cat, fn in (("helicase", helicase_cat, _assess_helicase),
                           ("nuclease", nuclease_cat, _assess_nuclease)):
        if cat not in present:
            mods[label] = "absent"
            ev.append(f"{label}=absent")
            continue
        if not by_cat.get(cat) and not (label == "helicase" and seq):
            mods[label] = "unassessed"
            ev.append(f"{label}=unassessed(no_alignment)")
            continue
        if label == "helicase":
            st, txt = fn(by_cat.get(cat, []), seq,
                         helicase_windows(len(seq), hits or [], spec, cat) if seq else None)
        else:
            st, txt = fn(by_cat[cat], model_sites)
        mods[label] = st
        ev.append(f"{label}={st}[{txt}]")
    bad = [m for m in mods.values() if m in ("absent", "disrupted")]
    if not bad:
        status = "intact" if all(m == "intact" for m in mods.values()) else "unknown"
    elif len(bad) == 2:
        status = "degenerate"
    else:
        status = "partial"
    return status, mods["helicase"], mods["nuclease"], ";".join(ev)


# --------------------------------------------------------------------------- #
# Phylogenetic placement (paralog-outgroup rooting + supported clade membership)
# --------------------------------------------------------------------------- #
def load_ref_labels(path: str) -> dict:
    """ref tip name -> 'Prometheus_ref' | 'Paralog_ref' (headerless TSV; names are
    canonicalised with utils.sanitize so they match IQ-TREE-renamed tips)."""
    labels = {}
    if not path or not os.path.exists(path):
        return labels
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            if len(parts) < 2:
                continue
            name, grp = parts[0].strip(), parts[1].strip()
            if grp in ("Prometheus_ref", "Paralog_ref"):
                key = utils.sanitize(name)
                if key in labels and labels[key] != grp:
                    raise InputError(f"reference '{name}' is labelled both Prometheus_ref "
                                     f"and Paralog_ref in {path}")
                labels[key] = grp
    return labels


def parse_support(label):
    """IQ-TREE internal-node label -> (SH-aLRT, UFBoot); NaN where absent.
    'a/b' = SH-aLRT/UFBoot (iqtree2 --alrt + -B); 'a/b/c' = SH-aLRT/aBayes/UFBoot;
    a single number is treated as a bootstrap (UFBoot-scale) value."""
    nan = float("nan")
    if label is None:
        return nan, nan
    s = str(label).strip()
    if not s:
        return nan, nan
    try:
        vals = [float(x) for x in s.split("/")]
    except ValueError:
        return nan, nan
    if len(vals) == 1:
        return nan, vals[0]
    if len(vals) == 2:
        return vals[0], vals[1]
    if len(vals) == 3:
        return vals[0], vals[2]
    return nan, nan


def is_supported(alrt: float, ufb: float, min_alrt: float, min_ufb: float) -> bool:
    have = False
    if not math.isnan(ufb):
        have = True
        if ufb < min_ufb:
            return False
    if not math.isnan(alrt):
        have = True
        if alrt < min_alrt:
            return False
    return have


def _fmt_support(alrt, ufb) -> str:
    def f(x):
        return "NA" if math.isnan(x) else (f"{x:g}")
    return f"{f(alrt)}/{f(ufb)}"


def assign_by_tree(tree_path: str, ref_labels: dict, min_alrt: float = 80.0,
                   min_ufb: float = 95.0, details: dict = None):
    """{candidate tip -> 'Prometheus'|'Paralog'|'Unplaced'} (see module docstring).

    When `details` (a dict) is given it is filled with
    {tip: {"support", "clade", "note"}} for every candidate tip in the tree, plus the
    key "__tree__" with a tree-level note. Tips are canonicalised with utils.sanitize.
    """
    details = details if details is not None else {}
    assignments = {}

    def tree_note(msg):
        details["__tree__"] = msg
        utils.log(f"placement: {msg}")

    if not tree_path or not os.path.exists(tree_path) or os.path.getsize(tree_path) == 0:
        tree_note("no_tree")
        return assignments
    txt = open(tree_path).read().strip()
    if not txt or re.fullmatch(r"\(\s*\)\s*;?", txt):
        tree_note("no_tree (placeholder: too few sequences for a tree)")
        return assignments
    import dendropy  # hard dependency: a missing package is an error, not 'Unplaced'
    try:
        tree = dendropy.Tree.get(data=txt, schema="newick", preserve_underscores=True,
                                 rooting="force-unrooted")
    except Exception as exc:  # malformed tree from a crashed tool must not pass silently
        raise InputError(f"could not parse validation tree {tree_path}: {exc}")

    leaves = [nd for nd in tree.leaf_node_iter() if nd.taxon is not None]
    labels = [utils.sanitize(nd.taxon.label) for nd in leaves]
    if len(set(labels)) != len(labels):
        raise InputError(f"duplicate tip labels in {tree_path}")
    index = {lab: i for i, lab in enumerate(labels)}
    bit = {nd: 1 << index[utils.sanitize(nd.taxon.label)] for nd in leaves}
    ALL = (1 << len(labels)) - 1
    prom = sum(1 << index[l] for l in labels if ref_labels.get(l) == "Prometheus_ref")
    para = sum(1 << index[l] for l in labels if ref_labels.get(l) == "Paralog_ref")
    cands = [l for l in labels if l not in ref_labels]
    missing_refs = [l for l in ref_labels if l not in index]
    if missing_refs:
        utils.log(f"WARNING: {len(missing_refs)} reference label(s) not found among tree "
                  f"tips (e.g. {missing_refs[:3]})")
    if not prom or not para:
        tree_note("no_reference_contrast (need >= 1 Prometheus_ref and >= 1 Paralog_ref in tree)")
        for c in cands:
            assignments[c] = "Unplaced"
            details[c] = {"support": "", "clade": "", "note": "no_reference_contrast"}
        return assignments

    # every edge as (one-side bitmask, alrt, ufb); terminal edges carry no support
    masks = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            masks[nd] = bit.get(nd, 0)
        else:
            m = 0
            for ch in nd.child_node_iter():
                m |= masks[ch]
            masks[nd] = m
    splits = []
    seen = {}
    for nd in tree.preorder_node_iter():
        if nd.parent_node is None:
            continue
        side = masks[nd]
        if side == 0 or side == ALL:
            continue
        alrt, ufb = (float("nan"), float("nan")) if nd.is_leaf() else parse_support(nd.label)
        key = min(side, ALL ^ side)
        if key in seen:      # a bifurcating root lists one edge twice: keep max support
            j = seen[key]
            a0, u0 = splits[j][1], splits[j][2]
            splits[j] = (splits[j][0], max(a0, alrt) if not math.isnan(alrt) else a0,
                         max(u0, ufb) if not math.isnan(ufb) else u0)
            continue
        seen[key] = len(splits)
        splits.append((side, alrt, ufb))

    # root edge: one side holds ALL paralog refs and NO Prometheus ref; closest to paralogs
    root = None
    for side, alrt, ufb in splits:
        for a in (side, ALL ^ side):
            if (a & para) == para and (a & prom) == 0:
                n = bin(a).count("1")
                if root is None or n < root[0]:
                    root = (n, a, alrt, ufb)
    if root is None:
        tree_note("references_not_separable (no edge separates every Paralog_ref from "
                  "every Prometheus_ref; check the reference sets)")
        for c in cands:
            assignments[c] = "Unplaced"
            details[c] = {"support": "", "clade": "", "note": "references_not_separable"}
        return assignments
    _, A, r_alrt, r_ufb = root
    B = ALL ^ A
    # clades of the paralog-rooted tree: for each edge, the side that is a subset of A or B
    clades = []
    for side, alrt, ufb in splits:
        for x in (side, ALL ^ side):
            if (x & ~A & ALL) == 0 or (x & ~B & ALL) == 0:
                clades.append((x, bin(x).count("1"), alrt, ufb,
                               is_supported(alrt, ufb, min_alrt, min_ufb)))
    # crown = smallest ingroup clade containing every Prometheus ref
    crown, crown_n = None, None
    for x, n, _, _, _ in clades:
        if (x & prom) == prom and (x & ~B & ALL) == 0:
            if crown is None or n < crown_n:
                crown, crown_n = x, n
    if crown is None:
        crown = prom
    details["__tree__"] = (f"rooted on paralog outgroup ({bin(A).count('1')} tips); "
                           f"root-edge support {_fmt_support(r_alrt, r_ufb)}")

    for c in cands:
        cb = 1 << index[c]
        ingroup = bool(cb & B)
        own, rootside = (prom, B) if ingroup else (para, A)
        best_sup, nearest = None, None
        root_clade = (rootside, bin(rootside).count("1"), r_alrt, r_ufb,
                      is_supported(r_alrt, r_ufb, min_alrt, min_ufb))
        for x, size, alrt, ufb, sup in clades + [root_clade]:
            if not (x & cb) or (x & ~rootside & ALL) or not (x & own):
                continue
            if nearest is None or size < nearest[0]:
                nearest = (size, alrt, ufb)
            if sup and (best_sup is None or size < best_sup[0]):
                best_sup = (size, alrt, ufb)
        side_name = "Prometheus" if ingroup else "Paralog"
        if ingroup:
            clade = "crown" if (cb & crown) else "stem"
        else:
            clade = "outgroup"
        if best_sup is not None:
            assignments[c] = side_name
            details[c] = {"support": _fmt_support(best_sup[1], best_sup[2]),
                          "clade": clade, "note": "supported"}
        else:
            assignments[c] = "Unplaced"
            sup = _fmt_support(nearest[1], nearest[2]) if nearest else ""
            details[c] = {"support": sup, "clade": clade,
                          "note": f"unsupported_{'ingroup' if ingroup else 'outgroup'}"
                                  f"(<{min_alrt:g}/{min_ufb:g})"}
    return assignments


# --------------------------------------------------------------------------- #
# Gene order (genes.tsv) and the +/-N-ORF defense-context test
# --------------------------------------------------------------------------- #
def load_genes_table(genes_dir: str, genome: str):
    """{local protein id -> dict(contig_index, gene_index, partial, length_aa)} from
    {genes_dir}/{g}/{g}.genes.tsv, or None when the file does not exist."""
    if not genes_dir:
        return None
    path = os.path.join(genes_dir, genome, f"{genome}.genes.tsv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    need = {"protein_id", "contig_index", "gene_index", "partial"}
    if not need.issubset(df.columns):
        raise InputError(f"{path} lacks columns {sorted(need - set(df.columns))} "
                         "(SCHEMA Contract v2 item 2)")
    out = {}
    for r in df.itertuples(index=False):
        try:
            ci, gi = int(r.contig_index), int(r.gene_index)
        except ValueError:
            raise InputError(f"{path}: non-integer contig/gene index for {r.protein_id}")
        out[str(r.protein_id)] = {"contig_index": ci, "gene_index": gi,
                                  "partial": str(r.partial).strip()}
    return out


def load_df_genes(df_dir: str, genome: str):
    """DefenseFinder genes table for one genome (DataFrame, str dtype) or None."""
    path = os.path.join(df_dir, genome, f"{genome}_defense_finder_genes.tsv")
    if not df_dir or not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    return df


def gene_coordinate(genome: str, local_id: str, genes_tbl):
    """(contig_index, gene_index, method) for a per-genome gene id, or None."""
    if genes_tbl is not None and local_id in genes_tbl:
        g = genes_tbl[local_id]
        return g["contig_index"], g["gene_index"], "genes_tsv"
    sp = utils.split_gene_id(local_id)
    if sp is not None and sp[0] == utils.sanitize(genome):
        return sp[1], sp[2], "gene_id"
    return None


def defense_context(genome: str, local_id: str, own_sys_ids: set, focal: str,
                    window: int, df_genes, genes_tbl):
    """(context_bool, method, n_systems, 'type:sys_id,...') for one candidate."""
    if df_genes is None or "hit_id" not in df_genes.columns or "sys_id" not in df_genes.columns:
        return False, "unavailable(no_defensefinder_genes)", 0, ""
    focal_l = focal.lower()
    other = df_genes[~df_genes["sys_id"].astype(str).isin(own_sys_ids)].copy()
    if "type" in other.columns:
        t = other["type"].astype(str).str.lower()
        st = other["subtype"].astype(str).str.lower() if "subtype" in other.columns else t
        other = other[(t != focal_l) & (st != focal_l)]
    me = gene_coordinate(genome, local_id, genes_tbl)
    hits = {}
    if me is not None:
        method = f"neighbourhood_+/-{window}_{me[2]}"
        for r in other.itertuples(index=False):
            co = gene_coordinate(genome, utils.sanitize(str(r.hit_id)), genes_tbl)
            if co is None:
                continue
            if co[0] == me[0] and abs(co[1] - me[1]) <= window:
                hits[str(r.sys_id)] = str(getattr(r, "type", "") or "")
    else:
        # DefenseFinder's own coordinates: same replicon, |hit_pos difference| <= N
        need = {"replicon", "hit_pos"}
        if not need.issubset(df_genes.columns):
            return False, "unavailable(no_gene_order)", 0, ""
        mine = df_genes[df_genes["hit_id"].map(lambda h: utils.sanitize(str(h))) == local_id]
        if mine.empty:
            return False, "unavailable(candidate_not_in_defensefinder_genes)", 0, ""
        rep, pos = str(mine.iloc[0]["replicon"]), int(float(mine.iloc[0]["hit_pos"]))
        method = f"neighbourhood_+/-{window}_replicon_hit_pos"
        for r in other.itertuples(index=False):
            try:
                p = int(float(r.hit_pos))
            except ValueError:
                continue
            if str(r.replicon) == rep and abs(p - pos) <= window:
                hits[str(r.sys_id)] = str(getattr(r, "type", "") or "")
    desc = ",".join(f"{t}:{s}" for s, t in sorted(hits.items()))
    return bool(hits), method, len(hits), desc


# --------------------------------------------------------------------------- #
# Tiering
# --------------------------------------------------------------------------- #
def assign_tier(detect_df, missing_required, canonical, tree_assignment, fragmented,
                catalytic, defense_ctx, length_aa, length_in_bounds, min_len, max_len,
                placement_note=""):
    """(confidence_tier, reason) -- see module docstring for the definitions."""
    if tree_assignment == "Paralog":
        return "Rejected", "placed_with_cellular_paralog"
    if length_aa and not fragmented:
        if length_aa < 0.5 * min_len:
            return "Rejected", f"complete_orf_far_below_min({length_aa}<0.5*{min_len})"
        if length_aa > 2.0 * max_len:
            return "Rejected", f"length_far_above_max({length_aa}>2*{max_len})"
    short = []
    if not detect_df:
        short.append("no_defensefinder_call")
    if missing_required:
        short.append("missing_required:" + ",".join(missing_required))
    elif not canonical:
        short.append("noncanonical_order")
    if tree_assignment != "Prometheus":
        short.append(f"placement:{placement_note or 'unplaced'}")
    if catalytic != "intact":
        short.append(f"catalytic:{catalytic}")
    if fragmented:
        short.append("fragmented")
    if not length_in_bounds:
        short.append("length_out_of_bounds")
    if not short:
        return "A", ""
    basis = []
    if tree_assignment == "Prometheus":
        basis.append("prometheus_placement")
    if not missing_required and defense_ctx:
        basis.append("core_architecture+defense_context")
    if basis:
        return "B", ";".join(basis + short)
    return "C", ";".join(short)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", required=True, help="candidate_sources.tsv")
    ap.add_argument("--raw", required=True, help="raw_hits.tsv (accepted for CLI compatibility)")
    ap.add_argument("--domtbl", required=True, help="hmmscan --domtblout vs Pfam")
    ap.add_argument("--tree", required=True, help="validation tree (newick, IQ-TREE supports)")
    ap.add_argument("--ref-labels", required=True, help="reference tip -> group TSV (headerless)")
    ap.add_argument("--candidates", required=True, help="candidate proteins FASTA")
    ap.add_argument("--meta", required=True, help="genome table with genome_id (the universe)")
    ap.add_argument("--reps", required=True,
                    help="optional restriction of the genome universe: TSV with genome_id "
                         "or a one-id-per-line list (representatives.txt); same file as "
                         "--meta is a no-op")
    ap.add_argument("--df-dir", required=True, help="DefenseFinder output root ({g}/{g}_*.tsv)")
    ap.add_argument("--context-window", type=int, default=10,
                    help="+/- N genes scanned for other defense systems")
    ap.add_argument("--domains", required=True, help="JSON domain spec (config prometheus.domains)")
    ap.add_argument("--min-len", type=int, required=True)
    ap.add_argument("--max-len", type=int, required=True)
    ap.add_argument("--out-catalog", required=True)
    ap.add_argument("--out-proteins", required=True)
    ap.add_argument("--out-sensitivity", required=True)
    # --- optional (Contract v2 additions; absent flags keep the old behaviour) ---
    ap.add_argument("--genes-dir", default="",
                    help="results/genes: genes/{g}/{g}.genes.tsv gives gene order + Prodigal "
                         "partial flags (else id-based order and a length proxy)")
    ap.add_argument("--hmmscan-out", default="",
                    help="hmmscan main output (-o, --notextw) for residue-level motifs")
    ap.add_argument("--min-alrt", type=float, default=80.0, help="SH-aLRT threshold (%%)")
    ap.add_argument("--min-ufboot", type=float, default=95.0, help="UFBoot threshold (%%)")
    ap.add_argument("--system", default="Prometheus",
                    help="focal DefenseFinder system type (excluded from context)")
    ap.add_argument("--helicase-category", default="DNA2",
                    help="spec category holding the SF1 helicase (Walker A / motif II)")
    ap.add_argument("--nuclease-category", default="REase",
                    help="spec category holding the nuclease (conserved acidic core)")
    args = ap.parse_args()

    for p in (args.out_catalog, args.out_proteins, args.out_sensitivity):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    try:
        spec = utils.load_domain_spec(args.domains)
    except ValueError as exc:   # malformed spec is a configuration error
        raise SystemExit(f"[validate_prometheus] FATAL: invalid --domains JSON: {exc}")
    if not spec:
        raise SystemExit("[validate_prometheus] FATAL: empty --domains spec")
    for cat in (args.helicase_category, args.nuclease_category):
        if cat not in spec:
            utils.log(f"WARNING: category '{cat}' not in the domain spec; its catalytic "
                      "module will always be 'absent'")
    columns = catalog_columns(spec)

    # ---- genome universe ------------------------------------------------- #
    universe = read_id_list(args.meta)
    reps = read_id_list(args.reps) if os.path.abspath(args.reps) != os.path.abspath(args.meta) else []
    if reps:
        keep = set(reps)
        universe = [g for g in universe if g in keep]

    # ---- candidate sources ----------------------------------------------- #
    sources = pd.read_csv(args.sources, sep="\t", dtype=str, keep_default_na=False)
    if "protein_id" not in sources.columns:
        raise SystemExit(f"[validate_prometheus] FATAL: {args.sources} has no protein_id column")

    def write_sensitivity(sens):
        rows = []
        for g in list(dict.fromkeys(universe + sorted(sens))):
            d, present = sens.get(g, (0, False))
            rows.append({"genome_id": g, "detect_defensefinder": d,
                         "prometheus_present": bool(present)})
        pd.DataFrame(rows, columns=SENSITIVITY_HEADER.split("\t")).to_csv(
            args.out_sensitivity, sep="\t", index=False)

    if sources.empty:
        utils.write_empty(args.out_catalog, "\t".join(columns))
        open(args.out_proteins, "w").close()
        write_sensitivity({})
        utils.log("No candidate sources; wrote header-only catalog.")
        return

    seqs = {pid: strip_stop(s) for pid, s in utils.read_fasta(args.candidates)}
    missing_seq = [p for p in sources["protein_id"] if p not in seqs]
    if missing_seq:
        raise SystemExit(f"[validate_prometheus] FATAL: {len(missing_seq)} candidate(s) in "
                         f"{args.sources} have no sequence in {args.candidates} "
                         f"(e.g. {missing_seq[:3]})")

    domtbl_hits = utils.parse_domtbl(args.domtbl)
    if args.hmmscan_out and not os.path.exists(args.hmmscan_out):
        raise SystemExit(f"[validate_prometheus] FATAL: --hmmscan-out {args.hmmscan_out} missing")
    alignments = parse_hmmscan_alignments(args.hmmscan_out) if args.hmmscan_out else {}
    have_aln = bool(args.hmmscan_out)
    model_sites = acidic_sites_by_model(alignments)
    n_dom = sum(len(v) for v in domtbl_hits.values())
    n_aln = sum(len(v) for v in alignments.values())
    if have_aln and n_dom and n_aln == 0:
        raise SystemExit("[validate_prometheus] FATAL: the domtbl has domain hits but no "
                         f"alignment could be parsed from {args.hmmscan_out} "
                         "(run hmmscan with -o FILE --notextw)")
    if domtbl_hits == {} and len(seqs) > 0:
        utils.log("WARNING: the domtbl holds no domain hits for any candidate")

    ref_labels = load_ref_labels(args.ref_labels)
    placement = {}
    tree_assign = assign_by_tree(args.tree, ref_labels, args.min_alrt, args.min_ufboot,
                                 details=placement)
    utils.log("validation tree: " + str(placement.get("__tree__", "")))

    # DefenseFinder system ids per candidate (its own system is excluded from context)
    raw_sys = {}
    rows = []
    proteins_out = []
    sens = {}
    genes_cache, dfg_cache = {}, {}
    for r in sources.itertuples(index=False):
        pid = str(r.protein_id)
        gid = utils.genome_of(pid)
        local = utils.local_of(pid)
        d_df = 1 if str(getattr(r, "detect_defensefinder", "1")).strip() not in ("0", "", "False") else 0
        seq = seqs[pid]
        length_aa = len(seq)

        hits = domtbl_hits.get(pid, [])
        acc_by_name = {h["domain"]: h.get("accession", "") for h in hits}
        alns = alignments.get(pid, [])
        for a in alns:
            a["accession"] = acc_by_name.get(a["model"], "")
        arch = architecture_record(hits, spec)
        cat_status, hel, nuc, cat_ev = catalytic_assessment(
            arch["present"], alns, spec, args.helicase_category, args.nuclease_category,
            have_aln, has_hits=bool(hits), model_sites=model_sites, seq=seq, hits=hits)

        assignment = tree_assign.get(pid, "Unplaced")
        pdet = placement.get(pid)
        if pdet is None:
            tn = placement.get("__tree__", "no_tree")
            note = "not_in_tree" if tn.startswith("rooted") else tn.split(" ")[0]
            pdet = {"support": "", "clade": "", "note": note}

        if gid not in genes_cache:
            genes_cache[gid] = load_genes_table(args.genes_dir, gid)
        genes_tbl = genes_cache[gid]
        partial = ""
        if genes_tbl is not None:
            if local not in genes_tbl:
                raise SystemExit(f"[validate_prometheus] FATAL: {local} missing from "
                                 f"{args.genes_dir}/{gid}/{gid}.genes.tsv")
            partial = genes_tbl[local]["partial"]
            fragmented = partial not in ("00", "")
            frag_ev = f"prodigal_partial={partial}"
        else:
            fragmented = bool(d_df and 0 < length_aa < args.min_len)
            frag_ev = "length_proxy(length<min_len)" if fragmented else "length_proxy"
        in_bounds = args.min_len <= length_aa <= args.max_len

        if gid not in dfg_cache:
            dfg_cache[gid] = load_df_genes(args.df_dir, gid)
        dfg = dfg_cache[gid]
        own = set()
        if dfg is not None and "hit_id" in dfg.columns and "sys_id" in dfg.columns:
            own = set(dfg.loc[dfg["hit_id"].map(lambda h: utils.sanitize(str(h))) == local,
                              "sys_id"].astype(str))
        ctx, ctx_method, ctx_n, ctx_desc = defense_context(
            gid, local, own, args.system, args.context_window, dfg, genes_tbl)

        tier, reason = assign_tier(
            d_df, arch["missing_required"], arch["canonical_order"], assignment,
            fragmented, cat_status, ctx, length_aa, in_bounds, args.min_len, args.max_len,
            placement_note=pdet["note"])

        row = {
            "locus_id": pid, "genome_id": gid, "protein_id": pid, "length_aa": length_aa,
            "detect_defensefinder": d_df, "architecture": arch["architecture"],
            "n_domains": arch["n_domains"], "catalytic_status": cat_status,
            "defense_context": bool(ctx), "tree_assignment": assignment,
            "fragmented": bool(fragmented), "confidence_tier": tier, "reason": reason,
            "missing_required": ",".join(arch["missing_required"]),
            "canonical_order": bool(arch["canonical_order"]),
            "other_domains": ",".join(arch["other_domains"]),
            "helicase_module": hel, "nuclease_module": nuc, "catalytic_evidence": cat_ev,
            "placement_support": pdet["support"], "placement_clade": pdet["clade"],
            "placement_note": pdet["note"],
            "context_method": ctx_method, "context_n_systems": ctx_n,
            "context_systems": ctx_desc,
            "partial": partial, "fragment_evidence": frag_ev,
            "length_in_bounds": bool(in_bounds),
        }
        for c in utils.canonical_order(spec):
            row[f"has_{c}"] = c in arch["present"]
        rows.append(row)

        s = sens.setdefault(gid, [0, False])
        s[0] = max(s[0], d_df)
        if tier != "Rejected":
            s[1] = True
            proteins_out.append((pid, seq))

    pd.DataFrame(rows, columns=columns).to_csv(args.out_catalog, sep="\t", index=False)
    if proteins_out:
        utils.write_fasta(proteins_out, args.out_proteins)
    else:
        open(args.out_proteins, "w").close()
    outside = sorted(set(sens) - set(universe))
    if universe and outside:
        utils.log(f"WARNING: {len(outside)} candidate genome(s) are not in the --meta/--reps "
                  f"universe (e.g. {outside[:3]}); they are added to the sensitivity table")
    write_sensitivity({g: tuple(v) for g, v in sens.items()})

    tiers = pd.Series([r["confidence_tier"] for r in rows]).value_counts().to_dict()
    utils.log(f"Catalog: {len(rows)} loci, tiers {tiers}; {len(proteins_out)} retained "
              f"proteins written; {len(set(universe) | set(sens))} genomes in sensitivity.")


if __name__ == "__main__":
    try:
        main()
    except InputError as exc:
        raise SystemExit(f"[validate_prometheus] FATAL: {exc}")
