"""Classify each genome as methanogen / anaerobic methanotroph (ANME) / other from its
credible Mcr genes, separating true methanogens from ANME and from divergent McrA-like
enzymes by McrA phylogenetic placement (design §3.3).

Definitions
-----------
credible Mcr  every subunit in --required (config methanogen.mcr_required; KOs or gene
              names, default K00399 = McrA) has an ACCEPTED hit in mcr_hits.tsv (the
              curated KOfam / GA threshold, or a contig-edge fragment; mcr_search.py).
mcrA          an accepted McrA hit exists.
n_mcr         number of distinct core subunits (mcrA / mcrB / mcrG) with accepted hits.
placement     (--anme-tree True) the representative McrA of every credible genome is
              added to the curated, labelled reference alignment with
              `mafft --add --keeplength` (reference columns kept, query insertions
              dropped); an ML tree is inferred (IQ-TREE 2, LG+G, -fast, fixed --seed)
              and midpoint-rooted; each candidate is assigned by the NEAREST LABELLED
              CLADE rule: walk up from the candidate to the first node whose subtree
              holds reference sequences; if those references share one label the
              candidate takes it, otherwise it is 'ambiguous'. Candidates with fewer
              than --min-aligned residues in reference columns are 'unplaced'. The
              patristically nearest reference is reported as secondary evidence only
              (it is sensitive to long branches).
              Reference labels (ref_id<TAB>group): 'methanogen', 'other' (e.g. an
              alkyl-CoM reductase) or any label starting with 'ANME' (ANME-1, ANME-2a,
              ...), case-insensitive. Any other label is an error.
methanogen    credible Mcr AND (placement disabled OR mcr_group == 'methanogen').
is_anme       mcr_group == 'ANME' (ANME run Mcr in reverse; never counted as methanogens).
              'other' and 'ambiguous' / 'unplaced' genomes are not methanogens either;
              they are counted in the log and flagged in mcr_group / evidence so that a
              sensitivity analysis can revisit them.
methanogenesis_subtype  'unresolved' for methanogens, 'NA' otherwise. The subtype
              (hydrogenotrophic / methylotrophic / acetoclastic / methyl-reducing) needs
              pathway markers (mta/mtb/mtt methyltransferases, Mtr, CODH/ACS,
              hydrogenases) that this step does not test; mcr genes cannot resolve it,
              so it is never asserted.

Output: one row per genome of the universe (--pmap genomes plus any genome with hits),
SCHEMA columns first, then additive columns (mcr_group, mcrA_protein, mcrA_call,
n_mcrA_copies, nearest_ref, nearest_ref_group). --out-placement: one row per placed
candidate.

Failure policy
--------------
With --anme-tree True the reference alignment and labels must exist and be consistent
(every aligned reference labelled, ANME and methanogen labels both present, equal row
lengths, >= 3 sequences), and MAFFT and IQ-TREE must succeed; otherwise the script
fails. A candidate McrA that is missing from the proteome is an error. Only zero
candidates (no credible McrA anywhere) skip the placement, and that is recorded.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import csv
import shutil
import subprocess

from mcr_search import KO_TO_GENE  # one KO -> subunit table for the stage

HEADER_COLS = ["genome_id", "mcrA", "n_mcr", "methanogen", "is_anme",
               "methanogenesis_subtype", "evidence", "mcr_group", "mcrA_protein",
               "mcrA_call", "n_mcrA_copies", "nearest_ref", "nearest_ref_group"]
HEADER = "\t".join(HEADER_COLS)
PLACEMENT_COLS = ["genome_id", "protein_id", "aligned_residues", "clade_ref_groups",
                  "mcr_group", "nearest_ref", "nearest_ref_group", "patristic_distance"]
CORE = ("mcrA", "mcrB", "mcrG")
GENES = ("mcrA", "mcrB", "mcrG", "mcrC", "mcrD")
ACCEPTED_CALLS = {"full_threshold", "partial_rescue"}


# --------------------------------------------------------------------------- #
def norm_group(label: str) -> str:
    g = (label or "").strip()
    low = g.lower()
    if low.startswith("anme"):
        return "ANME"
    if low == "methanogen":
        return "methanogen"
    if low == "other":
        return "other"
    raise ValueError(f"reference label {label!r} is not one of methanogen | other | ANME*")


def load_ref_labels(path: str) -> dict:
    """ref id -> raw label; TSV id<TAB>group, optional header (ref_id/id, group)."""
    labels = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}: line without a label: {line!r}")
            rid, grp = parts[0].strip(), parts[1].strip()
            if rid.lower() in ("ref_id", "id") and grp.lower() == "group":
                continue  # header
            labels[rid] = grp
    return labels


def parse_required(spec: str) -> list:
    out = []
    for tok in [t.strip() for t in (spec or "").split(",") if t.strip()]:
        if tok in KO_TO_GENE:
            out.append(KO_TO_GENE[tok])
            continue
        match = [g for g in GENES if g.lower() == tok.lower()]
        if not match:
            raise ValueError(f"--required token {tok!r} is neither an mcr KO nor a subunit name")
        out.append(match[0])
    if not out:
        raise ValueError("--required is empty")
    return sorted(set(out))


def read_mcr_hits(path: str):
    """genome -> {gene: {'protein_id', 'call', 'n_copies'}} for ACCEPTED hits."""
    per = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None:
            return per
        need = {"genome_id", "gene", "protein_id"}
        if not need <= set(rd.fieldnames):
            raise ValueError(f"{path} lacks columns {sorted(need - set(rd.fieldnames))}")
        has_call = "call" in rd.fieldnames
        for r in rd:
            gene = r["gene"].strip()
            if gene not in GENES:
                continue
            call = r["call"].strip() if has_call else "full_threshold"
            if call not in ACCEPTED_CALLS:
                continue
            n = int(r["n_copies"]) if r.get("n_copies", "").strip().isdigit() else 1
            per.setdefault(r["genome_id"].strip(), {})[gene] = {
                "protein_id": r["protein_id"].strip(), "call": call, "n_copies": n}
    return per


def read_universe(pmap: str) -> set:
    out = set()
    if not pmap:
        return out
    with open(pmap) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        if rd.fieldnames is None:
            return out
        if "genome_id" not in rd.fieldnames:
            raise ValueError(f"{pmap} has no genome_id column")
        for r in rd:
            g = r["genome_id"].strip()
            if g:
                out.add(g)
    return out


def fetch_sequences(proteome: str, wanted: set) -> dict:
    """Stream the proteome and keep only `wanted` ids (the combined proteome is large)."""
    out = {}
    for name, seq in utils.read_fasta(proteome):
        if name in wanted:
            out[name] = seq
    return out


def clean_protein(seq: str) -> str:
    return seq.strip().rstrip("*").replace("*", "X").upper()


def read_alignment(path: str) -> list:
    recs = [(n, s.replace(".", "-").upper()) for n, s in utils.read_fasta(path)]
    return recs


def validate_reference(ref_aln: str, ref_labels: str):
    """-> (records, {id: normalised group}, {id: raw label}); raise on any inconsistency."""
    if not ref_aln or not os.path.exists(ref_aln) or os.path.getsize(ref_aln) == 0:
        raise FileNotFoundError(f"McrA reference alignment missing or empty: {ref_aln!r}")
    if not ref_labels or not os.path.exists(ref_labels) or os.path.getsize(ref_labels) == 0:
        raise FileNotFoundError(f"McrA reference labels missing or empty: {ref_labels!r}")
    recs = read_alignment(ref_aln)
    if len(recs) < 3:
        raise ValueError(f"McrA reference alignment has {len(recs)} sequences (need >= 3)")
    ids = [n for n, _ in recs]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate ids in the McrA reference alignment")
    lengths = {len(s) for _, s in recs}
    if len(lengths) != 1:
        raise ValueError(f"McrA reference rows differ in length {sorted(lengths)[:5]}: "
                         "not an alignment")
    raw = load_ref_labels(ref_labels)
    unlabelled = [i for i in ids if i not in raw]
    if unlabelled:
        raise ValueError(f"{len(unlabelled)} reference sequence(s) have no label "
                         f"(e.g. {unlabelled[:5]}); ids must match {ref_labels}")
    groups = {i: norm_group(raw[i]) for i in ids}
    present = set(groups.values())
    if not {"ANME", "methanogen"} <= present:
        raise ValueError("the McrA reference must contain both ANME- and methanogen-labelled "
                         f"sequences (found {sorted(present)})")
    return recs, groups, {i: raw[i] for i in ids}


def run(cmd, log_prefix, **kw):
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        raise RuntimeError(f"{log_prefix} failed (exit {res.returncode}): "
                           f"{(res.stderr or res.stdout).strip()[-600:]}")
    return res


def place_candidates(cands: dict, recs, groups: dict, threads, seed, workdir,
                     min_aligned: int):
    """-> {protein_id: placement dict}. cands: protein_id -> sequence."""
    import dendropy
    import root_tree  # robust midpoint rooting + IQ-TREE name mapping

    for tool in ("mafft", "iqtree2"):
        if shutil.which(tool) is None:
            raise FileNotFoundError(f"{tool} not on PATH (needed for the McrA placement)")
    ref_ids = {n for n, _ in recs}
    clash = sorted(set(cands) & ref_ids)
    if clash:
        raise ValueError(f"candidate ids collide with reference ids: {clash[:5]}")
    os.makedirs(workdir, exist_ok=True)
    ref_path = os.path.join(workdir, "reference.afa")
    cand_path = os.path.join(workdir, "candidates.faa")
    utils.write_fasta(recs, ref_path)
    utils.write_fasta(sorted((k, clean_protein(v)) for k, v in cands.items()), cand_path)
    combined = os.path.join(workdir, "combined.afa")
    with open(combined, "w") as out:
        # --threadit 0: MAFFT's multithreaded iterative refinement differs run to run (MAFFT multithreading docs); this keeps it reproducible
        res = subprocess.run(["mafft", "--add", cand_path, "--keeplength",
                              "--thread", str(threads), "--threadit", "0", ref_path],
                             stdout=out, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"mafft --add failed (exit {res.returncode}): {res.stderr[-600:]}")
    aln = dict(read_alignment(combined))
    width = len(recs[0][1])
    missing = sorted((set(cands) | ref_ids) - set(aln))
    if missing:
        raise RuntimeError(f"mafft output lacks {len(missing)} sequence(s), e.g. {missing[:5]}")
    if {len(s) for s in aln.values()} != {width}:
        raise RuntimeError("mafft --keeplength output is not a fixed-width alignment")

    placements = {}
    aligned = {c: sum(1 for ch in aln[c] if ch not in "-X?") for c in cands}
    placeable = sorted(c for c in cands if aligned[c] >= min_aligned)
    for c in sorted(set(cands) - set(placeable)):
        placements[c] = {"aligned_residues": aligned[c], "clade_ref_groups": "",
                         "mcr_group": "unplaced", "nearest_ref": "",
                         "nearest_ref_group": "", "patristic_distance": ""}
    if not placeable:
        return placements
    tree_in = os.path.join(workdir, "placement.afa")
    utils.write_fasta([(n, s) for n, s in recs] + [(c, aln[c]) for c in placeable], tree_in)
    prefix = os.path.join(workdir, "mcrA")
    run(["iqtree2", "-s", tree_in, "-m", "LG+G", "-fast", "-T", str(threads),
         "--seed", str(seed), "--prefix", prefix, "-redo", "-quiet"], "iqtree2 (McrA)")
    tree = dendropy.Tree.get(path=prefix + ".treefile", schema="newick",
                             preserve_underscores=True)
    back = {}
    for n in list(ref_ids) + placeable:
        back.setdefault(root_tree.iqtree_name(n), []).append(n)
    for lf in tree.leaf_node_iter():
        orig = back.get(lf.taxon.label, [])
        if len(orig) != 1:
            raise RuntimeError(f"cannot map McrA tree tip {lf.taxon.label!r} to one input id")
        lf.taxon.label = orig[0]
    for c, res in classify_tree(tree, groups, set(placeable)).items():
        placements[c] = {"aligned_residues": aligned[c], **res}
    return placements


def classify_tree(tree, groups: dict, cands: set) -> dict:
    """Nearest-labelled-clade assignment on an (unrooted) dendropy tree.

    The tree is midpoint-rooted (root_tree.reroot_midpoint, always bifurcating); for
    each candidate leaf, walk up to the first ancestor whose subtree holds reference
    leaves: one label -> that group, several -> 'ambiguous'. Returns
    {candidate: {clade_ref_groups, mcr_group, nearest_ref, nearest_ref_group,
    patristic_distance}}.
    """
    import root_tree
    if len(tree.seed_node.child_nodes()) == 2:
        tree.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tree.is_rooted = True
    root_tree.reroot_midpoint(tree)
    below = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            lab = nd.taxon.label
            below[id(nd)] = {groups[lab]: 1} if lab in groups else {}
        else:
            acc = {}
            for ch in nd.child_node_iter():
                for g, k in below[id(ch)].items():
                    acc[g] = acc.get(g, 0) + k
            below[id(nd)] = acc
    nearest = nearest_reference(tree, set(groups))
    out = {}
    for lf in tree.leaf_node_iter():
        c = lf.taxon.label
        if c not in cands:
            continue
        node, clade = lf.parent_node, {}
        while node is not None:
            clade = below[id(node)]
            if clade:
                break
            node = node.parent_node
        grp = next(iter(clade)) if len(clade) == 1 else "ambiguous"
        nref, d = nearest.get(c, ("", float("nan")))
        out[c] = {"clade_ref_groups": ",".join(f"{g}:{k}" for g, k in sorted(clade.items())),
                  "mcr_group": grp, "nearest_ref": nref,
                  "nearest_ref_group": groups.get(nref, ""),
                  "patristic_distance": f"{d:.6f}" if nref else ""}
    return out


def nearest_reference(tree, refs: set) -> dict:
    """Patristically nearest reference for every leaf (two-pass, O(n))."""
    INF = (float("inf"), "")
    down = {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            down[id(nd)] = (0.0, nd.taxon.label) if nd.taxon.label in refs else INF
        else:
            best = INF
            for ch in nd.child_node_iter():
                d, r = down[id(ch)]
                cand = (d + (ch.edge.length or 0.0), r)
                if r and (cand < best):
                    best = cand
            down[id(nd)] = best
    up = {id(tree.seed_node): INF}
    for nd in tree.preorder_node_iter():
        kids = nd.child_nodes()
        for ch in kids:
            L = ch.edge.length or 0.0
            best = (up[id(nd)][0] + L, up[id(nd)][1]) if up[id(nd)][1] else INF
            for sib in kids:
                if sib is ch:
                    continue
                d, r = down[id(sib)]
                if r:
                    cand = (d + (sib.edge.length or 0.0) + L, r)
                    if cand < best:
                        best = cand
            up[id(ch)] = best
    out = {}
    for lf in tree.leaf_node_iter():
        a, b = down[id(lf)], up[id(lf)]
        lab = lf.taxon.label
        if lab in refs:
            continue
        best = min(a, b) if a[1] and b[1] else (a if a[1] else b)
        if best[1]:
            out[lab] = (best[1], best[0])
    return out


# --------------------------------------------------------------------------- #
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hits", required=True, help="mcr_hits.tsv (mcr_search.py)")
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--pmap", default="", help="protein_map.tsv: the genome universe")
    ap.add_argument("--ref-aln", default="", help="curated McrA reference alignment")
    ap.add_argument("--ref-labels", default="", help="ref_id<TAB>group TSV")
    ap.add_argument("--required", required=True,
                    help="comma list of KOs or subunits defining a credible Mcr (K00399 = mcrA)")
    ap.add_argument("--anme-tree", required=True, help="True/False: run McrA placement")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--seed", default="12345", help="IQ-TREE seed (reproducibility)")
    ap.add_argument("--min-aligned", type=int, default=50,
                    help="minimum residues in reference columns to place a candidate")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-placement", default="")
    args = ap.parse_args(argv)

    for p in (args.out, args.out_placement):
        if p:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)
    flag = str(args.anme_tree).strip().lower()
    if flag not in ("true", "false", "1", "0", "yes", "no"):
        raise ValueError(f"--anme-tree must be True/False (got {args.anme_tree!r})")
    run_anme = flag in ("true", "1", "yes")
    required = parse_required(args.required)

    per = read_mcr_hits(args.hits)
    universe = sorted(read_universe(args.pmap) | set(per))
    if run_anme:
        recs, groups, raw_labels = validate_reference(args.ref_aln, args.ref_labels)

    credible = {g for g in universe if set(required) <= set(per.get(g, {}))}
    cand_pid = {g: per[g]["mcrA"]["protein_id"] for g in sorted(credible) if "mcrA" in per[g]}
    placements = {}
    if run_anme and cand_pid:
        seqs = fetch_sequences(args.proteome, set(cand_pid.values()))
        lost = sorted(set(cand_pid.values()) - set(seqs))
        if lost:
            raise KeyError(f"{len(lost)} McrA protein(s) from mcr_hits are not in the proteome "
                           f"(e.g. {lost[:5]})")
        placements = place_candidates(seqs, recs, groups, args.threads, args.seed,
                                      args.workdir, args.min_aligned)
    elif run_anme:
        utils.log("No credible McrA in any genome: McrA placement skipped (nothing to place).")

    rows, counts = [], {}
    for g in universe:
        genes = per.get(g, {})
        has_mcrA = "mcrA" in genes
        n_mcr = sum(1 for x in CORE if x in genes)
        ev = []
        a = genes.get("mcrA")
        if a:
            ev.append("mcrA=" + ("contig_edge_fragment" if a["call"] == "partial_rescue"
                                 else "full_threshold"))
            if a["n_copies"] > 1:
                ev.append(f"n_mcrA_copies={a['n_copies']}")
        if genes:
            ev.append("subunits=" + ",".join(x for x in GENES if x in genes))
        else:
            ev.append("no_mcr")
        nref = nref_grp = ""
        if g not in credible:
            grp = "NA"
            if genes:
                miss = [x for x in required if x not in genes]
                ev.append("required_missing=" + ",".join(miss))
        elif not run_anme:
            grp = "not_placed"
            ev.append("anme_check_disabled")
        elif g not in cand_pid:
            grp = "unplaced"
            ev.append("no_mcrA_to_place")
        else:
            pl = placements[cand_pid[g]]
            grp = pl["mcr_group"]
            ev.append(f"mcrA_clade={grp}({pl['clade_ref_groups']})")
            nref, nref_grp = pl["nearest_ref"], pl["nearest_ref_group"]
            if nref:
                ev.append(f"nearest_ref={nref}({nref_grp},d={pl['patristic_distance']})")
        methanogen = g in credible and (not run_anme or grp == "methanogen")
        is_anme = grp == "ANME"
        counts[grp] = counts.get(grp, 0) + 1
        rows.append({
            "genome_id": g, "mcrA": has_mcrA, "n_mcr": n_mcr, "methanogen": methanogen,
            "is_anme": is_anme,
            "methanogenesis_subtype": "unresolved" if methanogen else "NA",
            "evidence": ";".join(ev), "mcr_group": grp,
            "mcrA_protein": a["protein_id"] if a else "",
            "mcrA_call": a["call"] if a else "",
            "n_mcrA_copies": a["n_copies"] if a else 0,
            "nearest_ref": nref, "nearest_ref_group": nref_grp,
        })

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in HEADER_COLS) + "\n")
    if args.out_placement:
        with open(args.out_placement, "w") as out:
            out.write("\t".join(PLACEMENT_COLS) + "\n")
            for g in sorted(cand_pid):
                pl = placements.get(cand_pid[g])
                if pl is None:
                    continue
                vals = {"genome_id": g, "protein_id": cand_pid[g], **pl}
                out.write("\t".join(str(vals[c]) for c in PLACEMENT_COLS) + "\n")
    n_meth = sum(1 for r in rows if r["methanogen"])
    utils.log(f"methanogen_calls: {len(rows)} genomes; {len(credible)} with credible Mcr; "
              f"{n_meth} methanogen; groups {dict(sorted(counts.items()))}; "
              f"required={required}; placement={'on' if run_anme else 'off'}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
