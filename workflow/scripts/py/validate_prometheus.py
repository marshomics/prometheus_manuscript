"""Integrate detection, domain architecture, phylogenetic separation, defense
context, length sanity and fragmentation into a confidence-tiered Prometheus
catalog (Stage 1.6 of the design).

Scientific framing
-------------------
A DefenseFinder Prometheus call is *not* accepted at face value. Prometheus is a
single ~2221-aa DNA2-family helicase-nuclease whose four/five canonical domains
(DUF4011, DUF3320, DNA2 helicase/nuclease, REase MTES) must be present in the
right N->C order, and which must place phylogenetically *away* from the cellular
DNA2 / NurA-HerA paralogs that share the helicase-nuclease core. This script
combines the available genomic evidence into per-locus confidence tiers (A/B/C)
plus an explicit Rejected class, mirroring §1.6/§1.115 of the design document.

Important method caveats (documented, not hidden):
  * CATALYTIC STATUS here is a *domain-architecture proxy*, not a residue-level
    check. A defensible catalytic call needs an MSA against the characterised
    Pseudomonas Prometheus and inspection of the nuclease active site, the
    Walker A/B helicase motifs and the REase catalytic residues (design §1.3).
    That requires the family alignment built in Stage 5, so here we only state
    whether the catalytic *modules* (DNA2 + a nuclease/REase) are present.
  * DEFENSE CONTEXT here is a *genome-level proxy* (does the genome carry other
    defense systems / >=2 systems total), not a true ±10-ORF neighbourhood test
    (design §1.5/§1.8). The neighbourhood test needs per-ORF coordinates that
    are not available at this step; the genome-level signal is conservative and
    only ever *raises* confidence, never rejects.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import json
import pandas as pd

CATALOG_HEADER = (
    "locus_id\tgenome_id\tprotein_id\tlength_aa\tdetect_defensefinder\t"
    "architecture\tn_domains\tcatalytic_status\tdefense_context\t"
    "tree_assignment\tfragmented\tconfidence_tier\treason"
)
SENSITIVITY_HEADER = "genome_id\tdetect_defensefinder\tprometheus_present"


# --------------------------------------------------------------------------- #
# Domain architecture
# --------------------------------------------------------------------------- #
def domain_matches(hit_name: str, hit_acc: str, spec: dict) -> bool:
    """True if an hmmscan hit matches a configured domain.

    Matching is deliberately permissive and release-robust (design §1.2 note):
      * by NAME: case-insensitive substring match against any configured name
        (Pfam clan/family names drift between releases, so substring not equals);
      * by ACCESSION: exact match against any configured accession (e.g. PF13195).
    Accession matching ignores the Pfam version suffix (PF13195.7 -> PF13195).
    """
    names = [str(n).lower() for n in spec.get("names", [])]
    accs = [str(a) for a in spec.get("accessions", [])]
    name_l = str(hit_name).lower()
    acc_base = str(hit_acc).split(".")[0]
    for n in names:
        if n and n in name_l:
            return True
    for a in accs:
        if a and (a == hit_acc or a.split(".")[0] == acc_base):
            return True
    return False


def architecture_for(hits: list, domains: dict):
    """Return (architecture_string, n_domains, present_set) for one protein.

    `hits` are the domtbl hit dicts for this protein (utils.parse_domtbl).
    Present domains are ordered N->C by alignment start coordinate. When two
    different domains overlap on the same residues, the lower start wins; ties
    are broken by the configured domain order so output is deterministic.
    """
    domain_order = list(domains.keys())
    # earliest ali start at which each configured domain is observed
    first_start = {}
    for h in hits:
        for dname, spec in domains.items():
            if domain_matches(h.get("domain", ""), h.get("accession", ""), spec):
                s = h.get("start", 10**9)
                if dname not in first_start or s < first_start[dname]:
                    first_start[dname] = s
    present = sorted(
        first_start.keys(),
        key=lambda d: (first_start[d], domain_order.index(d)),
    )
    arch = "-".join(present)
    return arch, len(present), set(present)


def catalytic_status(present: set, n_domains: int, has_domtbl: bool) -> str:
    """Domain-architecture proxy for catalytic potential (NOT residue-level).

    intact     : DNA2 module present AND (nuclease OR REase) present AND >=3 domains
    partial    : 1<=n_domains<3, OR the DNA2 module is missing (helicase-less)
    degenerate : n_domains<=1 (essentially a lone domain)
    unknown    : no domain information available for this protein
    See module docstring for why a true catalytic call needs an MSA.
    """
    if not has_domtbl:
        return "unknown"
    has_dna2 = "DNA2" in present
    has_cut = ("nuclease" in present) or ("REase" in present)
    if n_domains <= 1:
        return "degenerate"
    if has_dna2 and has_cut and n_domains >= 3:
        return "intact"
    if (1 <= n_domains < 3) or (not has_dna2):
        return "partial"
    # n_domains>=3 but missing a catalytic module -> still partial (architecture
    # present but catalytic potential not supported).
    return "partial"


# --------------------------------------------------------------------------- #
# Phylogenetic separation
# --------------------------------------------------------------------------- #
def assign_by_tree(tree_path: str, ref_labels: dict):
    """Map each candidate tip to 'Prometheus'/'Paralog' by nearest reference tip.

    Uses dendropy's patristic (phylogenetic) distance matrix: for every tip in
    the tree that is NOT a reference, find the reference tip at minimum patristic
    distance and inherit its group. Candidates not in the tree, or trees with
    <2 reference tips, yield 'Unplaced' for everything (the validation tree could
    not separate the clades). dendropy is imported lazily so this module imports
    even where dendropy is absent (e.g. during py_compile / unrelated tests).

    Returns dict {candidate_tip -> 'Prometheus'|'Paralog'|'Unplaced'}.
    """
    assignments = {}
    if not tree_path or not os.path.exists(tree_path) or os.path.getsize(tree_path) == 0:
        return assignments
    try:
        import dendropy
    except Exception:
        utils.log("WARNING: dendropy unavailable; all candidates -> Unplaced")
        return assignments
    try:
        tree = dendropy.Tree.get(path=tree_path, schema="newick",
                                 preserve_underscores=True)
    except Exception as exc:
        utils.log(f"WARNING: could not parse tree ({exc}); all -> Unplaced")
        return assignments

    taxa = [t for t in tree.taxon_namespace if t.label is not None]
    tip_labels = [t.label for t in taxa]
    # reference tips actually present in the tree, by group
    ref_present = {lbl: ref_labels[lbl] for lbl in tip_labels if lbl in ref_labels}
    if len(ref_present) < 2:
        utils.log("Validation tree has <2 reference tips; candidates -> Unplaced")
        return assignments

    try:
        pdm = tree.phylogenetic_distance_matrix()
    except Exception as exc:
        utils.log(f"WARNING: distance matrix failed ({exc}); all -> Unplaced")
        return assignments

    taxon_by_label = {t.label: t for t in taxa}
    for lbl in tip_labels:
        if lbl in ref_labels:
            continue  # references are not candidates
        cand_taxon = taxon_by_label[lbl]
        best_group, best_dist = None, None
        for ref_lbl, group in ref_present.items():
            ref_taxon = taxon_by_label[ref_lbl]
            try:
                d = pdm.patristic_distance(cand_taxon, ref_taxon)
            except Exception:
                continue
            if d is None:
                continue
            if best_dist is None or d < best_dist:
                best_dist, best_group = d, group
        if best_group == "Prometheus_ref":
            assignments[lbl] = "Prometheus"
        elif best_group == "Paralog_ref":
            assignments[lbl] = "Paralog"
        else:
            assignments[lbl] = "Unplaced"
    return assignments


def load_ref_labels(path: str) -> dict:
    """ref tip name -> 'Prometheus_ref' | 'Paralog_ref' (TSV, no header assumed)."""
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
            # tolerate a header row
            if grp in ("Prometheus_ref", "Paralog_ref"):
                labels[name] = grp
    return labels


# --------------------------------------------------------------------------- #
# Defense context (genome-level proxy)
# --------------------------------------------------------------------------- #
def defense_context_for_genomes(df_dir: str, genomes: set) -> dict:
    """Genome-level defense-island proxy (design §1.5/§1.8 approximation).

    True if a genome carries >=2 DefenseFinder systems total OR >=1 system other
    than Prometheus. Counts data rows of {df_dir}/{g}/{g}_defense_finder_systems.tsv.
    This is conservative: it can only raise confidence, never reject, and it
    stands in for the real ±10-ORF neighbourhood test which needs ORF coordinates.
    """
    out = {}
    for g in genomes:
        n_total, n_other = 0, 0
        f = os.path.join(df_dir, g, f"{g}_defense_finder_systems.tsv")
        if os.path.exists(f):
            try:
                df = pd.read_csv(f, sep="\t")
                df.columns = [c.strip() for c in df.columns]
                n_total = len(df)
                if "type" in df.columns:
                    for t in df["type"].dropna().astype(str):
                        if "prometheus" not in t.lower():
                            n_other += 1
            except Exception:
                pass
        out[g] = (n_total >= 2) or (n_other >= 1)
    return out


# --------------------------------------------------------------------------- #
# Tiering
# --------------------------------------------------------------------------- #
def assign_tier(detect_df, n_domains, tree_assignment, fragmented, catalytic,
                defense_context, length_aa, min_len, max_len):
    """Return (confidence_tier, reason).

    Detection is DefenseFinder-only, so tier A cannot rest on multi-method
    agreement. Instead A requires the strongest standalone evidence: the full
    canonical four-domain architecture, placement with the Prometheus references
    (not the cellular DNA2/NurA-HerA paralogs), intact catalytic modules, and a
    full-length (non-fragmented) ORF. B is strong-but-incomplete (Prometheus
    placement OR architecture+defense context); C is retained-but-flagged.

    Order of checks: Rejected gates first (paralog placement, or length far
    outside the sanity window), then A, then B, else C.
    """
    reasons = []
    # ---- Rejected gates -------------------------------------------------- #
    if tree_assignment == "Paralog":
        return "Rejected", "placed_with_cellular_paralog"
    # length far outside sanity bounds: a ~2.2 kb protein cannot be <0.5*min_len,
    # and a hit several-fold over max_len is a fusion/mis-call, not Prometheus.
    if length_aa is not None and length_aa > 0:
        if length_aa < 0.5 * min_len:
            return "Rejected", f"length_far_below_min({length_aa}<0.5*{min_len})"
        if length_aa > 2.0 * max_len:
            return "Rejected", f"length_far_above_max({length_aa}>2*{max_len})"

    # ---- Tier A ---------------------------------------------------------- #
    # Full canonical architecture (>=4 of the 5 configured domains) + correct
    # placement + intact catalytic modules + full length.
    if (detect_df and n_domains >= 4 and tree_assignment == "Prometheus"
            and not fragmented and catalytic == "intact"):
        return "A", ""

    # ---- Tier B ---------------------------------------------------------- #
    if tree_assignment == "Prometheus":
        reasons.append("prometheus_placement")
    if n_domains >= 3 and defense_context:
        reasons.append("architecture+defense_context")
    if reasons:
        if fragmented:
            reasons.append("fragmented")
        return "B", ";".join(reasons)

    # ---- Tier C (retained, labelled) ------------------------------------ #
    c_reasons = []
    if tree_assignment == "Unplaced":
        c_reasons.append("unplaced")
    if catalytic in ("partial", "degenerate"):
        c_reasons.append(catalytic)
    if n_domains < 3:
        c_reasons.append(f"n_domains={n_domains}")
    if fragmented:
        c_reasons.append("fragmented")
    if not c_reasons:
        c_reasons.append("divergent_or_partial_evidence")
    return "C", ";".join(c_reasons)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sources", required=True, help="candidate_sources.tsv")
    ap.add_argument("--raw", required=True, help="raw_hits.tsv")
    ap.add_argument("--domtbl", required=True, help="hmmscan domtblout vs Pfam")
    ap.add_argument("--tree", required=True, help="validation tree (newick)")
    ap.add_argument("--ref-labels", required=True, help="reference tip -> group TSV")
    ap.add_argument("--candidates", required=True, help="candidate proteins FASTA")
    ap.add_argument("--meta", required=True, help="genome_metadata.tsv")
    ap.add_argument("--reps", required=True, help="representatives.txt")
    ap.add_argument("--df-dir", required=True, help="DefenseFinder output root")
    ap.add_argument("--context-window", type=int, default=10)
    ap.add_argument("--domains", required=True, help="JSON domain definitions")
    ap.add_argument("--min-len", type=int, required=True)
    ap.add_argument("--max-len", type=int, required=True)
    ap.add_argument("--out-catalog", required=True)
    ap.add_argument("--out-proteins", required=True)
    ap.add_argument("--out-sensitivity", required=True)
    args = ap.parse_args()

    for p in (args.out_catalog, args.out_proteins, args.out_sensitivity):
        os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)

    try:
        domains = json.loads(args.domains) if args.domains else {}
    except json.JSONDecodeError as exc:
        utils.log(f"WARNING: could not parse --domains JSON ({exc}); using empty set")
        domains = {}
    domain_names = list(domains.keys())

    # ---- genome universe (from metadata) for the sensitivity table ------- #
    meta_genomes = []
    try:
        meta = pd.read_csv(args.meta, sep="\t")
        if "genome_id" in meta.columns:
            meta_genomes = [str(g) for g in meta["genome_id"].tolist()]
    except Exception:
        utils.log("WARNING: could not read --meta; sensitivity limited to candidate genomes")

    # ---- candidate sources (detection flags) ----------------------------- #
    try:
        sources = pd.read_csv(args.sources, sep="\t")
    except Exception:
        sources = pd.DataFrame()

    if sources.empty or "protein_id" not in sources.columns:
        # No candidates at all -> header-only catalog/proteins, but still emit a
        # sensitivity row per known genome (all absent) so downstream joins work.
        utils.write_empty(args.out_catalog, CATALOG_HEADER)
        open(args.out_proteins, "w").close()
        if meta_genomes:
            rows = [{
                "genome_id": g, "detect_defensefinder": 0,
                "prometheus_present": False,
            } for g in meta_genomes]
            pd.DataFrame(rows, columns=SENSITIVITY_HEADER.split("\t")).to_csv(
                args.out_sensitivity, sep="\t", index=False)
        else:
            utils.write_empty(args.out_sensitivity, SENSITIVITY_HEADER)
        utils.log("No candidate sources; wrote header-only catalog.")
        return

    def as_int01(v):
        try:
            return 1 if int(float(v)) != 0 else 0
        except (ValueError, TypeError):
            return 0

    # ---- sequences (for length_aa and the proteins FASTA) ---------------- #
    seqs = dict(utils.read_fasta(args.candidates))

    # ---- domain hits ----------------------------------------------------- #
    domtbl_hits = utils.parse_domtbl(args.domtbl)

    # ---- phylogenetic placement ----------------------------------------- #
    ref_labels = load_ref_labels(args.ref_labels)
    tree_assign = assign_by_tree(args.tree, ref_labels)

    # ---- defense context per genome ------------------------------------- #
    candidate_genomes = set(utils.genome_of(str(p)) for p in sources["protein_id"])
    ctx = defense_context_for_genomes(args.df_dir, candidate_genomes)

    # ---- per-locus catalog ---------------------------------------------- #
    rows = []
    proteins_out = []  # (locus_id, seq) for non-Rejected loci
    # per-genome accumulator: genome -> [detect_defensefinder(0/1), present(bool)]
    # (DefenseFinder is the sole detector; no padloc/hmm columns in the schema).
    sens = {}

    for _, r in sources.iterrows():
        pid = str(r["protein_id"])
        gid = utils.genome_of(pid)
        d_df = as_int01(r.get("detect_defensefinder", 0))
        detected = d_df >= 1

        seq = seqs.get(pid, "")
        length_aa = len(seq) if seq else 0

        hits = domtbl_hits.get(pid, [])
        arch, n_domains, present = architecture_for(hits, domains)
        # Catalytic status is 'unknown' for a protein with NO domain hits of its
        # own (conservative): a global domtbl with hits for OTHER proteins does
        # not let us assert anything about this one. So key off this protein's
        # own hits, not has_any_domtbl.
        cat = catalytic_status(present, n_domains, has_domtbl=(pid in domtbl_hits))

        assignment = tree_assign.get(pid, "Unplaced")

        # FRAGMENTED: a long (~2.2 kb) protein split by assembly -> short ORF
        # that was nevertheless detected. Only meaningful for detected loci.
        fragmented = bool(detected and length_aa > 0 and length_aa < args.min_len)

        defense_ctx = bool(ctx.get(gid, False))

        tier, reason = assign_tier(
            d_df, n_domains, assignment, fragmented, cat, defense_ctx,
            length_aa if length_aa > 0 else None, args.min_len, args.max_len)

        rows.append({
            "locus_id": pid, "genome_id": gid, "protein_id": pid,
            "length_aa": length_aa,
            "detect_defensefinder": d_df,
            "architecture": arch, "n_domains": n_domains,
            "catalytic_status": cat, "defense_context": defense_ctx,
            "tree_assignment": assignment, "fragmented": fragmented,
            "confidence_tier": tier, "reason": reason,
        })

        if tier != "Rejected":
            if seq:
                proteins_out.append((pid, seq))
            # a genome counts as 'present' only on a non-Rejected locus
            s = sens.setdefault(gid, [0, False])
            s[1] = True

        # DefenseFinder detection flag contributes regardless of tier
        s = sens.setdefault(gid, [0, False])
        s[0] = max(s[0], d_df)

    # ---- write catalog --------------------------------------------------- #
    if rows:
        pd.DataFrame(rows, columns=CATALOG_HEADER.split("\t")).to_csv(
            args.out_catalog, sep="\t", index=False)
    else:
        utils.write_empty(args.out_catalog, CATALOG_HEADER)

    # ---- write proteins FASTA (non-Rejected) ---------------------------- #
    if proteins_out:
        utils.write_fasta(proteins_out, args.out_proteins)
    else:
        open(args.out_proteins, "w").close()

    # ---- write detection sensitivity (every known genome) --------------- #
    all_genomes = list(dict.fromkeys(meta_genomes + sorted(candidate_genomes)))
    sens_rows = []
    for g in all_genomes:
        d_df, present = sens.get(g, [0, False])
        sens_rows.append({
            "genome_id": g,
            "detect_defensefinder": d_df,
            "prometheus_present": bool(present),
        })
    if sens_rows:
        pd.DataFrame(sens_rows, columns=SENSITIVITY_HEADER.split("\t")).to_csv(
            args.out_sensitivity, sep="\t", index=False)
    else:
        utils.write_empty(args.out_sensitivity, SENSITIVITY_HEADER)

    n_rej = sum(1 for r in rows if r["confidence_tier"] == "Rejected")
    utils.log(
        f"Catalog: {len(rows)} loci ({len(rows) - n_rej} retained, {n_rej} rejected); "
        f"{len(proteins_out)} proteins written; {len(sens_rows)} genomes in sensitivity."
    )


if __name__ == "__main__":
    main()
