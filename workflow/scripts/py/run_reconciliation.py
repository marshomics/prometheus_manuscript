"""Species-tree-aware DTL reconciliation of the Prometheus gene family (GeneRax; optional ALE).

Stage 5.3. Verified against the GeneRax 2.1.3 source (src/GeneRax/GeneRaxArguments.cpp,
ext/GeneRaxCore@8487673: IO/FamiliesFileParser.cpp, IO/GeneSpeciesMapping.cpp,
IO/Families.cpp, routines/Routines.cpp, util/Scenario.cpp) and by running it:

  * CLI: ``generax -f FAMILIES -s SPECIES_TREE -r UndatedDTL -p PREFIX --strategy SPR|EVAL
    --seed N`` (valid -r: UndatedDL, UndatedDTL, SimpleDS, ParsimonyD; 'Auto' is NOT valid
    in 2.1.3). The gene-tree root is optimised by default (rooted reconciliation likelihood
    with the root searched; --enforce-gene-tree-root off), so the unrooted IQ-TREE tree is
    the correct input. The species tree is read as ROOTED.
  * families file: ``[FAMILIES]`` / ``- name`` / ``starting_gene_tree = ...`` /
    ``alignment = ...`` / ``mapping = ...`` / ``subst_model = ...`` (whitespace stripped).
  * mapping: Treerecs format ``gene species`` per line (Phyldog format is chosen when the
    first line contains ':').
  * GeneRax DISCARDS A FAMILY (with a log warning, exit 0) when any gene maps to a species
    absent from the species tree, when alignment and gene-tree labels differ, or when the
    gene tree is not binary -- genes whose genome is not a species-tree tip are therefore
    pruned from the gene tree AND the alignment here, and the count is reported.
  * outputs: PREFIX/reconciliations/<family>_eventCounts.txt with lines S:, SL:, D:, T:,
    TL:, L:, Leaf: (Scenario::saveEventsCounts); <family>_speciesEventCounts.txt
    ('species_label, speciations, duplications, losses, transfers, presence, origination');
    <family>_transfers.txt ('donor recipient 1'). Per-species output dereferences every
    species-node label, so all internal nodes get unique labels before GeneRax runs
    (IQ-TREE support values are not labels).

Event accounting written to reconciliation_events.tsv (one row per event type per tool):
  Speciation = S + SL   (vertical transmissions through species-tree nodes)
  Transfer   = T + TL   (horizontal transfers; TL = transfer with loss of the donor copy)
  Duplication = D;  Loss = L + SL + TL;  Origination = 1 (node in the note)
  fraction (Speciation and Transfer rows) = VERTICAL fraction
            (S + SL) / (S + SL + T + TL); the Transfer row note also gives the transfer fraction.
ALE (optional, ALEobserve + ALEml_undated) re-estimates from the UFBoot distribution:
Total line of '<species>_<ale>.uml_rec' (D T L S means over sampled reconciliations).

NO SILENT FAILURE: a requested tool that is missing or fails, or a reconciliation whose
Leaf count differs from the number of genes, raises. Count-0 rows with a note are written
only when reconciliation is disabled or the trees are genuinely too small.
"""
import argparse
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

HEADER = "tool\tevent\tcount\tfraction\tnote"
EVENT_TYPES = ["Duplication", "Transfer", "Loss", "Speciation", "Origination"]
TRANSFER_HEADER = "tool\tdonor\trecipient\tcount\trecipient_methanogen_clade\trecipient_n_genomes"
MIN_TIPS = 4
GENERAX_MODELS = ("UndatedDL", "UndatedDTL", "SimpleDS", "ParsimonyD")
FAMILY = "prometheus"


def fail(msg):
    raise SystemExit(f"[run_reconciliation] {msg}")


def truthy(s):
    return str(s).strip().lower() in ("true", "1", "yes", "t")


def is_trivial(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return True
    with open(path) as fh:
        return fh.read().strip() in ("", "();", ";")


def read_tree(path, rooting):
    import dendropy
    return dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True, rooting=rooting)


def newick(tree, tip=lambda nd: nd.taxon.label, inner=lambda nd: "", lengths=True):
    def lab(nd):
        s = tip(nd) if nd.is_leaf() else inner(nd)
        bl = nd.edge.length
        return s + (f":{bl:.10g}" if lengths and bl is not None and nd.parent_node is not None else "")

    def rec(nd):
        if nd.is_leaf():
            return lab(nd)
        return "(" + ",".join(rec(c) for c in nd.child_nodes()) + ")" + lab(nd)
    return rec(tree.seed_node) + ";"


def species_tree_for_tools(path, out_path):
    """Rooted species tree with unique internal labels sp_n<k>; returns (tree, tip set, labels)."""
    t = read_tree(path, "force-rooted")
    root_children = t.seed_node.child_nodes()
    if len(root_children) != 2:
        fail(f"species tree {path} is not rooted (root has {len(root_children)} children); "
             "UndatedDTL needs the rooted tree from scaffolds/species_tree")
    poly = [nd for nd in t.preorder_internal_node_iter() if len(nd.child_nodes()) != 2]
    if poly:
        fail(f"species tree {path} is not binary ({len(poly)} internal node(s) without exactly "
             "two children); GeneRax/ALE need a rooted binary species tree")
    tips = [lf.taxon.label for lf in t.leaf_node_iter()]
    if len(set(tips)) != len(tips):
        fail("species tree has duplicated tip labels")
    # IQ-TREE 'SH-aLRT/UFBoot' support values arrive as internal node labels; GeneRax uses
    # internal labels as species-node NAMES (duplicates merge per-species counts, missing
    # ones are dereferenced), so every internal node gets a unique name instead.
    k = 0
    for nd in t.preorder_internal_node_iter():
        k += 1
        nd.label = f"sp_n{k}"
    with open(out_path, "w") as fh:
        fh.write(newick(t, inner=lambda nd: nd.label) + "\n")
    return t, set(tips)


def methanogen_clades(stree, master_path):
    """species-node label -> (all descendant tips methanogen?, n descendant genomes)."""
    import pandas as pd
    out = {}
    if not master_path:
        return out
    df = pd.read_csv(master_path, sep="\t", dtype=str, keep_default_na=False)
    if not {"genome_id", "methanogen"} <= set(df.columns):
        fail(f"master table {master_path} lacks genome_id/methanogen")
    meth = {g.strip(): m.strip().lower() in ("true", "1", "yes", "t")
            for g, m in zip(df["genome_id"], df["methanogen"])}
    for nd in stree.postorder_node_iter():
        tips = [lf.taxon.label for lf in nd.leaf_iter()]
        lab = nd.taxon.label if nd.is_leaf() else nd.label
        out[lab] = (all(meth.get(t, False) for t in tips), len(tips))
    return out


def prepare_family(args, sp_tips, workdir):
    """Prune genes whose genome is not a species-tree tip; write tree/aln/mapping/families.

    Returns (n_genes_kept, n_pruned, gene_tree_path, aln_path).
    """
    gt = read_tree(args.gene_tree, "force-unrooted")
    tips = [lf.taxon.label for lf in gt.leaf_node_iter()]
    if len(set(tips)) != len(tips):
        fail("gene tree has duplicated tip labels")
    keep = [t for t in tips if utils.genome_of(t) in sp_tips]
    pruned = sorted(set(tips) - set(keep))
    if pruned:
        utils.log(f"WARNING: {len(pruned)} gene(s) map to genomes absent from the species tree and "
                  f"are pruned (GeneRax would otherwise discard the whole family): {pruned[:10]}")
        gt.prune_taxa_with_labels(pruned, suppress_unifurcations=True)
    if len(gt.seed_node.child_nodes()) == 2 and len(keep) > 2:
        gt.collapse_basal_bifurcation(set_as_unrooted_tree=True)   # GeneRax reads it unrooted
    aln = dict(utils.read_fasta(args.aln))
    missing = [t for t in keep if t not in aln]
    if missing:
        fail(f"{len(missing)} gene-tree tip(s) absent from the alignment {args.aln}: {missing[:5]}")
    tree_path = os.path.join(workdir, "gene_tree.nwk")
    with open(tree_path, "w") as fh:
        fh.write(newick(gt) + "\n")
    aln_path = os.path.join(workdir, "gene_alignment.faa")
    utils.write_fasta([(t, aln[t]) for t in keep], aln_path)
    with open(os.path.join(workdir, "mapping.link"), "w") as fh:
        for t in keep:
            fh.write(f"{t} {utils.genome_of(t)}\n")
    return len(keep), len(pruned), tree_path, aln_path


def write_families(path, gene_tree, aln, mapping, subst_model):
    with open(path, "w") as out:
        out.write("[FAMILIES]\n")
        out.write(f"- {FAMILY}\n")
        out.write(f"starting_gene_tree = {os.path.abspath(gene_tree)}\n")
        out.write(f"alignment = {os.path.abspath(aln)}\n")
        out.write(f"mapping = {os.path.abspath(mapping)}\n")
        out.write(f"subst_model = {subst_model}\n")


def parse_event_counts(path):
    """GeneRax <family>_eventCounts.txt -> {S, SL, D, T, TL, L, Leaf}."""
    if not os.path.exists(path):
        fail(f"GeneRax produced no {path} (family discarded? see the GeneRax log)")
    counts = {}
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*([A-Za-z]+)\s*:\s*([0-9.]+)\s*$", line)
            if m:
                counts[m.group(1)] = float(m.group(2))
    need = ("S", "SL", "D", "T", "TL", "L", "Leaf")
    if any(k not in counts for k in need):
        fail(f"{path} lacks one of {need}: {counts}")
    return counts


def parse_species_events(path):
    """origination species label from <family>_speciesEventCounts.txt (or '')."""
    if not os.path.exists(path):
        return ""
    with open(path) as fh:
        header = fh.readline()
        if "origination" not in header:
            return ""
        for line in fh:
            f = [x.strip() for x in line.split(",")]
            if len(f) >= 7 and float(f[6] or 0) > 0:
                return f[0]
    return ""


def event_rows(tool, D, T, L, S, note, origin=""):
    denom = S + T
    vert = round(S / denom, 6) if denom > 0 else ""
    tfrac = round(T / denom, 6) if denom > 0 else ""
    num = lambda x: int(x) if float(x).is_integer() else round(x, 4)  # noqa: E731
    return [
        (tool, "Duplication", num(D), "", note),
        (tool, "Transfer", num(T), vert, f"{note}; transfer_fraction=T/(S+T)={tfrac}"),
        (tool, "Loss", num(L), "", note),
        (tool, "Speciation", num(S), vert, f"{note}; vertical_fraction=S/(S+T)={vert}"),
        (tool, "Origination", 1, "", f"{note}; origination species node={origin or 'NA'}"),
    ]


def run_generax(args, sp_path, sp_tree, sp_tips, rows, transfers):
    generax = shutil.which("generax")
    if generax is None:
        fail("generax not found on PATH (evolution.reconciliation.run_generax is true)")
    if args.model not in GENERAX_MODELS:
        fail(f"invalid GeneRax --rec-model '{args.model}' (2.1.3 accepts {GENERAX_MODELS})")
    wd = os.path.join(args.workdir, "generax_input")
    os.makedirs(wd, exist_ok=True)
    n_genes, n_pruned, gtree, aln = prepare_family(args, sp_tips, wd)
    if n_genes < MIN_TIPS:
        rows.extend((("generax", e, 0, "", f"too few genes on the species tree ({n_genes} < {MIN_TIPS}; "
                      f"{n_pruned} pruned)")) for e in EVENT_TYPES)
        return
    fam = os.path.join(wd, "families.txt")
    write_families(fam, gtree, aln, os.path.join(wd, "mapping.link"), args.subst_model)
    prefix = os.path.join(args.workdir, "generax")
    if os.path.isdir(prefix):
        shutil.rmtree(prefix)          # GeneRax resumes from checkpoints in an existing prefix
    cmd = [generax, "--families", fam, "--species-tree", sp_path, "--rec-model", args.model,
           "--prefix", prefix, "--strategy", args.strategy, "--seed", str(args.seed)]
    utils.log("running: " + " ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(os.path.join(args.workdir, "generax.stdout.log"), "w") as fh:
        fh.write(proc.stdout)
    if proc.returncode != 0:
        fail(f"generax exited {proc.returncode}: {proc.stdout[-2000:]}")
    # GeneRax always logs "Filtering invalid families..."; a DISCARDED family is reported as
    # "Error in family <name>: <reason>" and "WARNING!!! Found <n> invalid families".
    discard = [ln for ln in proc.stdout.splitlines()
               if re.search(r"Error in family|Found [1-9][0-9]* invalid famil", ln)]
    if discard:
        fail("GeneRax discarded the family: " + " | ".join(discard)[:1500])
    rdir = os.path.join(prefix, "reconciliations")
    c = parse_event_counts(os.path.join(rdir, f"{FAMILY}_eventCounts.txt"))
    if int(c["Leaf"]) != n_genes:
        fail(f"GeneRax reconciliation is inconsistent: Leaf={c['Leaf']:g} but {n_genes} genes "
             "(empty or truncated scenario)")
    origin = parse_species_events(os.path.join(rdir, f"{FAMILY}_speciesEventCounts.txt"))
    S, T = c["S"] + c["SL"], c["T"] + c["TL"]
    D, L = c["D"], c["L"] + c["SL"] + c["TL"]
    note = (f"GeneRax {args.model} --strategy {args.strategy}; {n_genes} genes ({n_pruned} pruned: genome "
            f"not in species tree); raw S={c['S']:g} SL={c['SL']:g} D={c['D']:g} T={c['T']:g} "
            f"TL={c['TL']:g} L={c['L']:g}")
    rows.extend(event_rows("generax", D, T, L, S, note, origin))
    tf = os.path.join(rdir, f"{FAMILY}_transfers.txt")
    if os.path.exists(tf):
        pairs = {}
        with open(tf) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 2:
                    key = (f[0], f[1])
                    pairs[key] = pairs.get(key, 0) + (float(f[2]) if len(f) > 2 else 1.0)
        mc = methanogen_clades(sp_tree, args.master)
        for (d, r), n in sorted(pairs.items()):
            is_m, ng = mc.get(r, ("", ""))
            transfers.append(("generax", d, r, n, is_m, ng))


def parse_uml_rec(path):
    """ALE '<species>_<ale>.uml_rec' -> (D, T, L, S) from the 'Total' line under '# of'."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    for i, line in enumerate(lines):
        if line.startswith("# of") and "Speciations" in line:
            hdr = [h.strip() for h in line.split("\t")[1:]]
            if i + 1 < len(lines) and lines[i + 1].startswith("Total"):
                vals = [float(x) for x in lines[i + 1].split("\t")[1:] if x.strip()]
                d = dict(zip(hdr, vals))
                try:
                    return d["Duplications"], d["Transfers"], d["Losses"], d["Speciations"]
                except KeyError:
                    break
    fail(f"could not parse the 'Total' event line of {path}")


def run_ale(args, sp_tree, sp_tips, rows):
    aleobs, aleml = shutil.which("ALEobserve"), shutil.which("ALEml_undated")
    if aleobs is None or aleml is None:
        fail("ALEobserve/ALEml_undated not on PATH (evolution.reconciliation.run_ale is true)")
    if is_trivial(args.ufboot):
        rows.extend(("ale", e, 0, "", "too few trees: UFBoot file trivial/empty") for e in EVENT_TYPES)
        return
    import dendropy
    adir = os.path.join(args.workdir, "ale")
    os.makedirs(adir, exist_ok=True)
    # ALE maps genes to species by splitting the gene name at a separator CHARACTER; genome ids
    # may contain '_' themselves, so species are renamed S<k> and genes S<k>.<j>.
    sp_name = {t: f"S{i}" for i, t in enumerate(sorted(sp_tips), 1)}
    st = dendropy.Tree(sp_tree)
    for lf in st.leaf_node_iter():
        lf.taxon.label = sp_name[lf.taxon.label]
    sp_file = os.path.join(adir, "species.nwk")
    with open(sp_file, "w") as fh:
        fh.write(newick(st, tip=lambda nd: nd.taxon.label, inner=lambda nd: "") + "\n")
    trees = dendropy.TreeList.get(path=args.ufboot, schema="newick", preserve_underscores=True,
                                  rooting="force-unrooted")
    gene_name, n_pruned = {}, 0
    out_trees = []
    for t in trees:
        tips = [lf.taxon.label for lf in t.leaf_node_iter()]
        drop = [x for x in tips if utils.genome_of(x) not in sp_tips]
        n_pruned = max(n_pruned, len(drop))
        if drop:
            t.prune_taxa_with_labels(drop, suppress_unifurcations=True)
        if len(t.seed_node.child_nodes()) == 2 and len(tips) - len(drop) > 2:
            t.collapse_basal_bifurcation(set_as_unrooted_tree=True)
        for lf in t.leaf_node_iter():
            g = lf.taxon.label
            if g not in gene_name:
                sp = sp_name[utils.genome_of(g)]
                gene_name[g] = f"{sp}.{sum(1 for v in gene_name.values() if v.startswith(sp + '.')) + 1}"
        out_trees.append(newick(t, tip=lambda nd: gene_name[nd.taxon.label]))
    if len(gene_name) < MIN_TIPS:
        rows.extend(("ale", e, 0, "", f"too few genes on the species tree ({len(gene_name)})") for e in EVENT_TYPES)
        return
    uf = os.path.join(adir, "gene.ufboot")
    with open(uf, "w") as fh:
        fh.write("\n".join(out_trees) + "\n")
    with open(os.path.join(adir, "gene_name_map.tsv"), "w") as fh:
        fh.write("gene_id\tale_name\n" + "".join(f"{g}\t{a}\n" for g, a in gene_name.items()))
    for cmd in ([aleobs, uf], [aleml, sp_file, uf + ".ale", f"sample={args.ale_samples}",
                                f"seed={args.seed}", "separators=."]):
        utils.log("running: " + " ".join(cmd))
        p = subprocess.run(cmd, cwd=adir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.returncode != 0:
            fail(f"{os.path.basename(cmd[0])} exited {p.returncode}: {p.stdout[-1500:]}")
    rec = os.path.join(adir, f"species.nwk_gene.ufboot.ale.uml_rec")
    if not os.path.exists(rec):
        fail(f"ALEml_undated produced no {rec}")
    D, T, L, S = parse_uml_rec(rec)
    note = (f"ALEml_undated over {len(out_trees)} UFBoot trees, {args.ale_samples} sampled "
            f"reconciliations; means; {n_pruned} genes pruned (genome not in species tree)")
    rows.extend(event_rows("ale", D, T, L, S, note))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--aln", required=True)
    ap.add_argument("--gene-tree", required=True)
    ap.add_argument("--ufboot", required=True)
    ap.add_argument("--species-tree", required=True)
    ap.add_argument("--model", default="UndatedDTL")
    ap.add_argument("--run-generax", default="True")
    ap.add_argument("--run-ale", default="False")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--strategy", default="SPR", choices=("SPR", "EVAL"),
                    help="SPR = joint gene-tree correction (GeneRax default); EVAL = reconcile the fixed tree")
    ap.add_argument("--subst-model", default="LG+G4", help="GeneRax per-family substitution model")
    ap.add_argument("--ale-samples", type=int, default=100)
    ap.add_argument("--master", default="", help="optional master table (methanogen transfer recipients)")
    ap.add_argument("--out-transfers", default="", help="optional per donor/recipient transfer table")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    rows, transfers = [], []

    def finish():
        with open(args.out, "w") as out:
            out.write(HEADER + "\n")
            for r in rows:
                out.write("\t".join(str(x) for x in r) + "\n")
        if args.out_transfers:
            with open(args.out_transfers, "w") as out:
                out.write(TRANSFER_HEADER + "\n")
                for r in transfers:
                    out.write("\t".join(str(x) for x in r) + "\n")
        utils.log(f"wrote {len(rows)} reconciliation row(s) to {args.out}")

    tools = [t for t, on in (("generax", args.run_generax), ("ale", args.run_ale)) if truthy(on)]
    if not tools:
        rows.extend(("generax", e, 0, "", "reconciliation disabled (run_generax/run_ale false)")
                    for e in EVENT_TYPES)
        return finish()
    if is_trivial(args.species_tree) or is_trivial(args.gene_tree):
        for t in tools:
            rows.extend((t, e, 0, "", "too few tips: trivial species/gene tree") for e in EVENT_TYPES)
        return finish()
    sp_path = os.path.join(args.workdir, "species_tree.labelled.nwk")
    sp_tree, sp_tips = species_tree_for_tools(args.species_tree, sp_path)
    if len(sp_tips) < MIN_TIPS:
        for t in tools:
            rows.extend((t, e, 0, "", f"too few tips: species tree has {len(sp_tips)}") for e in EVENT_TYPES)
        return finish()
    if "generax" in tools:
        run_generax(args, sp_path, sp_tree, sp_tips, rows, transfers)
    if "ale" in tools:
        run_ale(args, sp_tree, sp_tips, rows)
    finish()


if __name__ == "__main__":
    main()
