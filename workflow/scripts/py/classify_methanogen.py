"""Classify each genome as methanogen vs non-methanogen, and separate true
methanogens from anaerobic methanotrophs (ANME) by McrA phylogenetic placement
(design §3.3).

Logic
-----
methanogen      : McrA present (McrA is the diagnostic minimum for the Mcr
                  complex; design §3.3 / config methanogen.mcr_required).
n_mcr           : number of distinct mcr subunits (mcrA/mcrB/mcrG) detected.
is_anme         : ANME carry a *divergent* McrA that they run in reverse; they
                  must NOT be silently counted as methanogens. We test this by
                  adding each candidate McrA to a curated reference McrA alignment
                  and finding, on a fast tree, the nearest reference tip; if that
                  reference is labelled 'ANME' the genome is a methanotroph
                  (methanogen set to False, is_anme True).

methanogenesis_subtype: reported as 'unresolved' for methanogens and 'NA' for
non-methanogens. We deliberately do NOT call hydrogenotrophic / methylotrophic /
acetoclastic / methyl-reducing here: subtype assignment requires pathway-marker
analysis (Wood-Ljungdahl/methyl-branch markers, mta/mtb/mtt methyltransferases,
hydrogenases — design §3.3), which is out of scope for an mcr-only step. Calling
a subtype from McrA alone would overclaim.

Failure handling: every tool in the ANME check (mafft --add, iqtree2/FastTree,
dendropy) is optional. On any failure we set is_anme=False with
evidence='anme_check_skipped' rather than dropping or mis-labelling the genome.
dendropy is imported lazily so this module imports even where it is absent.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

import argparse
import shutil
import subprocess

HEADER = (
    "genome_id\tmcrA\tn_mcr\tmethanogen\tis_anme\t"
    "methanogenesis_subtype\tevidence"
)
MCR_GENES = ("mcrA", "mcrB", "mcrG")


def load_ref_labels(path: str) -> dict:
    """ref McrA id -> group ('ANME' | 'methanogen' | 'other'); TSV id<TAB>group."""
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
            rid, grp = parts[0].strip(), parts[1].strip()
            if rid.lower() == "ref_id" or grp.lower() == "group":
                continue  # header
            labels[rid] = grp
    return labels


def read_mcr_hits(path: str):
    """Return (per_genome_genes, mcrA_protein_ids).

    per_genome_genes : genome -> set of mcr genes detected
    mcrA_protein_ids : genome -> protein_id of its McrA hit (for sequence extraction)
    Tolerates a header-only / missing file (returns empties).
    """
    per_genome = {}
    mcrA_pid = {}
    if not path or not os.path.exists(path):
        return per_genome, mcrA_pid
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            gi = header.index("genome_id")
            ge = header.index("gene")
            pi = header.index("protein_id")
        except ValueError:
            return per_genome, mcrA_pid
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            f = line.split("\t")
            if len(f) <= max(gi, ge, pi):
                continue
            genome, gene, pid = f[gi].strip(), f[ge].strip(), f[pi].strip()
            if gene not in MCR_GENES:
                continue
            per_genome.setdefault(genome, set()).add(gene)
            if gene == "mcrA":
                mcrA_pid[genome] = pid
    return per_genome, mcrA_pid


def build_mcrA_tree(cand_pids, proteome_path, ref_aln, threads, workdir):
    """Add candidate McrA to the reference alignment and infer a fast tree.

    Returns the path to the resulting newick treefile, or None on any failure.
    Steps: extract candidate McrA seqs -> `mafft --add ... --keeplength ref_aln`
    -> `iqtree2 -fast` (or FastTree fallback).
    """
    os.makedirs(workdir, exist_ok=True)
    if shutil.which("mafft") is None:
        utils.log("mafft not on PATH; skipping ANME tree.")
        return None
    if not ref_aln or not os.path.exists(ref_aln) or os.path.getsize(ref_aln) == 0:
        utils.log("Reference McrA alignment missing/empty; skipping ANME tree.")
        return None

    # extract candidate McrA sequences from the proteome
    seqs = dict(utils.read_fasta(proteome_path))
    cand_records = [(pid, seqs[pid]) for pid in cand_pids if pid in seqs and seqs[pid]]
    if not cand_records:
        utils.log("No candidate McrA sequences found in proteome; skipping ANME tree.")
        return None
    cand_faa = os.path.join(workdir, "candidates.faa")
    utils.write_fasta(cand_records, cand_faa)

    combined = os.path.join(workdir, "combined.aln")
    # --keeplength preserves reference columns so candidates are placed in the
    # same coordinate frame; references already carry the ANME/methanogen labels.
    try:
        with open(combined, "w") as out:
            res = subprocess.run(
                ["mafft", "--add", cand_faa, "--keeplength", ref_aln],
                stdout=out, stderr=subprocess.PIPE, text=True)
        if res.returncode != 0:
            utils.log(f"WARNING: mafft --add returncode={res.returncode}; "
                      f"stderr tail: {res.stderr[-300:]}")
            return None
    except Exception as exc:
        utils.log(f"WARNING: mafft --add failed: {exc}")
        return None
    if not os.path.exists(combined) or os.path.getsize(combined) == 0:
        return None

    # Fast tree: prefer iqtree2, fall back to FastTree/fasttree.
    if shutil.which("iqtree2") is not None:
        prefix = os.path.join(workdir, "mcrA")
        try:
            res = subprocess.run(
                ["iqtree2", "-s", combined, "-m", "LG+G", "-fast",
                 "-T", str(threads), "--prefix", prefix, "-redo"],
                capture_output=True, text=True)
        except Exception as exc:
            utils.log(f"WARNING: iqtree2 failed to launch: {exc}")
            res = None
        treefile = prefix + ".treefile"
        if res is not None and os.path.exists(treefile) and os.path.getsize(treefile) > 0:
            return treefile
        utils.log("iqtree2 produced no tree; trying FastTree.")

    fasttree_bin = shutil.which("FastTree") or shutil.which("fasttree")
    if fasttree_bin is not None:
        treefile = os.path.join(workdir, "mcrA_fasttree.nwk")
        try:
            with open(treefile, "w") as out:
                res = subprocess.run(
                    [fasttree_bin, "-lg", combined],
                    stdout=out, stderr=subprocess.PIPE, text=True)
            if res.returncode == 0 and os.path.getsize(treefile) > 0:
                return treefile
        except Exception as exc:
            utils.log(f"WARNING: FastTree failed: {exc}")
    utils.log("No tree builder succeeded; skipping ANME placement.")
    return None


def anme_calls_from_tree(treefile, ref_labels, cand_pids):
    """For each candidate McrA tip, is its nearest reference tip an ANME?

    Returns dict {candidate_protein_id -> bool is_anme}. Empty on failure (then
    callers treat every candidate as is_anme=False / anme_check_skipped).
    """
    out = {}
    if not treefile or not os.path.exists(treefile):
        return out
    try:
        import dendropy
    except Exception:
        utils.log("dendropy unavailable; cannot place McrA -> ANME check skipped.")
        return out
    try:
        tree = dendropy.Tree.get(path=treefile, schema="newick",
                                 preserve_underscores=True)
        pdm = tree.phylogenetic_distance_matrix()
    except Exception as exc:
        utils.log(f"WARNING: could not analyse McrA tree ({exc}); ANME skipped.")
        return out

    taxa = [t for t in tree.taxon_namespace if t.label is not None]
    taxon_by_label = {t.label: t for t in taxa}
    ref_in_tree = {lbl: ref_labels[lbl] for lbl in taxon_by_label if lbl in ref_labels}
    if not ref_in_tree:
        utils.log("No reference McrA tips in tree; ANME check skipped.")
        return out

    cand_set = set(cand_pids)
    for lbl, taxon in taxon_by_label.items():
        if lbl not in cand_set:
            continue
        best_group, best_dist = None, None
        for ref_lbl, group in ref_in_tree.items():
            try:
                d = pdm.patristic_distance(taxon, taxon_by_label[ref_lbl])
            except Exception:
                continue
            if d is None:
                continue
            if best_dist is None or d < best_dist:
                best_dist, best_group = d, group
        if best_group is not None:
            out[lbl] = (str(best_group).strip().upper() == "ANME")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hits", required=True, help="mcr_hits.tsv")
    ap.add_argument("--proteome", required=True)
    ap.add_argument("--ref-aln", required=True, help="reference McrA alignment")
    ap.add_argument("--ref-labels", required=True, help="ref McrA id -> group TSV")
    ap.add_argument("--required", required=True,
                    help="comma KO/gene list (e.g. K00399 == mcrA)")
    ap.add_argument("--anme-tree", required=True,
                    help="'True'/'False': run McrA ANME placement")
    ap.add_argument("--threads", default="1")
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    os.makedirs(args.workdir, exist_ok=True)

    run_anme = str(args.anme_tree).strip().lower() in ("true", "1", "yes")

    per_genome, mcrA_pid = read_mcr_hits(args.hits)
    genomes = sorted(per_genome.keys())

    if not genomes:
        # Accept an empty/header-only hits file: nothing to classify.
        utils.write_empty(args.out, HEADER)
        utils.log("No mcr hits; wrote header-only methanogen_calls.")
        return

    # ANME placement (once, for all candidate McrA together)
    anme_by_pid = {}
    anme_attempted = False
    if run_anme:
        cand_pids = [mcrA_pid[g] for g in genomes if g in mcrA_pid]
        ref_labels = load_ref_labels(args.ref_labels)
        if cand_pids and ref_labels:
            anme_attempted = True
            treefile = build_mcrA_tree(
                cand_pids, args.proteome, args.ref_aln, args.threads, args.workdir)
            anme_by_pid = anme_calls_from_tree(treefile, ref_labels, cand_pids)

    rows = []
    for g in genomes:
        genes = per_genome.get(g, set())
        has_mcrA = "mcrA" in genes
        n_mcr = len([x for x in MCR_GENES if x in genes])
        methanogen = bool(has_mcrA)
        is_anme = False
        evidence_bits = []

        if has_mcrA:
            evidence_bits.append("mcrA_present")
        elif n_mcr >= 1:
            # mcr subunit(s) detected but the diagnostic McrA is absent: not a
            # methanogen call (could be a partial/fragmented complex).
            present = ",".join(x for x in MCR_GENES if x in genes)
            evidence_bits.append(f"mcrA_absent;mcr_subunits={present}")
        if has_mcrA and n_mcr >= 2:
            evidence_bits.append(f"n_mcr={n_mcr}")

        if run_anme and has_mcrA:
            pid = mcrA_pid.get(g)
            if pid is not None and pid in anme_by_pid:
                is_anme = bool(anme_by_pid[pid])
                if is_anme:
                    methanogen = False  # methanotroph runs Mcr in reverse
                    evidence_bits.append("mcrA_nearest_ref=ANME")
                else:
                    evidence_bits.append("mcrA_nearest_ref=methanogen_side")
            else:
                # tool chain did not yield a placement for this genome
                evidence_bits.append("anme_check_skipped")
        elif run_anme and not anme_attempted:
            # ANME requested but prerequisites absent
            pass

        subtype = "unresolved" if methanogen else "NA"
        evidence = ";".join(evidence_bits) if evidence_bits else "no_mcr"

        rows.append({
            "genome_id": g,
            "mcrA": has_mcrA,
            "n_mcr": n_mcr,
            "methanogen": methanogen,
            "is_anme": is_anme,
            "methanogenesis_subtype": subtype,
            "evidence": evidence,
        })

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        cols = HEADER.split("\t")
        for r in rows:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_meth = sum(1 for r in rows if r["methanogen"])
    n_anme = sum(1 for r in rows if r["is_anme"])
    utils.log(f"methanogen_calls: {len(rows)} genomes; {n_meth} methanogen, "
              f"{n_anme} ANME.")


if __name__ == "__main__":
    main()
