"""Per-domain alignments + ML trees, and Robinson-Foulds vs the family tree.

The four Prometheus domains (DUF4011, DUF3320, DNA2 helicase/nuclease, REase) may
not share one history (Stage 5.2: domain shuffling / fusion / independent
acquisition). For each domain we cut the sub-sequence spanning that domain's ali
coordinates out of every protein that carries it, build a domain alignment + tree,
and measure topological disagreement with the whole-protein ("family") tree via the
Robinson-Foulds distance on shared tips.

Pipeline per domain (only when >= MIN_SEQS sequences carry it):
  mafft --auto  ->  clipkit  ->  iqtree2 -m LG+G -fast

RF is computed with dendropy on a common taxon namespace restricted to tips shared
by the domain tree and the family tree; normalized RF = RF / (2*(n-3)) for n shared
tips (n>=4 needed for a non-degenerate normalization).

Never crashes the DAG: missing tools, empty inputs, <4 sequences, or trivial trees
all degrade to a logged note and (when nothing is computable) a header-only output.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils

HEADER = "domain\tn_tips\trf_to_family\tnormalized_rf\tnote"
MIN_SEQS = 4  # minimum sequences to attempt an alignment + tree for a domain


def _matches_domain(hit: dict, names: set, accs: set) -> bool:
    """True if a domtbl hit belongs to the requested domain (by name or accession)."""
    dom = str(hit.get("domain", ""))
    acc = str(hit.get("accession", ""))
    # accession match (strip PFAM version suffix, e.g. PF13195.8 -> PF13195)
    acc_base = acc.split(".")[0]
    if accs and (acc in accs or acc_base in accs):
        return True
    if names:
        dom_l = dom.lower()
        for nm in names:
            if nm and nm.lower() in dom_l:
                return True
    return False


def extract_domain_seqs(proteins: dict, hits: dict, names, accessions):
    """For one domain, return [(protein_id, subseq)] spanning the domain ali coords.

    If a protein has multiple hits for the domain, the widest span (min start ->
    max end across matching hits) is taken so a split domain is captured whole.
    """
    names = set(names or [])
    accs = set(accessions or [])
    out = []
    for pid, seq in proteins.items():
        prot_hits = hits.get(pid, [])
        spans = [
            (h["start"], h["end"]) for h in prot_hits if _matches_domain(h, names, accs)
        ]
        if not spans:
            continue
        start = max(1, min(s for s, _ in spans))
        end = min(len(seq), max(e for _, e in spans))
        if end <= start:
            continue
        sub = seq[start - 1:end]  # ali coords are 1-based, inclusive
        if sub:
            out.append((pid, sub))
    return out


def run_cmd(cmd, log_prefix):
    """Run an external command; return True on success, log stderr. Never raises."""
    utils.log(f"{log_prefix}: " + " ".join(str(c) for c in cmd))
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if proc.returncode != 0:
            err = proc.stderr.decode(errors="replace").strip()
            utils.log(f"{log_prefix}: exit {proc.returncode}: {err[:500]}")
            return False, proc.stdout
        return True, proc.stdout
    except Exception as e:  # noqa: BLE001
        utils.log(f"{log_prefix}: failed: {e}")
        return False, b""


def build_domain_tree(domain, records, workdir, threads):
    """mafft -> clipkit -> iqtree2; return path to the .treefile or None."""
    faa = os.path.join(workdir, f"{domain}.faa")
    utils.write_fasta(records, faa)
    if len(records) < MIN_SEQS:
        return None, f"{len(records)} seqs (<{MIN_SEQS})"

    mafft = shutil.which("mafft")
    iqtree = shutil.which("iqtree2") or shutil.which("iqtree")
    if mafft is None or iqtree is None:
        return None, "mafft/iqtree2 not on PATH"

    aln = os.path.join(workdir, f"{domain}.aln")
    ok, stdout = run_cmd([mafft, "--auto", "--thread", str(threads), faa], f"[{domain}] mafft")
    if not ok:
        return None, "mafft failed"
    with open(aln, "wb") as fh:
        fh.write(stdout)

    # ClipKIT trims; if it fails, fall back to the untrimmed alignment.
    trim = os.path.join(workdir, f"{domain}.trim.aln")
    clipkit = shutil.which("clipkit")
    aln_for_tree = aln
    if clipkit is not None:
        ok, _ = run_cmd([clipkit, aln, "-o", trim], f"[{domain}] clipkit")
        if ok and os.path.exists(trim) and os.path.getsize(trim) > 0:
            aln_for_tree = trim
        else:
            utils.log(f"[{domain}] clipkit unavailable result; using untrimmed alignment")
    else:
        utils.log(f"[{domain}] clipkit not on PATH; using untrimmed alignment")

    prefix = os.path.join(workdir, domain)
    ok, _ = run_cmd(
        [iqtree, "-s", aln_for_tree, "-m", "LG+G", "-fast",
         "-T", str(threads), "--prefix", prefix, "-redo"],
        f"[{domain}] iqtree2",
    )
    treefile = prefix + ".treefile"
    if ok and os.path.exists(treefile) and os.path.getsize(treefile) > 0:
        return treefile, "ok"
    return None, "iqtree2 failed"


def is_trivial_tree(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return True
    try:
        with open(path) as fh:
            txt = fh.read().strip()
    except OSError:
        return True
    return txt in ("", "();", ";")


def compute_rf(domain_tree_path, family_tree_path):
    """RF + normalized RF on shared tips. Returns (n_shared, rf, norm_rf, note)."""
    try:
        import dendropy
        from dendropy.calculate import treecompare
    except Exception as e:  # noqa: BLE001
        return 0, "", "", f"dendropy unavailable: {e}"

    if is_trivial_tree(domain_tree_path) or is_trivial_tree(family_tree_path):
        return 0, "", "", "trivial tree"

    tns = dendropy.TaxonNamespace()
    try:
        dt = dendropy.Tree.get(path=domain_tree_path, schema="newick",
                               taxon_namespace=tns, preserve_underscores=True)
        ft = dendropy.Tree.get(path=family_tree_path, schema="newick",
                               taxon_namespace=tns, preserve_underscores=True)
    except Exception as e:  # noqa: BLE001
        return 0, "", "", f"tree parse error: {e}"

    dt_tips = {lf.taxon.label for lf in dt.leaf_node_iter() if lf.taxon}
    ft_tips = {lf.taxon.label for lf in ft.leaf_node_iter() if lf.taxon}
    shared = dt_tips & ft_tips
    n = len(shared)
    if n < 4:
        return n, "", "", f"only {n} shared tips (<4)"

    # Restrict both trees to the shared taxon set, then compare on a common namespace.
    try:
        dt.retain_taxa_with_labels(shared)
        ft.retain_taxa_with_labels(shared)
        dt.encode_bipartitions()
        ft.encode_bipartitions()
        rf = treecompare.symmetric_difference(dt, ft)
    except Exception as e:  # noqa: BLE001
        return n, "", "", f"RF computation error: {e}"

    denom = 2 * (n - 3)
    norm = (rf / denom) if denom > 0 else ""
    return n, rf, (round(norm, 6) if norm != "" else ""), "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--proteins", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--family-tree", required=True)
    ap.add_argument("--domains", required=True, help="JSON {Domain:{names:[],accessions:[]}}")
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    try:
        domains = json.loads(args.domains) if args.domains else {}
    except json.JSONDecodeError as e:
        utils.log(f"could not parse --domains JSON: {e}")
        domains = {}

    proteins = dict(utils.read_fasta(args.proteins))
    hits = utils.parse_domtbl(args.domtbl)

    if not domains or not proteins:
        utils.write_empty(args.out, HEADER)
        utils.log("no domains or no proteins; header-only output")
        return

    rows = []
    for domain, spec in domains.items():
        spec = spec or {}
        records = extract_domain_seqs(
            proteins, hits, spec.get("names", []), spec.get("accessions", [])
        )
        n_seqs = len(records)
        if n_seqs == 0:
            rows.append((domain, 0, "", "", "no proteins carry this domain"))
            continue

        treefile, build_note = build_domain_tree(domain, records, args.workdir, args.threads)
        if treefile is None:
            rows.append((domain, n_seqs, "", "", f"no tree: {build_note}"))
            continue

        n_shared, rf, norm, rf_note = compute_rf(treefile, args.family_tree)
        rows.append((domain, n_shared if n_shared else n_seqs, rf, norm, rf_note))

    # Write the table (always at least a header).
    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for domain, n_tips, rf, norm, note in rows:
            out.write(f"{domain}\t{n_tips}\t{rf}\t{norm}\t{note}\n")

    utils.log(f"wrote {len(rows)} domain row(s) to {args.out}")


if __name__ == "__main__":
    main()
