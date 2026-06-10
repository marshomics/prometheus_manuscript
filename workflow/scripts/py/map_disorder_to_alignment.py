"""Map per-residue consensus disorder onto the Prometheus family alignment.

Comparing whole-protein %disorder across divergent sequences is confounded by length
and domain content. Anchoring disorder to the alignment lets the comparison be made on
HOMOLOGOUS positions: per alignment column we report the mean consensus-disorder and
the fraction of loci disordered there, and we tag each column by domain (vs linker)
using the reference (longest-ungapped) locus' domain hits. This is what feeds the
disorder-along-protein figure and the linker-vs-domain contrast.

Empty/missing inputs -> header-only output.
"""
from __future__ import annotations
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

HEADER = "alignment_col\tn_loci\tmean_disorder\tfrac_loci_disordered\tdomain"


def best_domain_for_hit(hit_name, hit_acc, domains):
    """Single most-specific domain category for an hmm hit (mirror validate_prometheus)."""
    name_l = str(hit_name).lower()
    acc_base = str(hit_acc).split(".")[0]
    for dname, spec in domains.items():
        for a in (spec or {}).get("accessions", []):
            if a and (a == hit_acc or a.split(".")[0] == acc_base):
                return dname
    order = list(domains.keys())
    best_key = best = None
    for dname, spec in domains.items():
        for nm in (spec or {}).get("names", []):
            nl = str(nm).lower()
            if nl and nl in name_l:
                key = (len(nl), -order.index(dname))
                if best_key is None or key > best_key:
                    best_key, best = key, dname
    return best


def reference_residue_domains(ref_id, hits, domains):
    """residue(1-based) -> domain for the reference protein, from its domtbl hits."""
    res_dom = {}
    for h in hits.get(ref_id, []):
        dom = best_domain_for_hit(h.get("domain", ""), h.get("accession", ""), domains)
        if dom:
            for r in range(int(h["start"]), int(h["end"]) + 1):
                res_dom.setdefault(r, dom)
    return res_dom


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--residue", required=True, help="disorder_per_residue.tsv")
    ap.add_argument("--alignment", required=True, help="family.aln (protein MSA)")
    ap.add_argument("--domtbl", default="", help="hmmscan domtbl for domain track")
    ap.add_argument("--domains", default="{}", help="domain JSON")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    aln = list(utils.read_fasta(args.alignment))
    if not aln or not os.path.exists(args.residue) or os.path.getsize(args.residue) == 0:
        utils.write_empty(args.out, HEADER)
        utils.log("empty alignment or residue table; header-only.")
        return

    # consensus disorder per (locus, residue position)
    cons = {}
    with open(args.residue) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            li, pi, ci = header.index("locus_id"), header.index("position"), header.index("consensus")
        except ValueError:
            utils.write_empty(args.out, HEADER); return
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(li, pi, ci):
                continue
            try:
                cons[(f[li], int(f[pi]))] = int(float(f[ci]))
            except ValueError:
                continue
    if not cons:
        utils.write_empty(args.out, HEADER)
        utils.log("no per-residue disorder values; header-only.")
        return

    try:
        domains = json.loads(args.domains) if args.domains else {}
    except json.JSONDecodeError:
        domains = {}
    hits = utils.parse_domtbl(args.domtbl) if args.domtbl else {}

    aln_len = max(len(s) for _, s in aln)
    # reference = the alignment record with the most ungapped residues that is in cons
    ref_id, ref_seq = max(aln, key=lambda kv: sum(c not in "-." for c in kv[1]))
    ref_res_dom = reference_residue_domains(ref_id, hits, domains) if hits else {}

    # walk columns; per locus track its ungapped residue index
    col_sum = [0.0] * aln_len
    col_n = [0] * aln_len
    for locus, seq in aln:
        ridx = 0
        for col in range(len(seq)):
            c = seq[col]
            if c in "-.":
                continue
            ridx += 1
            v = cons.get((locus, ridx))
            if v is None:
                continue
            col_sum[col] += v
            col_n[col] += 1
    # reference residue index per column (for the domain track)
    ref_dom_by_col = {}
    ridx = 0
    for col in range(len(ref_seq)):
        if ref_seq[col] in "-.":
            continue
        ridx += 1
        d = ref_res_dom.get(ridx)
        if d:
            ref_dom_by_col[col] = d

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for col in range(aln_len):
            if col_n[col] == 0:
                continue
            mean_dis = col_sum[col] / col_n[col]
            out.write(f"{col+1}\t{col_n[col]}\t{round(mean_dis,4)}\t"
                      f"{round(mean_dis,4)}\t{ref_dom_by_col.get(col,'NA')}\n")
    utils.log(f"disorder-by-column: {sum(1 for n in col_n if n)} columns over "
              f"{len(aln)} loci (reference {ref_id}).")


if __name__ == "__main__":
    main()
