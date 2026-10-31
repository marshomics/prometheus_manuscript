"""Synthetic inputs and stand-in tools for the preprocess_qc test suite.

* planted_genome(): a deterministic genome whose genes are known exactly. Genes are
  random sense codons (no internal stop and no ATG/GTG/TTG, so the planted ATG is the only
  possible start) behind an AGGAGG RBS; intergenic DNA is stop codons in all six frames.
  Real Prodigal 2.6.3 (-p single -g 11) calls exactly the planted genes, with exact
  coordinates, strands and partial codes (verified while writing this suite).
* write_fake_tools(): stand-ins for `gtdbtk`, `checkm2` and `statswrapper.sh` that emit
  the file layout / columns of the pinned real versions (verified from GTDB-Tk 2.4.0 and
  CheckM2 1.0.2 source and from a real BBMap 39.06 run, see REVIEW.md), so the
  Snakemake rules' shell code is exercised end to end without the 100-GB databases.
"""
from __future__ import annotations
import os
import random
import stat

STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"
          if a + b + c not in STOPS and a + b + c not in STARTS]
RBS = "AGGAGG"
_B = "TCAG"
_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11 = {a + b + c: _AA[i] for i, (a, b, c) in
          enumerate((a, b, c) for a in _B for b in _B for c in _B)}


def rc(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def translate(nt: str) -> str:
    """Table-11 translation as Prodigal prints it (stop -> '*')."""
    return "".join(CODE11.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def spacer(n: int) -> str:
    unit = "TAATTAGTTAACTAGT"   # stop codons in all three frames on both strands
    return (unit * (n // len(unit) + 1))[:n]


def planted_genome(seed=7, layout=((14, False), (0, False), (6, False), (1, True))):
    """layout = per contig (n_genes, last_gene_runs_off_right_end).

    Returns (sequences, genes) with genes = [(contig_index, start, end, strand, partial,
    n_codons_without_stop)], 1-based inclusive coordinates; every third gene of a contig
    is on the minus strand."""
    rnd = random.Random(seed)
    records, genes = [], []
    for ci, (ng, edge) in enumerate(layout, start=1):
        seq = spacer(30)
        for gi in range(ng):
            ncod = rnd.randint(250, 400)
            body = "ATG" + "".join(rnd.choice(CODONS) for _ in range(ncod - 1))
            if edge and gi == ng - 1:
                block = spacer(40) + RBS + "ATAACT" + body       # no stop codon
                s = len(seq) + len(block) - len(body) + 1
                seq += block
                genes.append((ci, s, len(seq), "+", "01", ncod))
                break
            g = body + "TAA"
            block = spacer(40) + RBS + "ATAACT" + g
            if gi % 3 == 2:
                start = len(seq) + 1
                seq += rc(block)
                genes.append((ci, start, start + len(g) - 1, "-", "00", ncod))
            else:
                s = len(seq) + len(block) - len(g) + 1
                seq += block
                genes.append((ci, s, len(seq), "+", "00", ncod))
        else:
            seq += spacer(60) if ng else spacer(600)
        records.append(seq)
    return records, genes


def write_fasta(path, names, seqs, width=70):
    with open(path, "w") as fh:
        for n, s in zip(names, seqs):
            fh.write(f">{n}\n")
            for i in range(0, len(s), width):
                fh.write(s[i:i + width] + "\n")


# --------------------------------------------------------------------------- #
# fake tools (formats of GTDB-Tk 2.4.0 / CheckM2 1.0.2 / BBMap 39.06)
# --------------------------------------------------------------------------- #
GTDBTK_HEADER = ["user_genome", "classification", "closest_genome_reference",
                 "closest_genome_reference_radius", "closest_genome_taxonomy", "closest_genome_ani",
                 "closest_genome_af", "closest_placement_reference", "closest_placement_radius",
                 "closest_placement_taxonomy", "closest_placement_ani", "closest_placement_af",
                 "pplacer_taxonomy", "classification_method", "note",
                 "other_related_references(genome_id,species_name,radius,ANI,AF)", "msa_percent",
                 "translation_table", "red_value", "warnings"]
CHECKM2_HEADER = ["Name", "Completeness", "Contamination", "Completeness_Model_Used",
                  "Translation_Table_Used", "Coding_Density", "Contig_N50", "Average_Gene_Length",
                  "Genome_Size", "GC_Content", "Total_Coding_Sequences", "Total_Contigs",
                  "Max_Contig_Length", "Additional_Notes"]
BBMAP_HEADER = ["n_scaffolds", "n_contigs", "scaf_bp", "contig_bp", "gap_pct", "scaf_N50", "scaf_L50",
                "ctg_N50", "ctg_L50", "scaf_N90", "scaf_L90", "ctg_N90", "ctg_L90", "scaf_max",
                "ctg_max", "scaf_n_gt50K", "scaf_pct_gt50K", "gc_avg", "gc_std", "filename"]

_COMMON = r'''
import os, sys, re
def read_fasta(p):
    recs, name, seq = [], None, []
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            if name is not None: recs.append((name, "".join(seq)))
            name, seq = line[1:], []
        elif line: seq.append(line)
    if name is not None: recs.append((name, "".join(seq)))
    return recs
def nl50(lengths, frac=0.5):
    """BBMap naming: returns (N = count, L = length) at the given fraction."""
    ls = sorted(lengths, reverse=True); tot = sum(ls); acc = 0
    for i, l in enumerate(ls, 1):
        acc += l
        if acc >= tot * frac: return i, l
    return 0, 0
'''

FAKE_GTDBTK = _COMMON + r'''
import argparse
a = sys.argv[1:]
if not a or a[0] != "classify_wf":
    sys.exit("fake gtdbtk: only classify_wf is supported")
ap = argparse.ArgumentParser()
ap.add_argument("--batchfile", required=True); ap.add_argument("--out_dir", required=True)
ap.add_argument("--cpus"); ap.add_argument("--skip_ani_screen", action="store_true")
o = ap.parse_args(a[1:])
if os.environ.get("FAKE_GTDBTK_FAIL"):
    sys.exit("fake GTDB-Tk: simulated crash")
tax = {}
for line in open(os.environ["FAKE_GTDBTK_TAX"]):
    g, c = line.rstrip("\n").split("\t")
    tax[g] = c
drop = set(os.environ.get("FAKE_GTDBTK_DROP", "").split(","))
rows = {"ar53": [], "bac120": []}
for line in open(o.batchfile):
    path, gid = line.rstrip("\n").split("\t")[:2]
    if not os.path.isfile(path): sys.exit(f"fake gtdbtk: missing {path}")
    if not gid.startswith("user_"): sys.exit(f"fake gtdbtk: GTDB-Tk id without prefix: {gid}")
    g = gid[len("user_"):]
    if g in drop: continue
    c = tax.get(g, "Unclassified")
    ms = "bac120" if (c.startswith("d__Bacteria") or c == "Unclassified") else "ar53"
    rows[ms].append((gid, c))
HDR = __HDR__
os.makedirs(os.path.join(o.out_dir, "classify"), exist_ok=True)
os.makedirs(os.path.join(o.out_dir, "identify"), exist_ok=True)
for ms, rr in rows.items():
    if not rr: continue
    rel = f"classify/gtdbtk.{ms}.summary.tsv"
    with open(os.path.join(o.out_dir, rel), "w") as fh:
        fh.write("\t".join(HDR) + "\n")
        for gid, c in rr:
            vals = [gid, c] + ["N/A"] * (len(HDR) - 2)
            vals[13] = "taxonomic classification defined by topology and ANI"
            fh.write("\t".join(vals) + "\n")
    # GTDB-Tk 2.4.0 also links each summary from the out_dir root and writes a
    # (different) markers summary in identify/; both match '*summary.tsv'.
    os.symlink(rel, os.path.join(o.out_dir, f"gtdbtk.{ms}.summary.tsv"))
    with open(os.path.join(o.out_dir, f"identify/gtdbtk.{ms}.markers_summary.tsv"), "w") as fh:
        fh.write("name\tnumber_unique_genes\tnumber_multiple_genes\n")
        for gid, _ in rr: fh.write(f"{gid}\t50\t0\n")
print("fake gtdbtk done")
'''.replace("__HDR__", repr(GTDBTK_HEADER))

FAKE_CHECKM2 = _COMMON + r'''
import argparse
if "--version" in sys.argv:
    print("1.0.2"); sys.exit(0)
a = sys.argv[1:]
if not a or a[0] != "predict": sys.exit("fake checkm2: only predict is supported")
ap = argparse.ArgumentParser()
ap.add_argument("--threads"); ap.add_argument("--force", action="store_true")
ap.add_argument("--database_path", required=True); ap.add_argument("--input", required=True)
ap.add_argument("-x", "--extension", default=".fna"); ap.add_argument("--output-directory", required=True)
ap.add_argument("--remove_intermediates", action="store_true")
o = ap.parse_args(a[1:])
if os.environ.get("FAKE_CHECKM2_FAIL"): sys.exit("fake CheckM2: simulated crash")
if not os.path.isfile(o.database_path): sys.exit("fake CheckM2: no database")
q = {}
for line in open(os.environ["FAKE_CHECKM2_TABLE"]):
    g, comp, cont = line.split()
    q[g] = (comp, cont)
skip = set(os.environ.get("FAKE_CHECKM2_SKIP", "").split(","))
os.makedirs(o.output_directory, exist_ok=True)
HDR = __HDR__
with open(os.path.join(o.output_directory, "quality_report.tsv"), "w") as fh:
    fh.write("\t".join(HDR) + "\n")
    for f in sorted(os.listdir(o.input)):
        if not f.endswith(o.extension): continue
        name = os.path.splitext(os.path.basename(f))[0]      # CheckM2 1.0.2 bin naming
        if name in skip: continue
        recs = read_fasta(os.path.join(o.input, f))
        lens = [len(s) for _, s in recs]
        comp, cont = q.get(name, ("50.00", "50.00"))
        n50 = nl50(lens)[1]
        fh.write("\t".join([name, comp, cont, "Neural Network (Specific Model)", "11", "0.9",
                            str(n50), "300", str(sum(lens)), "0.5", "10", str(len(lens)),
                            str(max(lens)), "None"]) + "\n")
'''.replace("__HDR__", repr(CHECKM2_HEADER))

FAKE_STATSWRAPPER = _COMMON + r'''
files = []
for a in sys.argv[1:]:
    if a.startswith("in="): files += a[3:].split(",")
if os.environ.get("FAKE_BBMAP_FAIL"): sys.exit(1)
omit = set(os.environ.get("FAKE_BBMAP_OMIT", "").split(","))
print("#" + "\t".join(__HDR__))
for f in files:
    if not os.path.exists(f):
        print(f"Exception in thread main: Input file does not appear to be valid: {f}", file=sys.stderr)
        sys.exit(1)
    if os.path.basename(f).rsplit(".", 1)[0] in omit: continue
    recs = read_fasta(f)
    scaf = [len(s) for _, s in recs]
    ctg = [len(p) for _, s in recs for p in re.split(r"N{10,}", s.upper()) if p]
    gc = sum(s.upper().count("G") + s.upper().count("C") for _, s in recs)
    acgt = sum(sum(s.upper().count(b) for b in "ACGT") for _, s in recs)
    sN, sL = nl50(scaf); cN, cL = nl50(ctg); sN9, sL9 = nl50(scaf, 0.9); cN9, cL9 = nl50(ctg, 0.9)
    gap = "%.3f" % ((sum(scaf) - sum(ctg)) * 100.0 / sum(scaf))
    row = [len(scaf), len(ctg), sum(scaf), sum(ctg), gap, sN, sL, cN, cL, sN9, sL9, cN9, cL9,
           max(scaf), max(ctg), 0, "0.000", "%.5f" % (gc / acgt), "0.00000", f]
    print("\t".join(str(x) for x in row))
'''.replace("__HDR__", repr(BBMAP_HEADER))


def write_fake_tools(bindir, python_exe):
    """Write executable fakes gtdbtk / checkm2 / statswrapper.sh into bindir."""
    os.makedirs(bindir, exist_ok=True)
    for name, body in (("gtdbtk", FAKE_GTDBTK), ("checkm2", FAKE_CHECKM2),
                       ("statswrapper.sh", FAKE_STATSWRAPPER)):
        p = os.path.join(bindir, name)
        with open(p, "w") as fh:
            fh.write(f"#!{python_exe}\n{body}")
        os.chmod(p, os.stat(p).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return bindir


def write_fake_gtdbtk_module(pydir, version="2.4.0"):
    """`python -c 'import gtdbtk'` in the gtdbtk_batch rule records the version."""
    d = os.path.join(pydir, "gtdbtk")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "__init__.py"), "w") as fh:
        fh.write(f"__version__ = {version!r}\n")
    return pydir
