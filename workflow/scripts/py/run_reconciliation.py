"""Species-tree-aware DTL reconciliation of the Prometheus gene tree (GeneRax; opt. ALE).

Reconciles the Prometheus gene tree against the archaeal species tree (Stage 5.3) to
infer Duplication / Transfer / Loss / Speciation events and a vertical-vs-transfer
fraction. GeneRax gives the ML point estimate under a DTL model; ALE (optional, often
absent) re-estimates from the bootstrap distribution to fold in gene-tree uncertainty.

Inputs are built defensively:
  - gene tip -> species mapping from utils.genome_of, written as a GeneRax per-family
    mapping file ("geneName speciesName", one pair per line);
  - a GeneRax families file pointing at the gene tree, alignment, mapping, and model.

Guards (never crash the DAG):
  - if run-generax != 'True', or the species/gene tree is trivial ("();" / <4 tips),
    a header-only events table is written and the function returns;
  - every external call (generax, ALEobserve, ALEml_undated) is wrapped; a tool that
    is missing or fails yields a count-0 row with an explanatory note, so the table
    always carries one row per tool attempted.
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
MIN_TIPS = 4


def tree_text(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return ""
    try:
        with open(path) as fh:
            return fh.read().strip()
    except OSError:
        return ""


def is_trivial(path):
    txt = tree_text(path)
    return txt in ("", "();", ";")


def count_tips(path):
    """Cheap tip count from a newick string (count leaf labels)."""
    txt = tree_text(path)
    if not txt or txt in ("();", ";"):
        return 0
    # leaves are tokens not immediately preceded by ')'; count via dendropy if available
    try:
        import dendropy
        t = dendropy.Tree.get(data=txt, schema="newick", preserve_underscores=True)
        return sum(1 for _ in t.leaf_node_iter())
    except Exception:  # noqa: BLE001 - fall back to a regex tip count
        # labels follow '(' or ',' and precede ':' ',' ')'
        return len(re.findall(r"[\(,]\s*([A-Za-z0-9_][^,:\(\)]*)", txt))


def write_header_only(path, note, tools=("generax",)):
    """Header + one count-0 Speciation row per tool, so downstream code has a shape."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as out:
        out.write(HEADER + "\n")
        for tool in tools:
            out.write(f"{tool}\tSpeciation\t0\t\t{note}\n")


def build_mapping(gene_tree_path, map_path):
    """Write a GeneRax per-family mapping: 'geneName speciesName' per line.

    species = utils.genome_of(geneTip). Returns the number of mapped tips.
    """
    try:
        import dendropy
        t = dendropy.Tree.get(path=gene_tree_path, schema="newick", preserve_underscores=True)
        tips = [lf.taxon.label for lf in t.leaf_node_iter() if lf.taxon]
    except Exception as e:  # noqa: BLE001
        utils.log(f"could not read gene-tree tips with dendropy: {e}; using regex")
        txt = tree_text(gene_tree_path)
        tips = re.findall(r"[\(,]\s*([A-Za-z0-9_][^,:\(\)]*)", txt)
    os.makedirs(os.path.dirname(os.path.abspath(map_path)), exist_ok=True)
    n = 0
    with open(map_path, "w") as out:
        for tip in tips:
            tip = tip.strip()
            if not tip:
                continue
            species = utils.genome_of(tip)
            out.write(f"{tip} {species}\n")
            n += 1
    return n


def write_families(families_path, gene_tree, aln, map_path, subst_model="LG+G"):
    os.makedirs(os.path.dirname(os.path.abspath(families_path)), exist_ok=True)
    with open(families_path, "w") as out:
        out.write("[FAMILIES]\n")
        out.write("- prometheus\n")
        out.write(f"starting_gene_tree = {os.path.abspath(gene_tree)}\n")
        out.write(f"alignment = {os.path.abspath(aln)}\n")
        out.write(f"mapping = {os.path.abspath(map_path)}\n")
        out.write(f"subst_model = {subst_model}\n")


# ---- GeneRax output parsing ------------------------------------------------

# event tags appearing in GeneRax NHX reconciliations
_NHX_EVENT = {
    "D=1": "Duplication", "D=Y": "Duplication",
    "H=Y": "Transfer", "T=Y": "Transfer",
}


def find_generax_outputs(rundir):
    """Return (events_files, nhx_files) discovered anywhere under the GeneRax run dir."""
    events_files, nhx_files = [], []
    if not os.path.isdir(rundir):
        return events_files, nhx_files
    for root, _dirs, files in os.walk(rundir):
        for f in files:
            p = os.path.join(root, f)
            low = f.lower()
            if "event" in low and low.endswith(".txt"):
                events_files.append(p)
            elif low.endswith((".nhx", ".newick", ".nwk")) and "reconcil" in root.lower():
                nhx_files.append(p)
            elif low.endswith(".nhx"):
                nhx_files.append(p)
    return events_files, nhx_files


def parse_events_txt(paths):
    """Sum event counts from GeneRax *_events.txt files.

    GeneRax writes lines like 'S = 12', 'D = 3', 'T = 2', 'L = 5' (or full words).
    Returns dict {EventType: count} for any recognised tokens.
    """
    counts = {e: 0 for e in EVENT_TYPES}
    found_any = False
    key_map = {
        "S": "Speciation", "SPECIATION": "Speciation", "SPECIATIONS": "Speciation",
        "D": "Duplication", "DUPLICATION": "Duplication", "DUPLICATIONS": "Duplication",
        "T": "Transfer", "TRANSFER": "Transfer", "TRANSFERS": "Transfer",
        "L": "Loss", "LOSS": "Loss", "LOSSES": "Loss",
        "O": "Origination", "ORIGINATION": "Origination", "ORIGINATIONS": "Origination",
    }
    for path in paths:
        try:
            with open(path) as fh:
                for line in fh:
                    m = re.match(r"\s*([A-Za-z]+)\s*[=:]\s*([0-9]+)", line)
                    if not m:
                        continue
                    key = m.group(1).strip().upper()
                    if key in key_map:
                        counts[key_map[key]] += int(m.group(2))
                        found_any = True
        except OSError:
            continue
    return (counts if found_any else None)


def parse_nhx(paths):
    """Count D/H/T tags in GeneRax NHX reconciliations as a fallback.

    Only Duplication and Transfer are counted (these have unambiguous NHX tags).
    Speciation is NOT inferred from `D=0` -- a `D=0` node can be a transfer or
    origination, so counting every `D=0` as a speciation badly overcounts and would
    corrupt the vertical fraction. We therefore leave Speciation at 0 on this fallback
    path; the caller flags the vertical fraction as unavailable for NHX-only parses.
    """
    counts = {e: 0 for e in EVENT_TYPES}
    found_any = False
    for path in paths:
        try:
            with open(path) as fh:
                txt = fh.read()
        except OSError:
            continue
        # NHX comments look like [&&NHX:S=...:D=1:H=...]
        for token, etype in _NHX_EVENT.items():
            k, v = token.split("=")
            n = len(re.findall(rf"{k}={v}\b", txt))
            if n:
                counts[etype] += n
                found_any = True
    return (counts if found_any else None)


def vertical_fraction(counts):
    """Speciation / (Speciation + Transfer); '' when undefined."""
    spec = counts.get("Speciation", 0)
    tr = counts.get("Transfer", 0)
    denom = spec + tr
    if denom <= 0:
        return ""
    return round(spec / denom, 6)


def run_generax(args, rows):
    """Run GeneRax, parse events, append rows. Returns the counts dict (or None)."""
    generax = shutil.which("generax")
    if generax is None:
        for e in EVENT_TYPES:
            rows.append(("generax", e, 0, "", "generax not on PATH"))
        return None

    gdir = os.path.join(args.workdir, "generax")
    os.makedirs(gdir, exist_ok=True)
    map_path = os.path.join(args.workdir, "map.link")
    n_map = build_mapping(args.gene_tree, map_path)
    if n_map < MIN_TIPS:
        for e in EVENT_TYPES:
            rows.append(("generax", e, 0, "", f"only {n_map} mapped tips (<{MIN_TIPS})"))
        return None

    families = os.path.join(args.workdir, "families.txt")
    write_families(families, args.gene_tree, args.aln, map_path, subst_model="LG+G")

    cmd = [
        generax, "-f", families, "-s", os.path.abspath(args.species_tree),
        "-r", args.model, "-p", gdir, "--strategy", "SPR",
    ]
    utils.log("running: " + " ".join(cmd))
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if proc.returncode != 0:
            err = proc.stderr.decode(errors="replace").strip()
            utils.log(f"generax exit {proc.returncode}: {err[:800]}")
            for e in EVENT_TYPES:
                rows.append(("generax", e, 0, "", f"generax failed (exit {proc.returncode})"))
            return None
    except Exception as e:  # noqa: BLE001
        utils.log(f"generax failed: {e}")
        for ev in EVENT_TYPES:
            rows.append(("generax", ev, 0, "", f"generax exception: {e}"))
        return None

    events_files, nhx_files = find_generax_outputs(gdir)
    counts = parse_events_txt(events_files)
    note = "parsed *_events.txt"
    nhx_only = False
    if counts is None:
        counts = parse_nhx(nhx_files)
        note = "parsed NHX reconciliation (D/T only; speciation not counted)"
        nhx_only = True
    if counts is None:
        utils.log("generax ran but no parseable reconciliation found")
        for e in EVENT_TYPES:
            rows.append(("generax", e, 0, "", "generax ran; no parseable output"))
        return None

    # vertical fraction needs a reliable Speciation count; on the NHX-only fallback we
    # do not count speciations, so the fraction is left blank rather than misleading.
    frac = "" if nhx_only else vertical_fraction(counts)
    for e in EVENT_TYPES:
        f = frac if e in ("Speciation", "Transfer") else ""
        rows.append(("generax", e, counts.get(e, 0), f, note))
    return counts


def run_ale(args, rows):
    """Try ALEobserve then ALEml_undated; many installs lack ALE -> log + skip."""
    aleobs = shutil.which("ALEobserve")
    aleml = shutil.which("ALEml_undated")
    if aleobs is None or aleml is None:
        for e in EVENT_TYPES:
            rows.append(("ale", e, 0, "", "ALE not on PATH (skipped)"))
        return

    adir = os.path.join(args.workdir, "ale")
    os.makedirs(adir, exist_ok=True)
    ufboot = os.path.abspath(args.ufboot)
    if is_trivial(ufboot):
        for e in EVENT_TYPES:
            rows.append(("ale", e, 0, "", "ufboot trivial/empty (skipped)"))
        return

    # ALEobserve builds the .ale object next to the ufboot file.
    try:
        proc = subprocess.run([aleobs, ufboot], cwd=adir,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        if proc.returncode != 0:
            err = proc.stderr.decode(errors="replace").strip()
            utils.log(f"ALEobserve exit {proc.returncode}: {err[:400]}")
            for e in EVENT_TYPES:
                rows.append(("ale", e, 0, "", "ALEobserve failed"))
            return
    except Exception as e:  # noqa: BLE001
        utils.log(f"ALEobserve failed: {e}")
        for ev in EVENT_TYPES:
            rows.append(("ale", ev, 0, "", f"ALEobserve exception: {e}"))
        return

    ale_obj = ufboot + ".ale"
    if not os.path.exists(ale_obj):
        # ALEobserve usually writes <basename>.ale in cwd
        cand = os.path.join(adir, os.path.basename(ufboot) + ".ale")
        ale_obj = cand if os.path.exists(cand) else ale_obj

    try:
        proc = subprocess.run(
            [aleml, os.path.abspath(args.species_tree), ale_obj], cwd=adir,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        )
        if proc.returncode != 0:
            err = proc.stderr.decode(errors="replace").strip()
            utils.log(f"ALEml_undated exit {proc.returncode}: {err[:400]}")
            for e in EVENT_TYPES:
                rows.append(("ale", e, 0, "", "ALEml_undated failed"))
            return
    except Exception as e:  # noqa: BLE001
        utils.log(f"ALEml_undated failed: {e}")
        for ev in EVENT_TYPES:
            rows.append(("ale", ev, 0, "", f"ALEml_undated exception: {e}"))
        return

    # ALE writes a *.uml_rec / *.uTs file with expected D/T/L/S counts.
    counts = parse_ale_results(adir)
    if counts is None:
        for e in EVENT_TYPES:
            rows.append(("ale", e, 0, "", "ALE ran; no parseable result"))
        return
    frac = vertical_fraction(counts)
    for e in EVENT_TYPES:
        f = frac if e in ("Speciation", "Transfer") else ""
        rows.append(("ale", e, counts.get(e, 0), f, "parsed ALE result"))


def parse_ale_results(adir):
    """Parse ALE *.uml_rec for the 'rate of' / total D T L S summary line."""
    counts = {e: 0 for e in EVENT_TYPES}
    found = False
    if not os.path.isdir(adir):
        return None
    rec_files = []
    for root, _d, files in os.walk(adir):
        for f in files:
            if f.endswith((".uml_rec", ".ucons_tree", ".uTs")):
                rec_files.append(os.path.join(root, f))
    for path in rec_files:
        try:
            with open(path) as fh:
                for line in fh:
                    # e.g. "# of	 Duplications	Transfers	Losses ..." then a value line,
                    # or "rate of\t..." ; capture the labelled totals if present.
                    low = line.lower()
                    m = re.search(r"duplications?\D+([0-9.]+)", low)
                    if m:
                        counts["Duplication"] += int(round(float(m.group(1))))
                        found = True
                    m = re.search(r"transfers?\D+([0-9.]+)", low)
                    if m:
                        counts["Transfer"] += int(round(float(m.group(1))))
                        found = True
                    m = re.search(r"loss(?:es)?\D+([0-9.]+)", low)
                    if m:
                        counts["Loss"] += int(round(float(m.group(1))))
                        found = True
                    m = re.search(r"speciations?\D+([0-9.]+)", low)
                    if m:
                        counts["Speciation"] += int(round(float(m.group(1))))
                        found = True
        except OSError:
            continue
    return (counts if found else None)


def main():
    ap = argparse.ArgumentParser()
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
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # ---- guards ----
    if str(args.run_generax) != "True":
        write_header_only(args.out, "run_generax disabled")
        utils.log("run_generax != True; header-only output")
        return
    if is_trivial(args.species_tree) or is_trivial(args.gene_tree):
        write_header_only(args.out, "trivial species/gene tree")
        utils.log("trivial species/gene tree; header-only output")
        return
    n_sp = count_tips(args.species_tree)
    n_gt = count_tips(args.gene_tree)
    if n_sp < MIN_TIPS or n_gt < MIN_TIPS:
        write_header_only(args.out, f"too few tips (species={n_sp}, gene={n_gt})")
        utils.log(f"too few tips (species={n_sp}, gene={n_gt}); header-only output")
        return

    rows = []
    run_generax(args, rows)

    if str(args.run_ale) == "True":
        run_ale(args, rows)

    with open(args.out, "w") as out:
        out.write(HEADER + "\n")
        for tool, event, count, frac, note in rows:
            out.write(f"{tool}\t{event}\t{count}\t{frac}\t{note}\n")

    utils.log(f"wrote {len(rows)} reconciliation row(s) to {args.out}")


if __name__ == "__main__":
    main()
