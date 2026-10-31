"""Preflight gate and output verification for the sole detector (DefenseFinder).

Why this exists
---------------
Prometheus entered the DefenseFinder model registry in defense-finder-models v2.0.0 and
those models need the defense-finder SOFTWARE >= 2.0.0. An older binary, an older model
set, a missing CasFinder package or a mismatched proteome all yield ZERO (or wrong)
Prometheus calls without any tool error, so the pipeline would otherwise finish "green"
on an empty catalog. Every check below FAILS (exit 1); nothing is downgraded to a warning.

Preflight (default mode) -- run mode (--run True)
  1. software: `defense-finder version` (defense_finder_cli/main.py prints
     "Using DefenseFinder version X"; there is NO --version flag) and the installed
     distribution `mdmparis-defense-finder` must agree, be >= --min-software and equal
     --software-version when given.
  2. models: located exactly as MacSyFinder 2.1.4 does (macsypy/config.py MacsyDefaults:
     --models-dir if given, else $VIRTUAL_ENV/share/macsyfinder/models or
     {/,/usr/local}/share/macsyfinder/models, plus ~/.macsyfinder/models, where
     `defense-finder update` = `macsydata install -U -u` puts them). Cross-checked against
     macsypy's own registry when macsypy is importable. Both `defense-finder-models` and
     `CasFinder` are required (DefenseFinder exits otherwise). The models version
     (metadata.yml `vers`) must be >= --models-version and equal --models-pin when given;
     CasFinder must equal --casfinder-pin when given. With macsyfinder 2.1.x every model
     definition must use grammar vers="2.0" (definition_parser.py); CasFinder 3.1.1 uses
     "2.1" and makes every DefenseFinder run crash.
  3. Prometheus model: definitions/DefenseFinder/<system>/<system>.xml must exist and
     declare --gene (default <system>__ProA) as mandatory; profiles/<gene>.hmm must exist.
     md5, LENG, GA and NSEQ are recorded for provenance.
  4. --positive-control true: `hmmemit -c` of that profile is searched with
     `defense-finder run --db-type gembase --skip-model-version-check`; a <system> call is
     required (proves the installed detector can call Prometheus at all).
  Ingest mode (--run False): the provided directory must exist and hold DefenseFinder
  output; versions cannot be verified and are recorded as UNVERIFIED.
  On success --out is written as a key<TAB>value provenance sentinel.

Ingest copy (--ingest-copy): copy one genome's provided tables, failing when they are
missing or ambiguous.

Output verification (--check-outputs): per genome, the systems/genes tables must have
the DefenseFinder v2/v3 columns, be mutually consistent, contain only protein ids present
in --faa (else: output produced on a different proteome), carry no v1 model paths
(DefenseFinder_1..5 = models < 2.0.0, no Prometheus), and -- with --expect-gembase --
replicon = hit_id up to its last '_' (contig-aware ids, Contract v2).
"""
from __future__ import annotations
import argparse
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SYSTEM_DEFAULT = "Prometheus"
SYSTEMS_COLUMNS = ["sys_id", "type", "subtype", "protein_in_syst"]
GENES_COLUMNS = ["replicon", "hit_id", "gene_name", "hit_pos", "model_fqn", "sys_id",
                 "hit_score"]


class PreflightError(RuntimeError):
    pass


def parse_semver(text: str):
    """First (major, minor, patch) in text, or None ('2.0.0+breaking...' -> (2, 0, 0))."""
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text or "")
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)), int(m.group(3) or 0))


def ge(a, b) -> bool:
    a = tuple(list(a) + [0, 0, 0])[:3]
    b = tuple(list(b) + [0, 0, 0])[:3]
    return a >= b


def vstr(v) -> str:
    return ".".join(map(str, v)) if v else "NA"


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------------- #
# 1. software version
# --------------------------------------------------------------------------- #
def software_version():
    """(cli_version|None, dist_version|None, raw_cli_text)."""
    cli, raw = None, ""
    exe = shutil.which("defense-finder")
    if exe:
        try:
            res = subprocess.run([exe, "version"], capture_output=True, text=True, timeout=120)
            raw = ((res.stdout or "") + (res.stderr or "")).strip()
            m = re.search(r"version\s+(\d+\.\d+(?:\.\d+)?)", raw, flags=re.I)
            if res.returncode == 0 and m:
                cli = parse_semver(m.group(1))
        except (OSError, subprocess.SubprocessError) as exc:
            raw = f"`defense-finder version` failed: {exc}"
    else:
        raw = "defense-finder not on PATH"
    dist = None
    try:
        from importlib.metadata import version, PackageNotFoundError
        try:
            dist = parse_semver(version("mdmparis-defense-finder"))
        except PackageNotFoundError:
            dist = None
    except ImportError:
        dist = None
    return cli, dist, raw


def macsyfinder_version() -> str:
    try:
        from importlib.metadata import version
        return version("macsyfinder")
    except Exception:  # noqa: BLE001 - provenance only
        return "NA"


# --------------------------------------------------------------------------- #
# 2. model packages (mirror of macsypy 2.1.4 Config/MacsyDefaults + ModelRegistry)
# --------------------------------------------------------------------------- #
def default_model_dirs() -> list:
    common = os.path.join("share", "macsyfinder", "models")
    venv = os.environ.get("VIRTUAL_ENV")
    if venv:
        system = os.path.join(venv, common)
    else:
        system = ""
        for prefix in ("/", os.path.join("/", "usr", "local")):
            p = os.path.join(prefix, common)
            if os.path.isdir(p):
                system = p          # the LAST existing prefix wins, as in macsypy
    user = os.path.join(os.path.expanduser("~"), ".macsyfinder", "models")
    return [p for p in (system, user) if p and os.path.exists(p)]


def read_vers(pkg_path: str):
    meta = os.path.join(pkg_path, "metadata.yml")
    if not os.path.exists(meta):
        return None
    with open(meta) as fh:
        for line in fh:
            m = re.match(r"^vers\s*:\s*(.+?)\s*$", line)
            if m:
                return m.group(1).strip().strip("'\"")
    return None


def installed_packages(models_dir: str = "") -> dict:
    """{package name: (path, vers)}; later directories override earlier ones."""
    dirs = [models_dir] if models_dir else default_model_dirs()
    if models_dir and not os.path.isdir(models_dir):
        raise PreflightError(f"--models-dir '{models_dir}' does not exist")
    pkgs = {}
    for d in dirs:
        for name in sorted(os.listdir(d)):
            p = os.path.join(d, name)
            if os.path.isdir(p) and os.path.isdir(os.path.join(p, "definitions")):
                pkgs[name] = (p, read_vers(p))
    return pkgs


def definition_versions(pkg_path: str) -> dict:
    """{model-definition grammar version (<model vers=...>): count} over a package.
    MacSyFinder 2.1.4 (definition_parser.py) only accepts vers="2.0"; CasFinder 3.1.1
    ships vers="2.1" definitions, so `defense-finder update` (latest CasFinder) breaks a
    macsyfinder 2.1.4 install -- the positive control below catches that."""
    out = {}
    for f in glob.glob(os.path.join(pkg_path, "definitions", "**", "*.xml"), recursive=True):
        try:
            v = ET.parse(f).getroot().get("vers", "missing")
        except ET.ParseError:
            v = "unparseable"
        out[v] = out.get(v, 0) + 1
    return out


def macsypy_packages(models_dir: str = ""):
    """Same listing from macsypy itself, or None when macsypy is not importable."""
    try:
        from macsypy.scripts.macsydata import _find_all_installed_packages
    except Exception:  # noqa: BLE001
        return None
    reg = _find_all_installed_packages(models_dir=models_dir or None)
    return {m.name: (m.path, m.version) for m in reg.models()}


# --------------------------------------------------------------------------- #
# 3. the Prometheus model definition + profile
# --------------------------------------------------------------------------- #
def find_definition(pkg_path: str, system: str):
    preferred = os.path.join(pkg_path, "definitions", "DefenseFinder", system, f"{system}.xml")
    if os.path.exists(preferred):
        return preferred
    found = glob.glob(os.path.join(pkg_path, "definitions", "**", f"{system}.xml"), recursive=True)
    return sorted(found)[0] if found else None


def inspect_model(pkg_path: str, system: str, gene: str) -> dict:
    xml = find_definition(pkg_path, system)
    if xml is None:
        raise PreflightError(f"no '{system}.xml' definition under {pkg_path}/definitions: the "
                             f"installed defense-finder-models cannot call {system}")
    root = ET.parse(xml).getroot()
    genes = {g.get("name"): g.get("presence") for g in root.iter("gene")}
    if gene not in genes:
        raise PreflightError(f"{xml} does not declare gene '{gene}' (genes: {sorted(genes)})")
    if genes[gene] != "mandatory":
        raise PreflightError(f"{xml}: gene '{gene}' is '{genes[gene]}', expected 'mandatory'")
    prof = os.path.join(pkg_path, "profiles", f"{gene}.hmm")
    if not os.path.exists(prof):
        raise PreflightError(f"profile {prof} for gene '{gene}' is missing")
    head = {}
    with open(prof) as fh:
        for line in fh:
            if line.startswith("HMM "):
                break
            k = line[:6].strip()
            if k in ("NAME", "LENG", "GA", "NSEQ", "DATE"):
                head[k] = line[6:].strip()
    return {"definition": xml, "definition_md5": md5(xml),
            "model_attrs": ";".join(f"{k}={v}" for k, v in sorted(root.attrib.items())),
            "genes": ",".join(f"{k}:{v}" for k, v in genes.items()),
            "profile": prof, "profile_md5": md5(prof),
            "profile_header": ";".join(f"{k}={v}" for k, v in head.items())}


# --------------------------------------------------------------------------- #
# 4. positive control
# --------------------------------------------------------------------------- #
def positive_control(profile: str, system: str, workdir: str, models_dir: str,
                     threads: int = 1) -> str:
    for tool in ("hmmemit", "defense-finder"):
        if shutil.which(tool) is None:
            raise PreflightError(f"positive control needs '{tool}' on PATH")
    os.makedirs(workdir, exist_ok=True)
    res = subprocess.run(["hmmemit", "-c", profile], capture_output=True, text=True)
    if res.returncode != 0 or not res.stdout.startswith(">"):
        raise PreflightError(f"hmmemit -c {profile} failed: {res.stderr[-500:]}")
    seq = "".join(res.stdout.splitlines()[1:])
    faa = os.path.join(workdir, "positive_control.faa")
    with open(faa, "w") as fh:          # gembase id: replicon 'pc_1', gene 1
        fh.write(">pc_1_1\n" + seq + "\n")
    out = os.path.join(workdir, "df_out")
    cmd = ["defense-finder", "run", "--db-type", "gembase", "--skip-model-version-check",
           "-w", str(max(1, int(threads))), "-o", out, faa]
    if models_dir:
        cmd[2:2] = ["--models-dir", models_dir]
    res = subprocess.run(cmd, capture_output=True, text=True)
    systems = os.path.join(out, "positive_control_defense_finder_systems.tsv")
    if res.returncode != 0 or not os.path.exists(systems):
        raise PreflightError("positive-control DefenseFinder run failed "
                             f"(exit {res.returncode}): {(res.stderr or res.stdout)[-1500:]}")
    with open(systems) as fh:
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    hdr = rows[0] if rows else []
    if "type" not in hdr:
        raise PreflightError(f"positive control: unexpected systems header {hdr}")
    ti, si = hdr.index("type"), hdr.index("subtype") if "subtype" in hdr else hdr.index("type")
    calls = [(r[ti], r[si]) for r in rows[1:] if len(r) > max(ti, si)]
    if not any(system.lower() in (t.lower(), s.lower()) for t, s in calls):
        raise PreflightError(f"positive control: DefenseFinder did NOT call {system} on the "
                             f"consensus of its own {os.path.basename(profile)} "
                             f"(calls: {calls})")
    return f"called {system} on hmmemit -c consensus (len {len(seq)}); calls={calls}"


# --------------------------------------------------------------------------- #
# output verification (--check-outputs)
# --------------------------------------------------------------------------- #
def check_outputs(systems: str, genes: str, faa: str, genome: str, expect_gembase: bool):
    import pandas as pd
    for p in (systems, genes, faa):
        if not os.path.exists(p):
            raise PreflightError(f"missing file {p}")
    sdf = pd.read_csv(systems, sep="\t", dtype=str, keep_default_na=False)
    gdf = pd.read_csv(genes, sep="\t", dtype=str, keep_default_na=False)
    for df, cols, p in ((sdf, SYSTEMS_COLUMNS, systems), (gdf, GENES_COLUMNS, genes)):
        miss = [c for c in cols if c not in df.columns]
        if miss:
            raise PreflightError(f"{p} lacks DefenseFinder v2 columns {miss} "
                                 f"(header: {list(df.columns)})")
    ids = set(utils.fasta_ids(faa))
    tag = f"genome {genome}: " if genome else ""
    bad = sorted(set(gdf["hit_id"]) - ids)
    in_sys = {p for v in sdf["protein_in_syst"] for p in str(v).split(",") if p}
    bad += sorted(in_sys - ids - set(bad))
    if bad:
        raise PreflightError(f"{tag}{len(bad)} DefenseFinder protein id(s) are absent from "
                             f"{faa} (e.g. {bad[:5]}): the output was produced on a different "
                             "proteome (re-run DefenseFinder on the pipeline's genes/*.faa)")
    s_ids, g_ids = set(sdf["sys_id"]), set(gdf["sys_id"])
    if s_ids != g_ids:
        raise PreflightError(f"{tag}systems and genes tables disagree on sys_id "
                             f"(only in systems: {sorted(s_ids - g_ids)[:3]}, only in genes: "
                             f"{sorted(g_ids - s_ids)[:3]})")
    v1 = sorted({m for m in gdf["model_fqn"] if re.search(r"/DefenseFinder_\d+/", m)})
    if v1:
        raise PreflightError(f"{tag}model paths {v1[:2]} are defense-finder-models v1 "
                             "(DefenseFinder_1..5): these models predate Prometheus (v2.0.0)")
    if expect_gembase and len(gdf):
        exp = gdf["hit_id"].str.rsplit("_", n=1).str[0]
        wrong = gdf.loc[exp != gdf["replicon"], ["hit_id", "replicon"]]
        if len(wrong):
            raise PreflightError(f"{tag}{len(wrong)} genes have replicon != contig prefix of "
                                 f"hit_id (e.g. {wrong.iloc[0].tolist()}): DefenseFinder was not "
                                 "run with --db-type gembase on contig-aware ids")
    return len(sdf), len(gdf)


def ingest_copy(provided_dir: str, genome: str, out_systems: str, out_genes: str) -> None:
    """Copy one genome's provided DefenseFinder tables (ingest mode). The exact names
    {g}/{g}_defense_finder_{systems,genes}.tsv are preferred; otherwise exactly one
    *_defense_finder_<kind>.tsv must exist in {provided}/{g}/. Missing or ambiguous
    files are errors (a silently empty table would read as 'no defense systems')."""
    gdir = os.path.join(provided_dir, genome)
    for kind, dst in (("systems", out_systems), ("genes", out_genes)):
        exact = os.path.join(gdir, f"{genome}_defense_finder_{kind}.tsv")
        if os.path.exists(exact):
            src = exact
        else:
            found = sorted(glob.glob(os.path.join(gdir, f"*_defense_finder_{kind}.tsv")))
            if len(found) != 1:
                raise PreflightError(f"ingest mode: expected {exact} (or exactly one "
                                     f"*_defense_finder_{kind}.tsv in {gdir}); found "
                                     f"{len(found)}: {found[:3]}")
            src = found[0]
        os.makedirs(os.path.dirname(os.path.abspath(dst)), exist_ok=True)
        shutil.copyfile(src, dst)
        utils.log(f"ingested {src} -> {dst}")


# --------------------------------------------------------------------------- #
def _write_sentinel(path: str, items) -> None:
    with open(path, "w") as fh:
        fh.write("# DefenseFinder preflight OK\n")
        for k, v in items:
            fh.write(f"{k}\t{v}\n")
    utils.log("DefenseFinder preflight passed: " + "; ".join(f"{k}={v}" for k, v in items))


def preflight(args) -> list:
    items = []
    run = str(args.run).strip().lower() in ("true", "1", "yes")
    gene = args.gene if args.gene is not None else f"{args.system}__ProA"
    if not run:
        pdir = args.provided_dir
        if not pdir or not os.path.isdir(pdir):
            raise PreflightError(f"run_defensefinder=False but provided dir '{pdir}' does "
                                 "not exist (prometheus.defensefinder_provided_dir)")
        sys_files = glob.glob(os.path.join(pdir, "*", "*_defense_finder_systems.tsv"))
        if not sys_files:
            raise PreflightError(f"no */*_defense_finder_systems.tsv under '{pdir}'")
        items += [("mode", "ingest"), ("provided_dir", pdir),
                  ("provided_systems_files", len(sys_files)),
                  ("defense_finder_version", "UNVERIFIED (ingest mode)"),
                  ("defense_finder_models_version", "UNVERIFIED (ingest mode)")]
        return items

    min_sw = parse_semver(args.min_software) or (2, 0, 0)
    cli, dist, raw = software_version()
    if cli is None and dist is None:
        raise PreflightError(f"cannot determine the defense-finder version ({raw})")
    if cli and dist and cli != dist:
        raise PreflightError(f"`defense-finder version` reports {vstr(cli)} but the installed "
                             f"mdmparis-defense-finder is {vstr(dist)}: PATH and env disagree")
    v = cli or dist
    if not ge(v, min_sw):
        raise PreflightError(f"defense-finder {vstr(v)} < required {args.min_software}: the "
                             "Prometheus model needs software >= 2.0.0")
    if args.software_version and parse_semver(args.software_version) != v:
        raise PreflightError(f"defense-finder {vstr(v)} != pinned {args.software_version} "
                             "(workflow/envs/defensefinder.yaml)")
    items += [("mode", "run"), ("defense_finder_version", vstr(v)),
              ("macsyfinder_version", macsyfinder_version())]

    pkgs = installed_packages(args.models_dir)
    mp = macsypy_packages(args.models_dir)
    if mp is not None:
        a = {k: (os.path.realpath(p), str(vv)) for k, (p, vv) in pkgs.items()}
        b = {k: (os.path.realpath(p), str(vv)) for k, (p, vv) in mp.items()}
        if a != b:
            raise PreflightError(f"model registry mismatch: this check sees {a}, macsypy sees {b}")
    items.append(("models_dirs", ",".join([args.models_dir] if args.models_dir
                                          else default_model_dirs()) or "NONE"))
    for name in ("defense-finder-models", "CasFinder"):
        if name not in pkgs:
            raise PreflightError(f"model package '{name}' is not installed (found: {sorted(pkgs)}). "
                                 "Install the pinned models, e.g. macsydata install --target "
                                 f"<dir> --org mdmparis 'defense-finder-models=={args.models_pin or args.models_version}' "
                                 "and macsydata install --target <dir> CasFinder, then set "
                                 "prometheus.defensefinder_models_dir")
    dpath, dvers = pkgs["defense-finder-models"]
    mv = parse_semver(dvers or "")
    if mv is None:
        raise PreflightError(f"{dpath}/metadata.yml has no parseable 'vers'")
    if not ge(mv, parse_semver(args.models_version) or (2, 0, 0)):
        raise PreflightError(f"defense-finder-models {dvers} < {args.models_version}: no "
                             f"{args.system} model before 2.0.0")
    if args.models_pin and parse_semver(args.models_pin) != mv:
        raise PreflightError(f"defense-finder-models {dvers} != pinned {args.models_pin}")
    cpath, cvers = pkgs["CasFinder"]
    if args.casfinder_pin and parse_semver(args.casfinder_pin) != parse_semver(cvers or ""):
        raise PreflightError(f"CasFinder {cvers} != pinned {args.casfinder_pin} (CasFinder 3.1.1 "
                             "definitions need a newer macsyfinder than 2.1.4)")
    msf = macsyfinder_version()
    dv = {"defense-finder-models": definition_versions(dpath), "CasFinder": definition_versions(cpath)}
    if msf.startswith("2.1."):
        bad = {k: v for k, v in dv.items() if set(v) - {"2.0"}}
        if bad:
            raise PreflightError(f"macsyfinder {msf} only parses model definitions with "
                                 f"vers=\"2.0\"; found {bad}. Install a compatible model "
                                 "version (e.g. CasFinder==3.1.0)")
    items += [("defense_finder_models_version", dvers), ("defense_finder_models_path", dpath),
              ("casfinder_version", cvers or "NA"), ("casfinder_path", cpath),
              ("definition_grammar_versions", ";".join(f"{k}:{v}" for k, v in dv.items()))]
    info = inspect_model(dpath, args.system, gene)
    items += [(f"{args.system.lower()}_{k}", val) for k, val in info.items()]
    if str(args.positive_control).strip().lower() in ("true", "1", "yes"):
        wd = args.workdir or os.path.join(os.path.dirname(os.path.abspath(args.out)),
                                          ".defensefinder_positive_control")
        items.append(("positive_control", positive_control(info["profile"], args.system, wd,
                                                           args.models_dir, threads=args.threads)))
    return items


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", default="True", help="run_defensefinder (True/False)")
    ap.add_argument("--system", default=SYSTEM_DEFAULT)
    ap.add_argument("--min-software", default="2.0.0")
    ap.add_argument("--models-version", default="2.0.0", help="minimum models version")
    ap.add_argument("--provided-dir", default="")
    ap.add_argument("--out", default="", help="sentinel (preflight mode)")
    # optional additions
    ap.add_argument("--software-version", default="", help="exact pinned software version")
    ap.add_argument("--models-pin", default="", help="exact pinned models version")
    ap.add_argument("--models-dir", default="", help="models root (defense-finder --models-dir)")
    ap.add_argument("--casfinder-pin", default="", help="exact pinned CasFinder version")
    ap.add_argument("--gene", default=None, help="mandatory gene (default <system>__ProA)")
    ap.add_argument("--positive-control", default="false")
    ap.add_argument("--workdir", default="")
    ap.add_argument("--threads", type=int, default=1,
                    help="DefenseFinder workers (-w) for the positive control")
    ap.add_argument("--check-outputs", action="store_true",
                    help="verify one genome's DefenseFinder tables instead of the preflight")
    ap.add_argument("--systems", default="")
    ap.add_argument("--genes", default="")
    ap.add_argument("--faa", default="")
    ap.add_argument("--genome", default="")
    ap.add_argument("--expect-gembase", action="store_true")
    ap.add_argument("--ingest-copy", action="store_true",
                    help="copy one genome's tables from --provided-dir to --systems/--genes")
    args = ap.parse_args()
    try:
        if args.ingest_copy:
            if not (args.provided_dir and args.genome and args.systems and args.genes):
                raise PreflightError("--ingest-copy needs --provided-dir --genome --systems --genes")
            ingest_copy(args.provided_dir, args.genome, args.systems, args.genes)
            return
        if args.check_outputs:
            n_s, n_g = check_outputs(args.systems, args.genes, args.faa, args.genome,
                                     args.expect_gembase)
            utils.log(f"DefenseFinder outputs OK for {args.genome or args.systems}: "
                      f"{n_s} systems, {n_g} genes")
            return
        if not args.out:
            raise PreflightError("--out is required in preflight mode")
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        _write_sentinel(args.out, preflight(args))
    except PreflightError as exc:
        utils.log(f"FATAL: {exc}")
        sys.exit(1)


if __name__ == "__main__":
    main()
