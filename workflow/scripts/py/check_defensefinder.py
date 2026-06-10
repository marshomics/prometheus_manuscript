"""Preflight gate for the sole detector: verify DefenseFinder can actually find Prometheus.

Why this exists
---------------
Prometheus entered the DefenseFinder model registry in defense-finder-models v2.0.0
(absent from <= 1.2.4), and those models require the defense-finder SOFTWARE at
>= 2.0.0. An older binary, or a models set without Prometheus, detects ZERO
Prometheus loci no matter the input -- and because every detection step degrades to
header-only output, the whole pipeline would otherwise complete "green" with an empty
catalog and null results everywhere downstream. This rule turns that silent failure
into a hard, early stop.

Checks (run mode)
  1. software version: parse `defense-finder --version` and FAIL if < --min-software.
  2. model presence: confirm a 'Prometheus' model is installed (macsydata, else a
     filesystem scan of the installed defense-finder-models definitions). If a models
     directory is found but carries no Prometheus definition -> FAIL. If the models
     cannot be located by any method, WARN (the version gate already protects).

Ingest mode (--run False): no tool is used, so versions are not checked; we only
confirm the provided DefenseFinder output directory exists and is non-empty.

On success a sentinel file (--out) is written recording what was verified; the
per-genome `defensefinder` rule depends on it, so detection cannot start until the
detector is known to be capable of calling Prometheus.
"""
from __future__ import annotations
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import utils  # noqa: E402

SYSTEM_DEFAULT = "Prometheus"


def parse_semver(text: str):
    """Return the first (major, minor, patch) tuple found in text, or None."""
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text or "")
    if not m:
        return None
    return (int(m.group(1)), int(m.group(2)), int(m.group(3) or 0))


def ge(a, b) -> bool:
    """Tuple >= comparison with missing components treated as 0."""
    a = tuple(list(a) + [0, 0, 0])[:3]
    b = tuple(list(b) + [0, 0, 0])[:3]
    return a >= b


def software_version():
    """(version_tuple|None, raw_string)."""
    exe = shutil.which("defense-finder")
    if exe is None:
        return None, "defense-finder not on PATH"
    for args in (["--version"], ["version"]):
        try:
            res = subprocess.run([exe] + args, capture_output=True, text=True)
            raw = (res.stdout or "") + (res.stderr or "")
            v = parse_semver(raw)
            if v:
                return v, raw.strip().splitlines()[0] if raw.strip() else str(v)
        except Exception as exc:  # noqa: BLE001
            return None, f"defense-finder --version failed: {exc}"
    return None, "could not parse defense-finder version"


def model_present_macsydata(system: str):
    """Return True/False if macsydata can decide, else None (undetermined)."""
    if shutil.which("macsydata") is None:
        return None
    # `macsydata definition <package> <name>` prints the model and exits 0 if present.
    for pkg in ("defense-finder-models", "mdmparis/defense-finder-models"):
        try:
            res = subprocess.run(["macsydata", "definition", pkg, system],
                                 capture_output=True, text=True)
        except Exception:  # noqa: BLE001
            continue
        out = (res.stdout or "") + (res.stderr or "")
        if res.returncode == 0 and system.lower() in out.lower():
            return True
        # an explicit "not found / no such definition" with a zero-ish listing => absent
        if re.search(r"(not found|no such|unknown definition|does not exist)", out, re.I):
            return False
    return None


def model_present_filesystem(system: str):
    """Scan installed defense-finder-models definitions for the system; None if no dir."""
    roots = []
    for env in ("MACSY_MODELS_DIR", "DEFENSE_FINDER_MODELS"):
        if os.environ.get(env):
            roots.append(os.environ[env])
    # common install locations
    roots += [
        os.path.join(sys.prefix, "share"),
        os.path.join(os.path.expanduser("~"), ".macsyfinder", "models"),
        os.path.join(os.path.expanduser("~"), ".macsydata"),
    ]
    found_models_dir = False
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        hits = glob.glob(os.path.join(root, "**", "*efense*inder*models*"),
                         recursive=True)
        # any directory that looks like the models package
        for h in hits + [root]:
            defin = os.path.join(h, "definitions") if os.path.isdir(h) else ""
            search_root = defin if defin and os.path.isdir(defin) else h
            if not os.path.isdir(search_root):
                continue
            found_models_dir = True
            for ext in ("*.xml", "*.def"):
                for f in glob.glob(os.path.join(search_root, "**", ext), recursive=True):
                    if system.lower() in os.path.basename(f).lower():
                        return True
    return False if found_models_dir else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default="True", help="run_defensefinder (True/False)")
    ap.add_argument("--system", default=SYSTEM_DEFAULT)
    ap.add_argument("--min-software", default="2.0.0")
    ap.add_argument("--models-version", default="2.0.0")
    ap.add_argument("--provided-dir", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    notes = []
    run = str(args.run).strip().lower() in ("true", "1", "yes")
    min_sw = parse_semver(args.min_software) or (2, 0, 0)

    if not run:
        # Ingest mode: verify the provided outputs exist and are non-empty.
        pdir = args.provided_dir
        sys_files = glob.glob(os.path.join(pdir, "*", "*_defense_finder_systems.tsv")) if pdir else []
        if not sys_files:
            utils.log(f"FATAL: run_defensefinder=False but no provided DefenseFinder "
                      f"systems files found under '{pdir}'. Set prometheus."
                      f"defensefinder_provided_dir or run_defensefinder=true.")
            sys.exit(1)
        notes.append(f"ingest mode: {len(sys_files)} provided systems file(s) under {pdir}")
        _write_sentinel(args.out, notes)
        return

    # ---- 1. software version gate (the reliable, high-value check) -------- #
    v, raw = software_version()
    if v is None:
        utils.log(f"FATAL: could not determine defense-finder version ({raw}). "
                  f"Install defense-finder >= {args.min_software} (the Prometheus model "
                  f"needs software >= 2.0.0).")
        sys.exit(1)
    if not ge(v, min_sw):
        utils.log(f"FATAL: defense-finder {'.'.join(map(str, v))} < required "
                  f"{args.min_software}. The Prometheus model only exists in "
                  f"defense-finder-models v2.0.0+, which needs software >= 2.0.0. "
                  f"Detection would silently return ZERO Prometheus loci. Upgrade the "
                  f"env (workflow/envs/defensefinder.yaml).")
        sys.exit(1)
    notes.append(f"defense-finder version {'.'.join(map(str, v))} (>= {args.min_software})")

    # ---- 2. Prometheus model presence ------------------------------------ #
    decided = model_present_macsydata(args.system)
    method = "macsydata"
    if decided is None:
        decided = model_present_filesystem(args.system)
        method = "filesystem-scan"
    if decided is True:
        notes.append(f"'{args.system}' model installed (verified via {method})")
    elif decided is False:
        utils.log(f"FATAL: defense-finder is >= {args.min_software} but the "
                  f"'{args.system}' model is NOT in the installed model set "
                  f"(checked via {method}). Install models >= "
                  f"{args.models_version}:  defense-finder update  (or "
                  f"`macsydata install 'defense-finder-models>={args.models_version}'`).")
        sys.exit(1)
    else:
        utils.log(f"WARNING: could not locate the defense-finder model set to confirm "
                  f"the '{args.system}' model is present; relying on the version gate "
                  f"(software >= {args.min_software}). Run `defense-finder update` and "
                  f"confirm models >= {args.models_version}.")
        notes.append(f"'{args.system}' model presence UNVERIFIED (models dir not found); "
                     f"version gate passed")

    _write_sentinel(args.out, notes)


def _write_sentinel(path: str, notes) -> None:
    with open(path, "w") as fh:
        fh.write("# DefenseFinder preflight OK\n")
        for n in notes:
            fh.write(n + "\n")
    utils.log("DefenseFinder preflight passed: " + "; ".join(notes))


if __name__ == "__main__":
    main()
