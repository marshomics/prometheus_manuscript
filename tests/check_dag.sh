#!/usr/bin/env bash
# Dry-run every DAG on the synthetic .test config. Prints DAG CHECK OK on success.
# Needs Snakemake >= 8 (PROM_SNAKEMAKE or `snakemake` on PATH). Run from the repo root.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; cd "$here/.."
if [ -f .unlazy/toolchain.env ]; then . .unlazy/toolchain.env; fi
SMK="${PROM_SNAKEMAKE:-snakemake}"; PY="${PROM_PYTHON:-python3}"
# One private scratch dir under an ABSOLUTE base: a relative TMPDIR (some sandboxes set
# ./.tmp) would resolve against the repo root and make mktemp fail.
base="${TMPDIR:-/tmp}"; case "$base" in /*) ;; *) base=/tmp ;; esac
scratch=$(mktemp -d "$base/check_dag.XXXXXX") || { echo "DAG CHECK FAILED: no scratch dir"; exit 1; }
trap 'rm -rf "$scratch"' EXIT
n_out=0; newout(){ n_out=$((n_out+1)); out="$scratch/dry_$n_out.log"; : > "$out"; }
cfg=.test/config/config.yaml; fail=0
jobs_total(){ awk '/^Job stats/{f=1} f && $1=="total" {print $2; exit}' "$1"; }
dry(){ # target expect_min_jobs label
  newout; "$SMK" "$1" -n --nolock -s workflow/Snakefile --configfile "$cfg" >"$out" 2>&1; rc=$?
  n=$(jobs_total "$out"); n=${n:-0}
  if [ $rc -ne 0 ] || [ "$n" -lt "$2" ]; then echo "FAIL $3: exit=$rc jobs=$n"; tail -n 25 "$out" | sed 's/^/    /'; fail=1
  else echo "ok   $3: $n jobs"; fi; }
"$PY" .test/setup_dryrun.py --remove
dry preprocess 5 "preprocess (Phase 1)"
dry all 5 "all before Phase 1 (falls back to Phase 1)"
newout; "$SMK" all -n --nolock -s workflow/Snakefile --configfile "$cfg" >"$out" 2>&1
if grep -qE "rule (validate_prometheus|species_tree):" "$out"; then echo "FAIL Phase-1 fallback: Phase-2 jobs scheduled before Phase 1"; fail=1; else echo "ok   Phase-1 fallback has no Phase-2 jobs"; fi
newout; "$SMK" preprocess -n --nolock -s workflow/Snakefile --configfile "$cfg" \
  --config 'plotting={"formats": ["tiff"]}' >"$out" 2>&1; rc=$?
if [ $rc -ne 0 ] && grep -qE "plotting.formats must be|.tiff. is not one of" "$out"; then echo "ok   unsupported plotting.formats rejected at parse time"
else echo "FAIL plotting.formats=[tiff] was not rejected (exit=$rc)"; fail=1; fi
"$PY" .test/setup_dryrun.py >/dev/null
dry all 40 "all (Phase 2)"
dry crossdomain 30 "crossdomain"
dry disorder 15 "disorder"
"$PY" .test/setup_dryrun.py --remove
[ $fail -eq 0 ] && echo "DAG CHECK OK" || { echo "DAG CHECK FAILED"; exit 1; }
