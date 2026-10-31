#!/usr/bin/env bash
# Dry-run solve conda env files for linux-64 (the cluster platform) and verify pip pins
# exist on PyPI.  usage: bash tests/check_env_solve.sh workflow/envs/a.yaml [b.yaml ...]
# Needs micromamba (MICROMAMBA or on PATH). Prints ENV SOLVE OK when every file solves.
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
if [ -f "$here/../.unlazy/toolchain.env" ]; then . "$here/../.unlazy/toolchain.env"; fi
MM="${MICROMAMBA:-$(command -v micromamba || true)}"
[ -n "$MM" ] || { echo "micromamba not found"; exit 2; }
export CONDA_OVERRIDE_GLIBC="${CONDA_OVERRIDE_GLIBC:-2.28}"
# micromamba needs a WRITABLE package cache and root prefix; when the configured ones are
# not writable (e.g. another user's workspace) it silently falls back to ./cache in the
# current directory (the repo). Pick the first writable location outside the repo instead.
_writable() { mkdir -p "$1" 2>/dev/null && [ -w "$1" ]; }
for _c in "${CONDA_PKGS_DIRS:-}" "${TMPDIR:+$(cd "${TMPDIR}" 2>/dev/null && pwd)/promsolve_pkgs}" "${HOME:-}/.promsolve_pkgs" "/tmp/promsolve_pkgs_$(id -u)"; do
  if [ -n "$_c" ] && _writable "$_c"; then export CONDA_PKGS_DIRS="$_c"; break; fi
done
for _r in "${MAMBA_ROOT_PREFIX:-}" "${CONDA_PKGS_DIRS%/*}/promsolve_root" "/tmp/promsolve_root_$(id -u)"; do
  if [ -n "$_r" ] && _writable "$_r"; then export MAMBA_ROOT_PREFIX="$_r"; break; fi
done
case "$CONDA_PKGS_DIRS" in "$(cd "$here/.." && pwd)"*) echo "refusing to use a package cache inside the repo: $CONDA_PKGS_DIRS"; exit 2;; esac
fail=0
for f in "$@"; do
  name="solve_$(basename "$f" .yaml)_$$"
  if out=$("$MM" create -y -n "$name" --platform linux-64 --dry-run -f "$f" 2>&1); then
    echo "solve ok : $f"
  else
    echo "SOLVE FAIL: $f"; echo "$out" | grep -E "does not exist|requires|conflict|incompatible|critical" | head -8 | sed 's/^/    /'; fail=1
  fi
  # pip pins (lines like '      - pkg==X.Y' under a pip: block)
  for spec in $(awk '/- pip:/{p=1;next} p && /^ *#/{next} p && /^ *- /{gsub(/^ *- */,""); sub(/ *#.*/,""); print; next} p && !/^ *$/{p=0}' "$f" | grep '==' || true); do
    pkg="${spec%%==*}"; ver="${spec##*==}"
    code=$(curl -s -o /dev/null -w "%{http_code}" -m 30 "https://pypi.org/pypi/$pkg/$ver/json")
    if [ "$code" = "200" ]; then echo "pip ok   : $spec"; else echo "PIP MISSING: $spec (PyPI http $code)"; fail=1; fi
  done
done
[ "$fail" = "0" ] && echo "ENV SOLVE OK" || { echo "ENV SOLVE FAILED"; exit 1; }
