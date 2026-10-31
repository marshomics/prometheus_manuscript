#!/bin/sh
# Post-deployment script of workflow/envs/rstats.yaml (POSIX sh).
#
# Snakemake (>= 6.14) runs it once, right after creating the conda env, as
#   source <conda>/bin/activate <env>; /bin/sh <env>.post-deploy.sh
# i.e. inside the activated env ($CONDA_PREFIX set, the env's R first on PATH), from a
# COPY of this file next to the env and with the workflow directory (the repository
# root) as the working directory.
#
# corHMM 2.8 (Beaulieu, O'Meara, Oliver & Boyko; CRAN, GPL >= 2, NeedsCompilation: no)
# provides fitCorrelationTest (Boyko & Beaulieu 2023 Syst Biol 72:476), used by
# methanogenesis_stats.R. conda-forge has no r-corhmm, so the CRAN source tarball is
# vendored at workflow/envs/vendor/corHMM_2.8.tar.gz and its dependencies are pinned in
# rstats.yaml. This script verifies the tarball's SHA-256, installs it into the env's R
# library with `R CMD INSTALL --no-docs --no-html`, and checks that the installed
# package loads and exports fitCorrelationTest. A missing tarball, a checksum mismatch,
# a missing or unwritable library, or any R error makes it exit non-zero, which fails
# the env deployment and therefore every rule that uses the env.
#
# Optional overrides (tests): PROM_CORHMM_TARBALL = tarball path; PROM_R_LIBRARY =
# target library (default $CONDA_PREFIX/lib/R/library).
set -eu

EXPECTED_SHA256=d6e46641f5b2b15720c8bb6c51d715220519dd8a6ba44171d13ef75bfa2eb6df
TARBALL_NAME=corHMM_2.8.tar.gz
PKG_VERSION=2.8

die() {
    echo "rstats.post-deploy: ERROR: $*" >&2
    exit 1
}
say() {
    echo "rstats.post-deploy: $*" >&2
}

# ---- the vendored source tarball ------------------------------------------------------
tarball=${PROM_CORHMM_TARBALL:-}
if [ -z "$tarball" ]; then
    here=$(CDPATH='' cd -- "$(dirname -- "$0")" 2>/dev/null && pwd) || here=''
    for cand in "${here:+$here/vendor/$TARBALL_NAME}" "workflow/envs/vendor/$TARBALL_NAME"; do
        if [ -n "$cand" ] && [ -f "$cand" ]; then
            tarball=$cand
            break
        fi
    done
fi
if [ -z "$tarball" ] || [ ! -f "$tarball" ]; then
    die "vendored $TARBALL_NAME not found (looked in vendor/ next to this script and in workflow/envs/vendor/ under $(pwd))"
fi

# ---- SHA-256 ---------------------------------------------------------------------------
if command -v sha256sum >/dev/null 2>&1; then
    got=$(sha256sum "$tarball" | awk '{print $1}')
elif command -v shasum >/dev/null 2>&1; then
    got=$(shasum -a 256 "$tarball" | awk '{print $1}')
elif command -v openssl >/dev/null 2>&1; then
    got=$(openssl dgst -sha256 "$tarball" | awk '{print $NF}')
else
    die "no SHA-256 tool (sha256sum, shasum or openssl) on PATH"
fi
if [ "$got" != "$EXPECTED_SHA256" ]; then
    die "SHA-256 mismatch for $tarball: expected $EXPECTED_SHA256, got ${got:-<empty>}"
fi

# ---- target library and R ----------------------------------------------------------------
if [ -n "${PROM_R_LIBRARY:-}" ]; then
    lib=$PROM_R_LIBRARY
else
    if [ -z "${CONDA_PREFIX:-}" ]; then
        die "CONDA_PREFIX is not set: run inside the activated rstats env (or set PROM_R_LIBRARY)"
    fi
    lib=$CONDA_PREFIX/lib/R/library
fi
[ -d "$lib" ] || die "R library directory not found: $lib"
[ -w "$lib" ] || die "R library directory is not writable: $lib"
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "$CONDA_PREFIX/bin/R" ]; then
    r_bin=$CONDA_PREFIX/bin/R
else
    r_bin=$(command -v R 2>/dev/null) || die "R not found (neither \$CONDA_PREFIX/bin/R nor R on PATH)"
fi

say "installing $tarball (SHA-256 verified) into $lib with $r_bin"
"$r_bin" CMD INSTALL --no-docs --no-html --library="$lib" "$tarball" ||
    die "R CMD INSTALL of $tarball into $lib failed"

# ---- the installed package must load and export fitCorrelationTest --------------------
PROM_POSTDEPLOY_LIB=$lib PROM_POSTDEPLOY_VERSION=$PKG_VERSION "$r_bin" --no-echo --no-save --no-restore -e '
lib <- Sys.getenv("PROM_POSTDEPLOY_LIB"); want <- Sys.getenv("PROM_POSTDEPLOY_VERSION")
v <- as.character(utils::packageVersion("corHMM", lib.loc = lib))
if (v != want) stop("installed corHMM ", v, " != ", want)
suppressPackageStartupMessages(library("corHMM", lib.loc = lib, character.only = TRUE))
if (!"fitCorrelationTest" %in% getNamespaceExports("corHMM"))
  stop("corHMM ", v, " does not provide fitCorrelationTest")
cat("rstats.post-deploy: corHMM", v, "installed and loadable from", lib, "\n")
' || die "corHMM installed into $lib but does not load (missing dependency?)"
