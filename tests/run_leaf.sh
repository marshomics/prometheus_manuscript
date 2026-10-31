#!/usr/bin/env bash
# Run one stage's synthetic test suite:  bash tests/run_leaf.sh <stage-dir-name>
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
if [ -f "$here/../.unlazy/toolchain.env" ]; then . "$here/../.unlazy/toolchain.env"; fi
stage="${1:?usage: tests/run_leaf.sh <stage>}"
exec "${PROM_PYTHON:-python3}" "$here/$stage/run_tests.py"
