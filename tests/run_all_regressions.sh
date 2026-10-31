#!/usr/bin/env bash
# Every stage suite, the DAG check, rules hygiene and the environment solves.
# Prints ALL REGRESSIONS PASSED only when every one of them passes.
#   bash tests/run_all_regressions.sh
set -uo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here/.."
[ -f .unlazy/toolchain.env ] && . .unlazy/toolchain.env
PY="${PROM_PYTHON:-python3}"
fail=0
run() {  # run <label> <expected-token> <command...>
  local label="$1" token="$2"; shift 2
  local out; out="$("$@" 2>&1)"; local rc=$?
  if [ $rc -eq 0 ] && printf '%s' "$out" | grep -qF "$token"; then
    echo "ok   $label"
  else
    echo "FAIL $label (exit=$rc, expected '$token')"
    printf '%s\n' "$out" | tail -15
    fail=1
  fi
}
for stage in preprocess_qc detection_diversity scaffolds distribution evolution ecology \
             synthesis crossdomain disorder methanogenesis; do
  token="$(printf '%s' "$stage" | tr '[:lower:]' '[:upper:]') TESTS PASSED"
  run "suite $stage" "$token" bash tests/run_leaf.sh "$stage"
done
if [ -d tests/methanogenesis_stats ]; then
  run "suite methanogenesis_stats" "METHANOGENESIS STATS TESTS PASSED" \
      bash tests/run_leaf.sh methanogenesis_stats
fi
if [ -d tests/methanogenesis_benchmark ]; then
  run "suite methanogenesis_benchmark" "METHANOGENESIS BENCHMARK TESTS PASSED" \
      bash tests/run_leaf.sh methanogenesis_benchmark
fi
run "DAG" "DAG CHECK OK" bash tests/check_dag.sh
run "rules hygiene" "RULES HYGIENE OK" "$PY" tests/check_rules.py \
    workflow/rules/preprocess.smk workflow/rules/qc.smk workflow/rules/prometheus.smk \
    workflow/rules/diversity.smk workflow/rules/scaffolds.smk workflow/rules/distribution.smk \
    workflow/rules/evolution.smk workflow/rules/ecology.smk workflow/rules/synthesis.smk \
    workflow/rules/plots_extra.smk workflow/rules/crossdomain.smk workflow/rules/disorder.smk
run "owned files" "OWNED FILES CLEAN" "$PY" tests/methanogenesis/check_owned.py
run "env solve" "ENV SOLVE OK" bash tests/check_env_solve.sh workflow/envs/*.yaml
[ $fail -eq 0 ] && echo "ALL REGRESSIONS PASSED" || { echo "REGRESSIONS FAILED"; exit 1; }
