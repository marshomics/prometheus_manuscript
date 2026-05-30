#!/usr/bin/env bash
# Grid Engine job-status probe for Snakemake's cluster-generic executor.
# Prints exactly one of: running | success | failed
set -uo pipefail
JOBID="${1}"

# Still in the queue/running?
if qstat -j "${JOBID}" >/dev/null 2>&1; then
    echo "running"; exit 0
fi

# Finished: consult accounting for the exit status (qacct may lag a few seconds).
EXIT=$(qacct -j "${JOBID}" 2>/dev/null | awk '/^exit_status/{print $2; exit}')
FAILED=$(qacct -j "${JOBID}" 2>/dev/null | awk '/^failed/{print $2; exit}')

if [ -z "${EXIT}" ]; then
    # Not yet in accounting; treat as still running so Snakemake polls again.
    echo "running"; exit 0
fi
if [ "${EXIT}" = "0" ] && [ "${FAILED:-0}" = "0" ]; then
    echo "success"
else
    echo "failed"
fi
