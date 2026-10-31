#!/usr/bin/env bash
# Grid Engine job-status probe for Snakemake's cluster-generic executor.
# Prints exactly one of: running | success | failed
set -uo pipefail
JOBID="${1}"

# Still known to the scheduler (pending / running / held / error)?
if QSTAT=$(qstat -j "${JOBID}" 2>/dev/null); then
    # A job in error state (Eqw) stays in qstat forever; report it as failed instead
    # of letting the workflow wait on it indefinitely.
    if printf '%s\n' "${QSTAT}" | grep -qiE '^error reason|job_state[^:]*:[[:space:]]*E'; then
        echo "failed"; exit 0
    fi
    echo "running"; exit 0
fi

# Finished: consult accounting (it may lag a few seconds behind qstat). If the job id
# appears in several accounting records (id wrap-around), the LAST record is the current one.
ACCT=$(qacct -j "${JOBID}" 2>/dev/null || true)
EXIT=$(printf '%s\n' "${ACCT}" | awk '/^exit_status/{v=$2} END{print v}')
FAILED=$(printf '%s\n' "${ACCT}" | awk '/^failed/{v=$2} END{print v}')

if [ -z "${EXIT}" ]; then
    # Not yet in accounting: report running so Snakemake polls again. NOTE: if your
    # site has accounting disabled, qacct never answers and jobs never complete here;
    # use a status command based on your scheduler's job-finished records instead.
    echo "running"; exit 0
fi
if [ "${EXIT}" = "0" ] && [ "${FAILED:-0}" = "0" ]; then
    echo "success"
else
    echo "failed"
fi
