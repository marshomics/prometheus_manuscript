#!/usr/bin/env bash
# Translate Snakemake resources into a Grid Engine qsub submission.
# Args (Snakemake appends the jobscript path as the final argument):
#   $1 mem_mb (TOTAL)  $2 runtime (MIN)  $3 threads  $4 rule  $5 jobid  $6 jobscript
set -euo pipefail

MEM_MB="${1}"; RUNTIME_MIN="${2}"; THREADS="${3}"; RULE="${4}"; JOBID="${5}"; JOBSCRIPT="${6}"

# On most SGE installs h_vmem is enforced PER SLOT, so divide the total by threads.
# If YOUR cluster treats h_vmem as per-job, drop the division (see README).
PERSLOT_MB=$(( (MEM_MB + THREADS - 1) / THREADS ))

# h_rt wants HH:MM:SS; Snakemake gives minutes.
H=$(( RUNTIME_MIN / 60 )); M=$(( RUNTIME_MIN % 60 ))
HRT=$(printf "%02d:%02d:00" "${H}" "${M}")

# Parallel environment name varies by site (qconf -spl): smp, make, mpi, ...
PE_NAME="${SGE_PE_NAME:-smp}"

mkdir -p logs/cluster

qsub -cwd -V -terse \
     -N "smk.${RULE}" \
     -pe "${PE_NAME}" "${THREADS}" \
     -l h_vmem="${PERSLOT_MB}M",h_rt="${HRT}" \
     -o "logs/cluster/${RULE}.${JOBID}.out" \
     -e "logs/cluster/${RULE}.${JOBID}.err" \
     "${JOBSCRIPT}"
