#!/usr/bin/env bash

set -euo pipefail

method1_shards="${METHOD1_SHARDS:-12}"
method2_shards="${METHOD2_SHARDS:-4}"
max_concurrent="${MAX_CONCURRENT:-8}"
total=$((3 * (method1_shards + method2_shards)))

mkdir -p validation/results/lsf_logs +         validation/results/lsf_shards +         validation/results/lsf_merged

submission=$(bsub +  -J "beam_validation[1-${total}]%${max_concurrent}" +  -oo "validation/results/lsf_logs/%J_%I.out" +  -eo "validation/results/lsf_logs/%J_%I.err" +  -env "all,METHOD1_SHARDS=${method1_shards},METHOD2_SHARDS=${method2_shards}" +  < validation/lsf/topology_study.lsf)
printf '%s\n' "${submission}"

job_id=$(printf '%s\n' "${submission}" |
  sed -n 's/Job <\([0-9][0-9]*\)>.*/\1/p')
test -n "${job_id}" || {
  echo "Could not parse the LSF job id; submit the merge job manually." >&2
  exit 1
}

bsub -w "done(${job_id})" +  -oo "validation/results/lsf_logs/merge_%J.out" +  -eo "validation/results/lsf_logs/merge_%J.err" +  < validation/lsf/merge_topology_study.lsf
