#!/usr/bin/env bash

set -euo pipefail

method1_shards="${METHOD1_SHARDS:-12}"
method2_shards="${METHOD2_SHARDS:-4}"
max_concurrent="${MAX_CONCURRENT:-8}"
total=$((3 * (method1_shards + method2_shards)))
run_id="$(date +%Y%m%d_%H%M%S)"
shard_output="${BEAM_STUDY_OUTPUT:-validation/results/lsf_shards_${run_id}}"
merged_output="${BEAM_STUDY_MERGED_OUTPUT:-${shard_output}_merged}"

mkdir -p validation/results/lsf_logs "${shard_output}" "${merged_output}"

queue_args=()
if [[ -n "${LSF_QUEUE:-}" ]]; then
  queue_args=(-q "$LSF_QUEUE")
fi

job_environment="all,METHOD1_SHARDS=${method1_shards},METHOD2_SHARDS=${method2_shards},BEAM_STUDY_OUTPUT=${shard_output},BEAM_STUDY_MERGED_OUTPUT=${merged_output}"
submission=$(bsub "${queue_args[@]}" -J "beam_validation[1-${total}]%${max_concurrent}" -oo "validation/results/lsf_logs/%J_%I.out" -eo "validation/results/lsf_logs/%J_%I.err" -env "${job_environment}" < validation/lsf/topology_study.lsf)
printf '%s\n' "${submission}"
printf 'Shard output: %s\nMerged output: %s\n' \
  "${shard_output}" "${merged_output}"

job_id=$(printf '%s\n' "${submission}" |
  sed -n 's/Job <\([0-9][0-9]*\)>.*/\1/p')
test -n "${job_id}" || {
  echo "Could not parse the LSF job id; submit the merge job manually." >&2
  exit 1
}

bsub "${queue_args[@]}" -w "done(${job_id})" -oo "validation/results/lsf_logs/merge_%J.out" -eo "validation/results/lsf_logs/merge_%J.err" -env "${job_environment}" < validation/lsf/merge_topology_study.lsf
