#!/usr/bin/env bash

set -euo pipefail

method1_shards="${METHOD1_SHARDS:-12}"
method2_shards="${METHOD2_SHARDS:-5}"
workers="${PARAM_WORKERS:-96}"
max_concurrent="${MAX_CONCURRENT:-${workers}}"
configurations=27
logical_tasks=$((configurations * method1_shards + 2 * configurations * method2_shards))

(( workers > 0 && workers <= logical_tasks )) || {
  echo "PARAM_WORKERS must be between 1 and ${logical_tasks}" >&2
  exit 2
}

export METHOD1_SHARDS="${method1_shards}"
export METHOD2_SHARDS="${method2_shards}"
export PARAM_WORKERS="${workers}"
export PARAM_STUDY_CASE="${PARAM_STUDY_CASE:-linear_progressive}"
export PARAM_STUDY_OUTPUT="${PARAM_STUDY_OUTPUT:-validation/results/parameter_study}"
export BEAM_METHOD1_SEEDS="${BEAM_METHOD1_SEEDS:-11,29,47,71,97}"
export BEAM_METHOD2_SEEDS="${BEAM_METHOD2_SEEDS:-1,2,3,4,5}"

log_directory="${PARAM_STUDY_LOGS:-validation/results/parameter_study_logs}"
mkdir -p "${log_directory}" "${PARAM_STUDY_OUTPUT}"

submission=$(bsub -q Batch24 \
  -J "beam_parameter_workers[1-${workers}]%${max_concurrent}" \
  -oo "${log_directory}/%J_%I.out" \
  -eo "${log_directory}/%J_%I.err" \
  -env all \
  < validation/lsf/parameter_study.lsf)

printf '%s\n' "${submission}"
printf 'Worker jobs: %s, maximum concurrent: %s\n' "${workers}" "${max_concurrent}"
printf 'Logical tasks distributed across workers: %s\n' "${logical_tasks}"
printf 'Case: %s, output: %s\n' "${PARAM_STUDY_CASE}" "${PARAM_STUDY_OUTPUT}"
