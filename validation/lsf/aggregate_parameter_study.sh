#!/usr/bin/env bash

set -euo pipefail

root="${1:-validation/results/parameter_study}"
output="${2:-${root}/aggregated}"

julia --project=validation validation/aggregate_parameter_study.jl   "${root}" "${output}"
