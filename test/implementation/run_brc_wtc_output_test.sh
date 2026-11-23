#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_core=$(cd "${this_dir}/../.." && pwd)
task_build=$(mktemp -d)
trap 'rm -rf "${task_build}"' EXIT
python3 "${this_dir}/brc_wtc_split_test.py" --extract-native "${task_build}/actual_wtc.f90"
"${FC:-gfortran}" -Wall -Wextra -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow \
    -fopenmp "${task_build}/actual_wtc.f90" \
    "${task_core}/Headers/brc_wtc_output_mod.F90" \
    "${this_dir}/brc_wtc_output_driver.F90" \
    -J "${task_build}" -I "${task_build}" -o "${task_build}/wtc_output_test"
OMP_NUM_THREADS=1 "${task_build}/wtc_output_test"
OMP_NUM_THREADS=2 "${task_build}/wtc_output_test"
