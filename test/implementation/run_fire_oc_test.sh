#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_build=$(mktemp -d)
trap 'rm -rf "${task_build}"' EXIT
python3 "${this_dir}/fire_oc_native_test.py" "${task_build}/actual_fire_oc.f90"
"${FC:-gfortran}" -Wall -Wextra -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow \
    -fopenmp "${task_build}/actual_fire_oc.f90" -J "${task_build}" -o "${task_build}/fire_oc_test"
OMP_NUM_THREADS=1 "${task_build}/fire_oc_test"
OMP_NUM_THREADS=2 "${task_build}/fire_oc_test"
