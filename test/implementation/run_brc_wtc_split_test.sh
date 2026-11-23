#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
test_tmp=$(mktemp -d)
trap 'rm -rf "${test_tmp}"' EXIT
python3 "${this_dir}/brc_wtc_split_test.py" --extract-native "${test_tmp}/actual_wtc.f90"
"${FC:-gfortran}" -Wall -Wextra -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow \
    -fopenmp "${test_tmp}/actual_wtc.f90" "${this_dir}/brc_wtc_split_driver.F90" \
    -J "${test_tmp}" -I "${test_tmp}" -o "${test_tmp}/brc_wtc_split_test"
OMP_NUM_THREADS=1 "${test_tmp}/brc_wtc_split_test"
OMP_NUM_THREADS=2 "${test_tmp}/brc_wtc_split_test"
