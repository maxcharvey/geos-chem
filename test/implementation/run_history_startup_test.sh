#!/usr/bin/env bash
# Tiny standalone regression: compile the actual History_Write with observed I/O.
# No NetCDF, model executable, model data, or full-model build is involved.
set -euo pipefail
task_core=$(cd "$(dirname "$0")/../.." && pwd)
task_source=${HISTORY_WRITE_SOURCE:-$task_core/History/history_mod.F90}
task_build=$(mktemp -d)
trap 'rm -rf "$task_build"' EXIT
awk '
  /^  SUBROUTINE History_Write\(/ { copying=1 }
  copying { print }
  /^  END SUBROUTINE History_Write$/ { found=1; exit }
  END { if (!found) exit 1 }
' "$task_source" > "$task_build/history_write_under_test.inc"
cd "$task_build"
timeout 30 "${FC:-gfortran}" -O0 -cpp -fopenmp -fcheck=all \
  -fbacktrace -ffree-line-length-none -Wall -Wno-unused-dummy-argument \
  -I"$task_build" \
  "$task_core/test/implementation/history_startup_stubs.F90" \
  "$task_core/test/implementation/history_startup_test.F90" \
  -o history_startup_test
OMP_NUM_THREADS=2 timeout 10 ./history_startup_test
