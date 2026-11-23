#!/usr/bin/env bash
# Small standalone NetCDF/kernel checks; invoke in the project's build/test allocation.
set -euo pipefail
task_repo=$(cd "$(dirname "$0")/../.." && pwd)
if [[ -d "$task_repo/GeosCore" ]]; then
  task_core="$task_repo/GeosCore"
else
  task_core="$task_repo/src/GEOS-Chem/GeosCore"
fi
task_build=$(mktemp -d)
trap 'rm -rf "$task_build"' EXIT
read -r -a task_nf_flags <<< "$("${NF_CONFIG:-nf-config}" --fflags)"
read -r -a task_nf_libs <<< "$("${NF_CONFIG:-nf-config}" --flibs)"
cd "$task_build"
"${FC:-gfortran}" -O2 -cpp -fopenmp -fcheck=all -ffree-line-length-none \
  "${task_nf_flags[@]}" -I"$task_core" -c "$task_core/restart_validation_mod.F90"
"${FC:-gfortran}" -O2 -cpp -fopenmp -fcheck=all -ffree-line-length-none \
  -I"$task_core" "$task_repo/test/implementation/restart_pressure_cache_test.F90" \
  -o pressure_cache_test
./pressure_cache_test
"${FC:-gfortran}" -O2 -cpp -fopenmp -fcheck=all -ffree-line-length-none \
  "${task_nf_flags[@]}" -I"$task_core" "$task_repo/test/implementation/restart_wrapper_test.F90" \
  restart_validation_mod.o "${task_nf_libs[@]}" -o wrapper_test
./wrapper_test "$@"
