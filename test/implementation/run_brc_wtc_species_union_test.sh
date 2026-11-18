#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_core=$(cd "${this_dir}/../.." && pwd)
task_build=$(mktemp -d)
trap 'rm -rf "${task_build}"' EXIT
awk '
 /^  SUBROUTINE Unique_Species_Names\(/ { copying=1 }
 copying { print }
 /^  END SUBROUTINE Unique_Species_Names$/ { found=1; exit }
 END { if (!found) exit 1 }
' "${task_core}/Headers/species_database_mod.F90" > "${task_build}/wtc_species_union_under_test.inc"
"${FC:-gfortran}" -Wall -Wextra -fcheck=all -ffree-line-length-none \
    -I "${task_build}" -J "${task_build}" \
    "${this_dir}/brc_wtc_species_union_test.F90" -o "${task_build}/wtc_union_test"
"${task_build}/wtc_union_test"
