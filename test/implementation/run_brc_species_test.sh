#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
test_tmp=$(mktemp -d)
trap 'rm -rf "${test_tmp}"' EXIT
"${FC:-gfortran}" -Wall -Wextra -fcheck=all \
    "${geos_root}/GeosCore/brc_species_mod.F90" \
    "${this_dir}/brc_species_test.F90" \
    -J "${test_tmp}" -I "${test_tmp}" -o "${test_tmp}/brc_species_test"
"${test_tmp}/brc_species_test"
