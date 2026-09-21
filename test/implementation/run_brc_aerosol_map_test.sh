#!/usr/bin/env bash

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

"${FC:-gfortran}" \
    "${geos_root}/GeosCore/brc_aerosol_map_mod.F90" \
    "${this_dir}/brc_aerosol_map_test.F90" \
    -J "${tmp_dir}" \
    -o "${tmp_dir}/brc_aerosol_map_test"

"${tmp_dir}/brc_aerosol_map_test"
