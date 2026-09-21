#!/usr/bin/env bash

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
"${FC:-gfortran}" "${geos_root}/GeosCore/brc_optics_mod.F90" \
    "${this_dir}/brc_optics_math_test.F90" \
    -o "${tmp_dir}/brc_optics_math_test"
"${tmp_dir}/brc_optics_math_test"
