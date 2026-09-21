#!/usr/bin/env bash

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
tmp_dir=$(mktemp -d)
trap 'rm -rf "${tmp_dir}"' EXIT

"${FC:-gfortran}" "${this_dir}/brc_optics_math_test.F90" \
    -o "${tmp_dir}/brc_optics_math_test"
"${tmp_dir}/brc_optics_math_test"
