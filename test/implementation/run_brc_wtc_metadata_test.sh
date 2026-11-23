#!/usr/bin/env bash
set -euo pipefail
this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
task_core=$(cd "${this_dir}/../.." && pwd)
task_build=$(mktemp -d)
trap 'rm -rf "${task_build}"' EXIT
python3 - "${task_core}/Headers/species_database_mod.F90" "${task_build}/wtc_metadata_under_test.inc" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]).read_text()
a=s.index('       ! Is_Gas and Is_Aero tags cannot both be TRUE')
b=s.index('       ! Make sure the molecular weight',a)
Path(sys.argv[2]).write_text(s[a:b])
PY
"${FC:-gfortran}" -Wall -Wextra -fcheck=all -ffree-line-length-none \
    -I "${task_build}" -J "${task_build}" \
    "${this_dir}/brc_wtc_metadata_test.F90" -o "${task_build}/wtc_metadata_test"
"${task_build}/wtc_metadata_test"
