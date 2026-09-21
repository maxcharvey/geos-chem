#!/usr/bin/env bash
set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
config_fullchem="${geos_root}/run/GCClassic/HEMCO_Config.rc.templates/HEMCO_Config.rc.fullchem"
config_aerosol="${geos_root}/run/GCClassic/HEMCO_Config.rc.templates/HEMCO_Config.rc.aerosol"
hco_interface="${geos_root}/GeosCore/hco_interface_gc_mod.F90"
diagn_fullchem="${geos_root}/run/GCClassic/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.fullchem"
diagn_aerosol="${geos_root}/run/GCClassic/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.aerosol"

require() {
    local pattern="$1" file="$2"
    rg -q --fixed-strings "$pattern" "$file" || { echo "FAIL: ${pattern} missing from ${file}" >&2; exit 1; }
}

for config in "$config_fullchem" "$config_aerosol"; do
    require "QFED2_BRC_HARMONIZED_SENSITIVITY : false" "$config"
    require "(((QFED2_BRC_HARMONIZED_SENSITIVITY" "$config"
    require "QFED_FSOAP_BRC_HS_PBL" "$config"
    require "QFED_FSOAP_BRC_HS_FT" "$config"
    require "QFED_DBRCPOA_BRC_HS_PBL" "$config"
    require "QFED_DBRCPOA_BRC_HS_FT" "$config"
    require "QFED_NPBRCPOA_BRC_HS_PBL" "$config"
    require "QFED_PBRCPOA_BRC_HS_FT" "$config"
    require "QFED_OCPI_BRC_HS_PBL" "$config"
    require "QFED_OCPO_BRC_HS_FT" "$config"
    require "311 QFED_PBL_FRAC 0.65" "$config"
    require "312 QFED_FT_FRAC  0.35" "$config"
done

for diagn in "$diagn_fullchem" "$diagn_aerosol"; do
    for species in FSOAP DBRCPOA NPBRCPOA PBRCPOA; do
        require "#InvQFED2_${species} ${species} 0 5 2 3" "$diagn"
        require "#InvQFED2_${species}Column ${species} 0 5 2 2" "$diagn"
    done
done

require "'QFED2', OptValBool=LTMP" "$hco_interface"
require "QFED2_BRC_HARMONIZED_SENSITIVITY" "$hco_interface"
require "Select only one biomass-burning inventory: GFED, GFAS, QFED2, or FINNv25!" "$hco_interface"
require "QFED2_BRC_HARMONIZED_SENSITIVITY requires QFED2: true!" "$hco_interface"
python3 "${this_dir}/qfed_brc_config_test.py" "$config_fullchem" "$config_aerosol"
if python3 "${this_dir}/qfed_brc_config_test.py" >/dev/null 2>&1; then
    echo "FAIL: qfed config checker accepted no templates" >&2
    exit 1
fi
echo "PASS: QFED2 BrC wiring static regression checks"
