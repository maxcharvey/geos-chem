#!/usr/bin/env bash
set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
config_fullchem="${geos_root}/run/GCClassic/HEMCO_Config.rc.templates/HEMCO_Config.rc.fullchem"
config_aerosol="${geos_root}/run/GCClassic/HEMCO_Config.rc.templates/HEMCO_Config.rc.aerosol"
hco_interface="${geos_root}/GeosCore/hco_interface_gc_mod.F90"

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
    require "QFED_POG1_BRC_HS_PBL" "$config"
    require "311 QFED_PBL_FRAC 0.65" "$config"
    require "312 QFED_FT_FRAC  0.35" "$config"
done

require "L_QFED2 = GetExtNr" "$hco_interface"
require "QFED2_BRC_HARMONIZED_SENSITIVITY" "$hco_interface"
require "Select only one biomass-burning inventory: GFED, GFAS, QFED2, or FINNv25!" "$hco_interface"
require "QFED2_BRC_HARMONIZED_SENSITIVITY requires QFED2: true!" "$hco_interface"
echo "PASS: QFED2 BrC wiring static regression checks"
