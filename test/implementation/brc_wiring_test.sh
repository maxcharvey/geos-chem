#!/usr/bin/env bash
#------------------------------------------------------------------------------
# Static regression checks for optional BrC wiring.
#
# This test deliberately does not build or run GEOS-Chem.  It guards the
# source-level contracts that let brown_carbon: false retain main-branch
# aerosol, photolysis, GFED-emission, and legacy-OC behaviour.
#------------------------------------------------------------------------------

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
repo_root=$(cd "${geos_root}/../.." && pwd)
hemco_root="${repo_root}/src/HEMCO"

require() {
    local pattern="${1}"
    local file="${2}"
    if ! rg -q --fixed-strings "${pattern}" "${file}"; then
        echo "FAIL: Missing '${pattern}' in ${file}" >&2
        exit 1
    fi
}

forbid() {
    local pattern="${1}"
    local file="${2}"
    if rg -q --fixed-strings "${pattern}" "${file}"; then
        echo "FAIL: Found stale '${pattern}' in ${file}" >&2
        exit 1
    fi
}

cldj="${geos_root}/GeosCore/cldj_interface_mod.F90"
carbon="${geos_root}/GeosCore/carbon_mod.F90"
aerosol="${geos_root}/GeosCore/aerosol_mod.F90"
photolysis="${geos_root}/GeosCore/photolysis_mod.F90"
hco_interface="${geos_root}/GeosCore/hco_interface_gc_mod.F90"
gfed="${hemco_root}/src/Extensions/hcox_gfed_mod.F90"

# Cloud-J stratospheric aerosol fields must follow the expanded NRHAER layout.
require "I_STRAT_AER_FIRST = 10 + NRHAER * NRH + 1" "${cldj}"
require "AERSP(L,I_STRAT_AER_FIRST)" "${cldj}"
require "AERSP(L,I_STRAT_AER_FIRST+1)" "${cldj}"
forbid "AERSP(L,41)" "${cldj}"
forbid "AERSP(L,42)" "${cldj}"

# brown_carbon must control chemistry, optical inputs, GFED partitioning, and
# return fossil-fuel OC to the main-branch OC tracers.
require "IF ( Input_Opt%LBRC ) THEN" "${carbon}"
require "IF ( Input_Opt%LBRC .AND. id_FFOCPO > 0 ) THEN" "${carbon}"
require "IF ( .NOT. LBRC ) SPECFIL(6:11) = \"org.dat  \"" "${aerosol}"
require "IF ( .NOT. Input_Opt%LBRC ) IND(6:NRHAER) = 36" "${photolysis}"

# The ordinary brown_carbon:false state has only the five legacy hygroscopic
# species.  Keep all 11 distinct optical bins valid; do not alias or drop
# inactive BrC slots.
require "Map_NRHAER(:) = (/ ( N, N = 1, NRHAER ) /)" "${aerosol}"
require "INTEGER :: Map_NRHAER(NRHAER)" "${aerosol}"
require "IF ( State_Chm%nHygGrth > NRHAER ) THEN" "${aerosol}"
require "IF ( Seen_NRHAER(Map_NRHAER(N)) ) THEN" "${aerosol}"
require "brown_carbon has duplicate hygroscopic species" "${aerosol}"
require "brown_carbon is missing canonical aerosol bin" "${aerosol}"
require "DO N = 1, State_Chm%nHygGrth" "${aerosol}"
require "DO NA = 1, State_Chm%nHygGrth" "${aerosol}"

# A longer Cloud-J table is usable for dry DBRC only when all five optional
# records are present and explicitly carry DBRC identities.
require "NAA >= DBRC_FJX_FIRST .AND. NAA < DBRC_FJX_LAST" "${photolysis}"
require "(/ 'DB00', 'DB50', 'DB70', 'DB80', 'DB90' /)" "${photolysis}"
require "MieTitle(1:4) /=" "${photolysis}"
require "Cloud-J DBRC optics: dedicated records" "${photolysis}"
require "Cloud-J DBRC optics: wet-BrC fallback records 57-61" "${photolysis}"

require "GEOSCHEM_BROWN_CARBON" "${hco_interface}"
require "GEOSCHEM_BROWN_CARBON" "${gfed}"
require "IF ( .NOT. Inst%UseBrC ) THEN" "${gfed}"
require "CALL Restore_Legacy_OC( State_Chm, HMRC )" "${hco_interface}"
require "SUBROUTINE Restore_Legacy_OC" "${hco_interface}"
require "CALL HCO_ArrAssert( HcoState%Spc(Hco_OCPI)%Emis" "${hco_interface}"
require "HcoState%Spc(Hco_OCPI)%Emis%Val +" "${hco_interface}"
require "HcoState%Spc(Hco_FFOCPI)%Emis%Val = 0.0_hp" "${hco_interface}"

echo "PASS: BrC wiring static regression checks"
