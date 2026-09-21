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
# An isolated GEOS-Chem worktree is not nested below the GCClassic wrapper.
# Keep the wrapper-relative path for ordinary CI, but allow its HEMCO sibling
# to be supplied explicitly when reviewing a core-only worktree.
hemco_root="${GEOSCHEM_HEMCO_ROOT:-${repo_root}/src/HEMCO}"

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
brc_map="${geos_root}/GeosCore/brc_aerosol_map_mod.F90"
brc_cloudj_map="${geos_root}/GeosCore/brc_cloudj_map_mod.F90"
brc_optics="${geos_root}/GeosCore/brc_optics_mod.F90"
photolysis="${geos_root}/GeosCore/photolysis_mod.F90"
fastjx="${geos_root}/GeosCore/fast_jx_mod.F90"
fullchem_config="${geos_root}/run/GCClassic/geoschem_config.yml.templates/geoschem_config.yml.fullchem"
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
require "BrC_Aerosol_Optics" "${aerosol}"
require "Dedicated BrC optical" "${aerosol}"
require "IF ( LBRC .AND. TRIM(Input_Opt%BrC_Aerosol_Optics) == 'DEDICATED' ) THEN" "${aerosol}"
require "CALL BUILD_BRC_CLOUDJ_MAP(" "${photolysis}"

# The ordinary brown_carbon:false state has only the five legacy hygroscopic
# species.  Keep all 11 distinct optical bins valid; do not alias or drop
# inactive BrC slots.
require "Map_NRHAER(:) = (/ ( N, N = 1, NRHAER ) /)" "${aerosol}"
require "INTEGER :: Map_NRHAER(NRHAER)" "${aerosol}"
require "IF ( State_Chm%nHygGrth > NRHAER ) THEN" "${aerosol}"
require "CALL VALIDATE_BRC_AEROSOL_MAP(" "${aerosol}"
require "IF ( Seen(Bins(N)) ) THEN" "${brc_map}"
require "PURE INTEGER FUNCTION BRC_FASTJX_RH_BIN" "${brc_map}"
require "brown_carbon has duplicate hygroscopic species" "${aerosol}"
require "brown_carbon is missing canonical aerosol bin" "${aerosol}"
require "DO N = 1, State_Chm%nHygGrth" "${aerosol}"
require "DO NA = 1, State_Chm%nHygGrth" "${aerosol}"

# BrC PM accounting covers every BrC carrier.  DBRCPOA stays dry; WTC is
# reported with primary OC because it is the bleached primary-BrC product.
require "State_Chm%AerMass%BRCPI(I,J,L) * ORG_GROWTH" "${aerosol}"
require "State_Chm%AerMass%NPBRC(I,J,L) * ORG_GROWTH" "${aerosol}"
require "State_Chm%AerMass%WTCPI(I,J,L) * ORG_GROWTH" "${aerosol}"
require "State_Chm%AerMass%FSOAS(I,J,L) * ORG_GROWTH" "${aerosol}"
require "State_Chm%AerMass%PBRC(I,J,L)  * ORG_GROWTH" "${aerosol}"
require "+ State_Chm%AerMass%BRCPO(I,J,L)" "${aerosol}"
require "WTC is classified as primary OC here" "${aerosol}"

# The DBRC optical carrier is dry through all RRTMG properties, and its
# diagnostic is emitted after wavelength interpolation.
require "IF ( N == 11 ) THEN" "${aerosol}"
require "RTSSAER(I,J,L,IWV,NRT)  = SSAA(IWV,1,N,State_Chm%Phot%DRg)" "${aerosol}"
require "RTASYMAER(I,J,L,IWV,NRT) = ASYMAA(IWV,1,N,State_Chm%Phot%DRg)" "${aerosol}"
require "USE BRC_Optics_Mod, ONLY : BRC_DRY_AOD_AT_WAVELENGTH" "${aerosol}"
require "BrCDryAOD = BRC_DRY_AOD_AT_WAVELENGTH" "${aerosol}"
require "State_Diag%BrCDryAODWL1" "${aerosol}"
require "Returning zero also" "${brc_optics}"
require "BRC_DRY_AOD_AT_WAVELENGTH = 0.0D0" "${brc_optics}"

# Fast-JX's online-LUT solver reads the aerosol/RH slot directly.  Route the
# dry DBRC carrier to IR=1 before it reaches that solver.
fjx_interface="${geos_root}/GeosCore/fjx_interface_mod.F90"
require "USE BRC_AEROSOL_MAP_MOD, ONLY : BRC_FASTJX_RH_BIN" "${fjx_interface}"
require "IR_OPT = BRC_FASTJX_RH_BIN( N, IRHARR(NLON,NLAT,L) )" "${fjx_interface}"
require "IOPT,    IR_OPT, J" "${fjx_interface}"

# The family mass diagnostics are consistently carbon mass: FSOAS is OM and
# must be converted before aggregation; flux diagnostics retain their own
# stoichiometric units.
brc="${geos_root}/GeosCore/brc_mod.F90"
state_diag="${geos_root}/Headers/state_diag_mod.F90"
require "USE BRC_Species_Mod, ONLY : BRC_REQUIRED_SPECIES_ERROR" "${brc}"
require "ErrMsg = BRC_REQUIRED_SPECIES_ERROR( (/ id_FSOAS, id_BRCSOA, id_WTC /) )" "${brc}"
require "IF ( LEN_TRIM( ErrMsg ) > 0 ) THEN" "${brc}"
forbid "IF ( id_FSOAS  <= 0 ) RETURN" "${brc}"
forbid "IF ( id_BRCSOA <= 0 ) RETURN" "${brc}"
forbid "IF ( id_WTC    <= 0 ) RETURN" "${brc}"
forbid "tau ~ 0.25 day" "${brc}"
forbid "brc.dat" "${cldj}"
require "portable organic placeholder" "${aerosol}"
require "State_Chm%Species(id_FSOAS )%Conc(I,J,L) / OMOC_BBOA" "${brc}"
require "Absorbing BrC-family carbon mass" "${state_diag}"
require "Units = 'kgC'" "${state_diag}"

# The stratospheric aerosol records must follow NRHAER, not the legacy five
# hygroscopic bins.
require "NRHAER+1" "${cldj}"
require "NRHAER+2" "${cldj}"
require "LSTRATOD controls every optical consumer" "${aerosol}"
forbid "ODAER(:,:,:,:,NRHAER+1) = 0.d0" "${aerosol}"

# Cloud-J organic mode maps wet BrC to OC records and dry DBRC to OC00.
# Dedicated mode requires identified wet, persistent, and dry BrC records.
require "CASE ( 'ORGANIC' )" "${brc_cloudj_map}"
require "AerMap(11,J) = DRY_BRC_RECORD" "${brc_cloudj_map}"
require "(/ 'WB00', 'WB50', 'WB70', 'WB80', 'WB90' /)" "${brc_cloudj_map}"
require "(/ 'PB00', 'PB50', 'PB70', 'PB80', 'PB90' /)" "${brc_cloudj_map}"
require "IF ( MieTitle(1:4) /= 'DB00' ) THEN" "${brc_cloudj_map}"
require "Cloud-J BrC optics: organic-equivalence records" "${photolysis}"
forbid "wet-BrC fallback records 57-61" "${photolysis}"

# FAST-JX receives every expanded hygroscopic bin.  Its v2024-05 table
# reads only 56 records, so BrC must retain the portable OC-equivalence map.
require "DO M=1,NRHAER" "${fastjx}"
forbid "DO M=1,5" "${fastjx}"
forbid "LBRC" "${fastjx}"
require "FAST-JX BrC optics: organic-equivalence OC records 36-40" "${photolysis}"
require "fast-jx:" "${fullchem_config}"
require "FAST_JX/v2024-05/" "${fullchem_config}"

require "GEOSCHEM_BROWN_CARBON" "${hco_interface}"
require "GEOSCHEM_BROWN_CARBON" "${gfed}"
require "IF ( .NOT. Inst%UseBrC ) THEN" "${gfed}"
require "CALL Restore_Legacy_OC( State_Chm, HMRC )" "${hco_interface}"
require "SUBROUTINE Restore_Legacy_OC" "${hco_interface}"
require "CALL HCO_ArrAssert( HcoState%Spc(Hco_OCPI)%Emis" "${hco_interface}"
require "HcoState%Spc(Hco_OCPI)%Emis%Val +" "${hco_interface}"
require "HcoState%Spc(Hco_FFOCPI)%Emis%Val = 0.0_hp" "${hco_interface}"
require "IF ( .NOT. ASSOCIATED( HcoState ) ) RETURN" "${hco_interface}"
require "Select only one biomass-burning inventory" "${hco_interface}"
require "FINNv25_Inject requires FINNv25: true" "${hco_interface}"
require "A BrC harmonized fire sensitivity requires brown_carbon: true" "${hco_interface}"

echo "PASS: BrC wiring static regression checks"
