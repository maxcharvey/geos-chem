#!/usr/bin/env bash
# Static contracts for the default-off FINNv2.5 BrC proxy sensitivity.

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
template_dir="${geos_root}/run/GCClassic"
fullchem_config="${template_dir}/HEMCO_Config.rc.templates/HEMCO_Config.rc.fullchem"
aerosol_config="${template_dir}/HEMCO_Config.rc.templates/HEMCO_Config.rc.aerosol"
fullchem_diagn="${template_dir}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.fullchem"
aerosol_diagn="${template_dir}/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.aerosol"
aerosol_geoschem="${template_dir}/geoschem_config.yml.templates/geoschem_config.yml.aerosol"

require() {
    local pattern="${1}"
    local file="${2}"
    if ! rg -q --fixed-strings -- "${pattern}" "${file}"; then
        echo "FAIL: missing '${pattern}' in ${file}" >&2
        exit 1
    fi
}

require_route() {
    local file="${1}" name="${2}" species="${3}" factors="${4}" mode="${5}"
    awk -v name="${name}" -v species="${species}" -v factors="${factors}" -v mode="${mode}" '
        /^\(\(\(\.not\.FINNV25_BRC_HARMONIZED_SENSITIVITY$/ { branch="off" }
        /^\(\(\(FINNV25_BRC_HARMONIZED_SENSITIVITY$/ { branch="on" }
        /^\)\)\)(\.not\.)?FINNV25_BRC_HARMONIZED_SENSITIVITY$/ { branch="" }
        $2 == name && branch == mode {
            count++
            if ($(NF-3) != species || $(NF-2) != factors ||
                (factors == "-" ? $3 != "0" : $3 == "0")) bad=1
        }
        END { if (count != 1 || bad) exit 1 }
    ' "${file}" || {
        echo "FAIL: ${mode} route ${name} to ${species} with scales ${factors} in ${file}" >&2
        exit 1
    }
}

factor() {
    local file="${1}"
    local id="${2}"
    awk -v id="${id}" '$1 == id { print $3; exit }' "${file}"
}

check_sensitivity_condition_balance() {
    local file="${1}"
    awk '
        /^\(\(\(FINNv25$/ { in_finn=1 }
        in_finn && /^\(\(\(\.not\.FINNv25_Inject$/ { surface_open++ }
        in_finn && /^\)\)\)\.not\.FINNv25_Inject$/ { surface_close++ }
        in_finn && /^\(\(\(FINNv25_Inject$/ { inject_open++ }
        in_finn && /^\)\)\)FINNv25_Inject$/ { inject_close++ }
        in_finn && /^\(\(\(\.not\.FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            off_open++
        }
        in_finn && /^\)\)\)\.not\.FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            off_close++
        }
        in_finn && /^\(\(\(FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            on_open++
        }
        in_finn && /^\)\)\)FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            on_close++
        }
        in_finn && /^\)\)\)FINNv25$/ {
            # Each surface and injection path gates SOAP and residual OC
            # in separate sensitivity pairs.
            if (surface_open != 1 || surface_close != 1 ||
                inject_open != 1 || inject_close != 1 ||
                off_open != 4 || off_close != 4 ||
                on_open != 4 || on_close != 4) bad=1
            in_finn=0
        }
        END {
            if (in_finn || surface_open != 1 || surface_close != 1 ||
                inject_open != 1 || inject_close != 1 ||
                off_open != 4 || off_close != 4 ||
                on_open != 4 || on_close != 4 || bad) exit 1
        }
    ' "${file}" || {
        echo "FAIL: unbalanced FINNv2.5 sensitivity conditions in ${file}" >&2
        exit 1
    }
}

check_config() {
    local config="${1}"
    local name name_scales scales

    # Opt-in only; ordinary FINNv2.5 OC mapping remains explicit.
    require "--> FINNV25_BRC_HARMONIZED_SENSITIVITY : false" "${config}"
    require "165     FINNv25_Inject          : off" "${config}"
    require "(((.not.FINNV25_BRC_HARMONIZED_SENSITIVITY" "${config}"
    require_route "${config}" FINNv25_OCPI OCPI 75/72 off
    require_route "${config}" FINNv25_OCPO OCPO 75/73 off
    require_route "${config}" FINNV25_INJECT_OCPI OCPI 75/72 off
    require_route "${config}" FINNV25_INJECT_OCPO OCPO 75/73 off

    require_route "${config}" FINNv25_SOAP SOAP 75/281 off
    require_route "${config}" FINNV25_INJECT_SOAP SOAP 75/281 off
    require_route "${config}" FINNv25_SOAP SOAP - on
    require_route "${config}" FINNV25_INJECT_SOAP SOAP - on

    # The controlled branches route residual OC to the fire-specific tracers.
    require_route "${config}" FINNv25_FOCPI_BRC_HS FOCPI 75/282/72 on
    require_route "${config}" FINNv25_FOCPO_BRC_HS FOCPO 75/282/73 on
    require_route "${config}" FINNV25_INJECT_FOCPI FOCPI 75/282/72 on
    require_route "${config}" FINNV25_INJECT_FOCPO FOCPO 75/282/73 on
    for name in FOCPI FOCPO FSOAP DBRCPOA NPBRCPOA PBRCPOA; do
        require_route "${config}" "FINNV25_INJECT_${name}" "${name}" - off
    done
    for name in OCPI OCPO; do
        require_route "${config}" "FINNV25_INJECT_${name}" "${name}" - on
    done

    # Both controlled branches use the same proxy factors.
    for name_scales in FSOAP:75/286/287 DBRCPOA:75/285 NPBRCPOA:75/283 PBRCPOA:75/284; do
        name="${name_scales%%:*}"
        scales="${name_scales#*:}"
        require_route "${config}" "FINNv25_${name}_BRC_HS" "${name}" "${scales}" on
        require_route "${config}" "FINNV25_INJECT_${name}" "${name}" "${scales}" on
    done

    for id in 282 283 284 285 286 287; do
        if [[ $(awk -v id="${id}" '$1 == id { n++ } END { print n+0 }' "${config}") -ne 1 ]]; then
            echo "FAIL: scale ID ${id} not unique in ${config}" >&2
            exit 1
        fi
        require "/${id}" "${config}"
    done

    oc_remain=$(factor "${config}" 282)
    oc_np=$(factor "${config}" 283)
    oc_p=$(factor "${config}" 284)
    dbrc=$(factor "${config}" 285)
    fsoap=$(factor "${config}" 286)
    mw_corr=$(factor "${config}" 287)

    [[ ${oc_remain} == "0.5" ]]
    [[ ${oc_np} == "0.375" ]]
    [[ ${oc_p} == "0.125" ]]
    [[ ${dbrc} == "4.0" ]]
    [[ ${fsoap} == "0.013" ]]
    [[ ${mw_corr} == "0.186733333333" ]]

    awk -v a="${oc_remain}" -v b="${oc_np}" -v c="${oc_p}" \
        'BEGIN { x=a+b+c; if (x < 0.999999999 || x > 1.000000001) exit 1 }' || {
        echo "FAIL: FINNv2.5 OC partition does not close in ${config}" >&2
        exit 1
    }

    # HEMCO molecule inputs use target MW. Check emitted mass ratio.
    awk -v r="${fsoap}" -v c="${mw_corr}" \
        'BEGIN { x=r*c*150.0/28.01; if (x < 0.012999999 || x > 0.013000001) exit 1 }' || {
        echo "FAIL: FINNv2.5 CO-to-FSOAP mass ratio in ${config}" >&2
        exit 1
    }

    check_sensitivity_condition_balance "${config}"
}

check_config "${fullchem_config}"
check_config "${aerosol_config}"

for species in BCPI BCPO OCPI OCPO FSOAP_BRC_HS DBRCPOA_BRC_HS \
               NPBRCPOA_BRC_HS PBRCPOA_BRC_HS; do
    require "InvFINNv25_${species}" "${aerosol_diagn}"
done
require "InvFINNv25_SOAP" "${aerosol_diagn}"
require "FINNv25_CO_field_scaled_to_SOAP_mass_flux_existing_path" "${aerosol_diagn}"

# Both templates expose all-source BrC diagnostics plus optional injection
# examples.  Detailed surface FINN diagnostics are aerosol-template-specific.
for diagn in "${fullchem_diagn}" "${aerosol_diagn}"; do
    require "EmisFSOAP_Total" "${diagn}"
    require "EmisFOCPI_TotalColumn" "${diagn}"
    require "EmisFOCPO_TotalColumn" "${diagn}"
    require "InvFINNv25Inject_FSOAP" "${diagn}"
done

for species in FSOAS BRCSOA PWTC SWTC DBRCPOA FSOAP NPBRCPOA PBRCPOA FOCPI FOCPO; do
    require "      - ${species}" "${aerosol_geoschem}"
done

echo "PASS: FINNv2.5 BrC config contracts"
