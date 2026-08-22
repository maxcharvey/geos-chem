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

entry_scale() {
    local file="${1}"
    local name="${2}"
    awk -v name="${name}" '$2 == name { print $(NF-2); exit }' "${file}"
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
        in_finn && /^\(\(\(\.not\.FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            if (off_open || on_open || off_close || on_close) bad=1
            off_open++
        }
        in_finn && /^\)\)\)\.not\.FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            if (off_open != 1 || off_close || on_open || on_close) bad=1
            off_close++
        }
        in_finn && /^\(\(\(FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            if (off_close != 1 || on_open || on_close) bad=1
            on_open++
        }
        in_finn && /^\)\)\)FINNV25_BRC_HARMONIZED_SENSITIVITY$/ {
            if (on_open != 1 || on_close) bad=1
            on_close++
        }
        in_finn && /^\)\)\)FINNv25$/ {
            if (off_open != 1 || off_close != 1 || on_open != 1 || on_close != 1) bad=1
            in_finn=0
        }
        END {
            if (in_finn || off_open != 1 || off_close != 1 ||
                on_open != 1 || on_close != 1 || bad) exit 1
        }
    ' "${file}" || {
        echo "FAIL: unbalanced FINNv2.5 sensitivity conditions in ${file}" >&2
        exit 1
    }
}

check_config() {
    local config="${1}"

    # Opt-in only; ordinary FINNv2.5 OC mapping remains explicit.
    require "--> FINNV25_BRC_HARMONIZED_SENSITIVITY : false" "${config}"
    require "(((.not.FINNV25_BRC_HARMONIZED_SENSITIVITY" "${config}"
    [[ $(entry_scale "${config}" FINNv25_OCPI) == "75/72" ]]
    [[ $(entry_scale "${config}" FINNv25_OCPO) == "75/73" ]]

    # Controlled mappings and molecule-to-mass correction are explicit.
    for name in FSOAP DBRCPOA NPBRCPOA PBRCPOA; do
        require "0 FINNv25_${name}_BRC_HS " "${config}"
    done
    [[ $(entry_scale "${config}" FINNv25_OCPI_BRC_HS) == "75/282/72" ]]
    [[ $(entry_scale "${config}" FINNv25_OCPO_BRC_HS) == "75/282/73" ]]

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

for diagn in "${fullchem_diagn}" "${aerosol_diagn}"; do
    for species in BCPI BCPO OCPI OCPO FSOAP_BRC_HS DBRCPOA_BRC_HS \
                   NPBRCPOA_BRC_HS PBRCPOA_BRC_HS; do
        require "InvFINNv25_${species}" "${diagn}"
    done
done
require "InvFINNv25_CO" "${fullchem_diagn}"
require "InvFINNv25_SOAP" "${aerosol_diagn}"
require "FINNv25_CO_field_scaled_to_SOAP_mass_flux_existing_path" "${aerosol_diagn}"

for species in FSOAS BRCSOA WTC DBRCPOA FSOAP NPBRCPOA PBRCPOA; do
    require "      - ${species}" "${aerosol_geoschem}"
done

echo "PASS: FINNv2.5 BrC config contracts"
