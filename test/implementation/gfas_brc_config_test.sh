#!/usr/bin/env bash
# Static checks for the opt-in GFAS BrC harmonized sensitivity.

set -euo pipefail

this_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
geos_root=$(git -C "${this_dir}" rev-parse --show-toplevel)
config_dir="${geos_root}/run/GCClassic/HEMCO_Config.rc.templates"
diagn_dir="${geos_root}/run/GCClassic/HEMCO_Diagn.rc.templates"
configs=("${config_dir}/HEMCO_Config.rc.fullchem" \
         "${config_dir}/HEMCO_Config.rc.aerosol")
fullchem_diagn="${diagn_dir}/HEMCO_Diagn.rc.fullchem"
aerosol_diagn="${diagn_dir}/HEMCO_Diagn.rc.aerosol"
aerosol_geos_config="${geos_root}/run/GCClassic/geoschem_config.yml.templates/geoschem_config.yml.aerosol"

require() {
    local pattern="${1}"
    local file="${2}"
    if ! rg -q --fixed-strings "${pattern}" "${file}"; then
        echo "FAIL: Missing '${pattern}' in ${file}" >&2
        exit 1
    fi
}

require_diag() {
    local name="${1}"
    local spec="${2}"
    local dim="${3}"
    local file="${4}"
    if ! awk -v name="${name}" -v spec="${spec}" -v dim="${dim}" \
        '$1 == name && $2 == spec && $3 == 0 && $4 == 5 && $5 == 3 && $6 == dim { found = 1 }
         END { exit !found }' "${file}"; then
        echo "FAIL: Missing diagnostic selector '${name}' in ${file}" >&2
        exit 1
    fi
}

# Check both simulation templates.
for config in "${configs[@]}"; do
    # Default-off gate and balanced explicit off/on OC paths.
    require "GFAS_BRC_HARMONIZED_SENSITIVITY : false" "${config}"
    for marker in "(((.not.GFAS_BRC_HARMONIZED_SENSITIVITY" \
                  "))).not.GFAS_BRC_HARMONIZED_SENSITIVITY" \
                  "(((GFAS_BRC_HARMONIZED_SENSITIVITY" \
                  ")))GFAS_BRC_HARMONIZED_SENSITIVITY"; do
        if [[ $(rg -c --fixed-strings "${marker}" "${config}") -ne 1 ]]; then
            echo "FAIL: Unbalanced or duplicate '${marker}' in ${config}" >&2
            exit 1
        fi
    done
    require "GFAS_OCPI  \$ROOT/GFAS/v2018-09/\$YYYY/GFAS_\$YYYY\$MM.nc ocfire" "${config}"
    require "GFAS_OCPO  \$ROOT/GFAS/v2018-09/\$YYYY/GFAS_\$YYYY\$MM.nc ocfire" "${config}"

    # Controlled proxy pathways and runtime scale factors.
    for source in GFAS_FSOAP_HARM GFAS_DBRCPOA_HARM GFAS_OCPI_HARM \
                  GFAS_OCPO_HARM GFAS_NPBRCPOA_HARM GFAS_PBRCPOA_HARM; do
        require "${source}" "${config}"
    done
    require "282 GFAS_BRC_OC_RESIDUAL  0.500" "${config}"
    require "283 GFAS_BRC_NPBRC_FROM_OC 0.375" "${config}"
    require "284 GFAS_BRC_PBRC_FROM_OC  0.125" "${config}"
    require "285 GFAS_BRC_DBRC_FROM_BC   4.000" "${config}"
    require "286 GFAS_BRC_FSOAP_FROM_CO  0.013" "${config}"

    # Every scale ID is defined once and used by the intended pathways only.
    for id in 282 283 284 285 286; do
        expected=1
        [[ "${id}" == 282 ]] && expected=2
        awk -v id="${id}" -v expected="${expected}" '
            $1 == id { definitions++ }
            $1 == 0 {
                for (field = 1; field <= NF; field++) {
                    count = split($field, values, "/")
                    for (value = 1; value <= count; value++) {
                        if (values[value] == id) uses++
                    }
                }
            }
            END { exit !(definitions == 1 && uses == expected) }
        ' "${config}" || {
            echo "FAIL: Scale ID ${id} definition/use count in ${config}" >&2
            exit 1
        }
    done

    # OC mass closure: residual + NPBRCPOA + PBRCPOA = native OC.
    awk '
        $2 == "GFAS_BRC_OC_RESIDUAL"  { residual = $3 }
        $2 == "GFAS_BRC_NPBRC_FROM_OC" { npbrc = $3 }
        $2 == "GFAS_BRC_PBRC_FROM_OC"  { pbrc = $3 }
        END {
            total = residual + npbrc + pbrc
            if (total < 0.999999 || total > 1.000001) exit 1
        }
    ' "${config}" || {
        echo "FAIL: GFAS OC closure in ${config}" >&2
        exit 1
    }
done

# Inventory-labelled diagnostics expose native and proxy-added fluxes with
# exact ExtNr=0, Cat=5, Hier=3 selectors in both dimensions.
common_diags=("EmisGFAS_BCPI_Native:BCPI" "EmisGFAS_BCPO_Native:BCPO" \
              "EmisGFAS_OCPI:OCPI" "EmisGFAS_OCPO:OCPO" \
              "EmisGFAS_FSOAP_HarmSens:FSOAP" \
              "EmisGFAS_DBRCPOA_HarmSens:DBRCPOA" \
              "EmisGFAS_NPBRCPOA_HarmSens:NPBRCPOA" \
              "EmisGFAS_PBRCPOA_HarmSens:PBRCPOA")
for diagn in "${fullchem_diagn}" "${aerosol_diagn}"; do
    for entry in "${common_diags[@]}"; do
        name=${entry%%:*}
        spec=${entry#*:}
        require_diag "${name}" "${spec}" 2 "${diagn}"
        require_diag "${name}3D" "${spec}" 3 "${diagn}"
    done
done
require_diag EmisGFAS_CO_Native CO 2 "${fullchem_diagn}"
require_diag EmisGFAS_CO_Native3D CO 3 "${fullchem_diagn}"
require_diag EmisGFAS_SOAP SOAP 2 "${aerosol_diagn}"
require_diag EmisGFAS_SOAP3D SOAP 3 "${aerosol_diagn}"

# Aerosol-only must transport the complete BrC chain emitted by its template.
for species in FSOAS BRCSOA WTC DBRCPOA FSOAP NPBRCPOA PBRCPOA; do
    if ! awk -v species="${species}" '
        /transported_species:/ { in_transport = 1; next }
        in_transport && /^[^ ]/ { in_transport = 0 }
        in_transport && $1 == "-" && $2 == species { found = 1 }
        END { exit !found }
    ' "${aerosol_geos_config}"; then
        echo "FAIL: Missing transported aerosol BrC species '${species}'" >&2
        exit 1
    fi
done

echo "PASS: GFAS BrC harmonized-sensitivity configuration checks"
