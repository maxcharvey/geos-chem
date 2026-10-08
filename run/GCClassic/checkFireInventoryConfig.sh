#!/usr/bin/env bash
# Fail fast on inconsistent GFAS/FINNv2.5 BrC sensitivity switches.

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 geoschem_config.yml HEMCO_Config.rc" >&2
  exit 2
fi

GC_CONFIG=$1
HEMCO_CONFIG=$2

[[ -f $GC_CONFIG ]] || { echo "ERROR: missing $GC_CONFIG" >&2; exit 2; }
[[ -f $HEMCO_CONFIG ]] || { echo "ERROR: missing $HEMCO_CONFIG" >&2; exit 2; }

yaml_true() {
  awk -F: -v key="$2" '
    $1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
      value=tolower($2); sub(/^[[:space:]]*/, "", value)
      sub(/[[:space:]#].*$/, "", value)
      found=(value == "true")
    }
    END { exit(found ? 0 : 1) }
  ' "$1"
}

hemco_true() {
  awk -v key="$2" '
    $0 ~ "-->[[:space:]]*" key "[[:space:]]*:" {
      value=tolower($0); sub(/^.*:[[:space:]]*/, "", value)
      sub(/^[[:space:]]*/, "", value); sub(/[[:space:]#].*$/, "", value)
      found=(value == "true")
    }
    END { exit(found ? 0 : 1) }
  ' "$1"
}

brc=false; gfas=false; finn=false; gfas_brc=false; finn_brc=false
yaml_true "$GC_CONFIG" brown_carbon && brc=true
hemco_true "$HEMCO_CONFIG" GFAS && gfas=true
hemco_true "$HEMCO_CONFIG" FINNv25 && finn=true
hemco_true "$HEMCO_CONFIG" GFAS_BRC_HARMONIZED_SENSITIVITY && gfas_brc=true
hemco_true "$HEMCO_CONFIG" FINNV25_BRC_HARMONIZED_SENSITIVITY && finn_brc=true

fail() { echo "ERROR: $*" >&2; exit 1; }

[[ $gfas == false || $finn == false ]] || fail "GFAS and FINNv25 both enabled"
[[ $gfas_brc == false || $finn_brc == false ]] || fail "both BrC sensitivities enabled"
[[ $gfas_brc == false || $gfas == true ]] || fail "GFAS BrC sensitivity requires GFAS"
[[ $finn_brc == false || $finn == true ]] || fail "FINNv2.5 BrC sensitivity requires FINNv25"
[[ $gfas_brc == false && $finn_brc == false || $brc == true ]] || \
  fail "BrC sensitivity requires brown_carbon: true"

echo "PASS: brown_carbon=$brc GFAS=$gfas FINNv25=$finn GFAS_BrC=$gfas_brc FINNv25_BrC=$finn_brc"
