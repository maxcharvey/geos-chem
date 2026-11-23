#!/usr/bin/env python3
"""Extract the production bleaching routines for a standalone Fortran test."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "GeosCore/brc_mod.F90"
DIAGNOSTICS = (
    "BrCTauBleach", "BrCKBleach", "BrCEtaBBOA", "BrCTemp", "BrCRH",
    "BrCO3ppbv", "BrCFluxBRC2WTC", "BrCFluxNPBRC2WTC",
)

# Keep the model-facing test state small. Chemistry storage is REAL8 and
# diagnostic storage is REAL4, as in a full-precision Classic build.
INTERFACES = """module actual_bleaching
use iso_fortran_env, only: fp=>real64, f4=>real32
implicit none
integer, parameter :: nx=2, ny=1, nz=3, GC_SUCCESS=0
integer :: ozone_index=3
type OptInput
   integer :: BrC_Bleach_Scheme
end type
type SpcConc
   real(fp), pointer :: Conc(:,:,:) => null()
end type
type ChmState
   type(SpcConc) :: Species(3)
end type
type GrdState
   integer :: NX=2, NY=1, NZ=3
end type
type MetState
   real(fp) :: T(nx,ny,nz), RH(nx,ny,nz), BXHEIGHT(nx,ny,nz)
   real(fp) :: AD(nx,ny,nz), PMID(nx,ny,nz)
end type
type DgnState
"""
DIAGNOSTIC_MEMBERS = "\n".join(
    f"   logical :: Archive_{name}=.true.\n   real(f4) :: {name}(nx,ny,nz)"
    for name in DIAGNOSTICS
)
STATE_AND_HELPERS = """
end type
real(fp) :: FSOAS_CONV(nx,ny,nz), BRCSOA_CONV(nx,ny,nz), NPBRC_CONV(nx,ny,nz)
"""
HELPERS = """
contains
integer function Ind_(name)
   character(*), intent(in) :: name
   Ind_=0
   if (name == 'O3') Ind_=ozone_index
end function
integer function GET_TS_CHEM()
   GET_TS_CHEM=1200
end function
"""


def extract_routine(source, kind, name):
    routine = re.search(
        rf"^ {kind} {name}\(.*?^ END {kind} {name}$", source, re.M | re.S,
    ).group()
    # Replace only the model-module imports with the small interfaces above.
    return re.sub(r"^ +USE .*\n", "", routine, flags=re.M)


def write_fixture(target):
    source = SOURCE.read_text()
    declarations, _ = source.split("\nCONTAINS\n", 1)
    constants = re.findall(r"^  REAL\(fp\), PARAMETER :: .*$", declarations, re.M)
    routines = [
        extract_routine(source, kind, name)
        for kind, name in (
            ("FUNCTION", "VISC_WATER_FUNC"), ("FUNCTION", "VISC_BBOA_FUNC"),
            ("FUNCTION", "CALC_TAU_BRC"), ("SUBROUTINE", "CHEM_BRCSOA"),
            ("SUBROUTINE", "CHEM_NPBRCPOA"),
        )
    ]
    target.write_text(
        INTERFACES + DIAGNOSTIC_MEMBERS + STATE_AND_HELPERS
        + "\n".join(constants) + HELPERS + "\n".join(routines) + "\nend module\n"
    )


if __name__ == "__main__":
    write_fixture(Path(sys.argv[1]))
