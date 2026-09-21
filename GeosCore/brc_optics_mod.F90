!------------------------------------------------------------------------------
! BrC requested-wavelength optical diagnostics.
!------------------------------------------------------------------------------
MODULE BRC_OPTICS_MOD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: BRC_DRY_AOD_AT_WAVELENGTH

CONTAINS

  PURE REAL(KIND=8) FUNCTION BRC_DRY_AOD_AT_WAVELENGTH( AOD_LO, AOD_HI, &
                                                         LINTERP, ACOEF, BCOEF )

    REAL(KIND=8), INTENT(IN) :: AOD_LO, AOD_HI, ACOEF, BCOEF
    LOGICAL,      INTENT(IN) :: LINTERP

    ! A zero endpoint has no defined Angstrom exponent.  Returning zero also
    ! clears an archived diagnostic rather than retaining a previous timestep.
    IF ( .NOT. LINTERP ) THEN
       BRC_DRY_AOD_AT_WAVELENGTH = AOD_LO
    ELSE IF ( AOD_LO > 0.0D0 .AND. AOD_HI > 0.0D0 ) THEN
       BRC_DRY_AOD_AT_WAVELENGTH = AOD_HI * ACOEF**( BCOEF * LOG(AOD_LO / AOD_HI) )
    ELSE
       BRC_DRY_AOD_AT_WAVELENGTH = 0.0D0
    ENDIF

  END FUNCTION BRC_DRY_AOD_AT_WAVELENGTH

END MODULE BRC_OPTICS_MOD
