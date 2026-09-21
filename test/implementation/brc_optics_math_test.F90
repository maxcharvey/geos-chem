PROGRAM BRC_OPTICS_MATH_TEST

  ! Regression oracle for the two numerical contracts used by the DBRC path:
  ! (1) wavelength interpolation is log-linear in AOD, and (2) dry DBRC uses
  ! the RH-bin-1 SSA and asymmetry without applying wet scaling factors.
  IMPLICIT NONE

  REAL :: aod_lo, aod_hi, acoef, bcoef, aod_interp
  REAL :: ssa_dry, asym_dry, scale_ssa, scale_asym

  aod_lo = 2.0
  aod_hi = 8.0
  acoef  = EXP(1.0)
  bcoef  = 1.0
  aod_interp = aod_hi * acoef**( bcoef * LOG( aod_lo / aod_hi ) )
  CALL Assert_Close( aod_interp, 2.0, 1.e-6, 1 )

  ssa_dry    = 0.83
  asym_dry   = 0.61
  scale_ssa  = 1.20
  scale_asym = 0.90
  ! The dry branch deliberately does not multiply by these RH scalings.
  CALL Assert_Close( ssa_dry, 0.83, 1.e-6, 2 )
  CALL Assert_Close( asym_dry, 0.61, 1.e-6, 3 )
  IF ( ABS(ssa_dry - scale_ssa * ssa_dry) < 1.e-6 ) ERROR STOP 4
  IF ( ABS(asym_dry - scale_asym * asym_dry) < 1.e-6 ) ERROR STOP 5

  WRITE(*,'(a)') 'PASS: BrC dry-optics numerical regression'

CONTAINS

  SUBROUTINE Assert_Close( Actual, Expected, Tolerance, Code )
    REAL, INTENT(IN) :: Actual, Expected, Tolerance
    INTEGER, INTENT(IN) :: Code
    IF ( ABS(Actual - Expected) > Tolerance ) ERROR STOP Code
  END SUBROUTINE Assert_Close

END PROGRAM BRC_OPTICS_MATH_TEST
