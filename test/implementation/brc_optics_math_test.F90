PROGRAM BRC_OPTICS_MATH_TEST

  USE BRC_OPTICS_MOD, ONLY : BRC_DRY_AOD_AT_WAVELENGTH

  ! Exercise the exact production helper used by the DBRC diagnostic.
  IMPLICIT NONE

  REAL(KIND=8) :: aod

  ! Exact requested wavelength (e.g., 550 nm).
  aod = BRC_DRY_AOD_AT_WAVELENGTH( 3.D-1, 9.D-1, .FALSE., 1.D0, 1.D0 )
  CALL Assert_Close( aod, 3.D-1, 1.D-12, 1 )

  ! Interpolated requested wavelengths 527.1 and 693.5 nm use their
  ! precomputed ACOEF/BCOEF pair; test two distinct nontrivial pairs.
  aod = BRC_DRY_AOD_AT_WAVELENGTH( 2.D-1, 8.D-1, .TRUE., EXP(1.D0), 5.D-1 )
  CALL Assert_Close( aod, 4.D-1, 1.D-12, 2 )
  aod = BRC_DRY_AOD_AT_WAVELENGTH( 8.D-1, 2.D-1, .TRUE., EXP(1.D0), 5.D-1 )
  CALL Assert_Close( aod, 4.D-1, 1.D-12, 3 )

  ! A zero endpoint must clear the diagnostic, not retain stale data.
  aod = BRC_DRY_AOD_AT_WAVELENGTH( 0.D0, 8.D-1, .TRUE., EXP(1.D0), 5.D-1 )
  CALL Assert_Close( aod, 0.D0, 1.D-12, 4 )
  aod = BRC_DRY_AOD_AT_WAVELENGTH( 8.D-1, 0.D0, .TRUE., EXP(1.D0), 5.D-1 )
  CALL Assert_Close( aod, 0.D0, 1.D-12, 5 )

  WRITE(*,'(a)') 'PASS: BrC requested-wavelength regression'

CONTAINS

  SUBROUTINE Assert_Close( Actual, Expected, Tolerance, Code )
    REAL(KIND=8), INTENT(IN) :: Actual, Expected, Tolerance
    INTEGER, INTENT(IN) :: Code
    IF ( ABS(Actual - Expected) > Tolerance ) ERROR STOP Code
  END SUBROUTINE Assert_Close

END PROGRAM BRC_OPTICS_MATH_TEST
