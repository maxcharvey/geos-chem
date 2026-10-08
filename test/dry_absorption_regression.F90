!------------------------------------------------------------------------------
! Default-off regression for exact-band dry absorption diagnostic arithmetic.
!------------------------------------------------------------------------------
PROGRAM Dry_Absorption_Regression

  USE Precision_Mod, ONLY : fp
  USE Dry_Absorption_Diag_Mod

  IMPLICIT NONE

  INTEGER  :: Failures, Status
  REAL(fp) :: Actual, Expected, Lower, Upper, Long_Abs, Short_Abs
  REAL(fp) :: Components(3)

  Failures = 0

  CALL Dry_Absorption_Optical_Depth( 0.2_fp, 0.75_fp, Actual, Status )
  CALL Require_Status( 'AAOD status', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'AAOD value', Actual, 0.05_fp, 1.0e-12_fp, Failures )

  CALL Dry_Absorption_Optical_Depth( 0.2_fp, 1.01_fp, Actual, Status )
  CALL Require_Status( 'invalid SSA', Status, DRYABS_INVALID_INPUT, Failures )

  CALL Layer_Absorption_Mm1( 2.0e-4_fp, 100.0_fp, Actual, Status )
  CALL Require_Status( 'layer coefficient status', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'layer coefficient', Actual, 2.0_fp, 1.0e-12_fp, Failures )

  CALL Interpolate_Dry_Absorption( 400.0_fp, 400.0_fp, 500.0_fp, &
                                   0.0_fp, 2.0_fp, Actual, Status )
  CALL Require_Status( 'exact zero node', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'exact zero value', Actual, 0.0_fp, 0.0_fp, Failures )

  Lower = 4.0_fp
  Upper = 4.0_fp * ( 500.0_fp / 400.0_fp )**(-2.0_fp)
  Expected = 4.0_fp * ( 450.0_fp / 400.0_fp )**(-2.0_fp)
  CALL Interpolate_Dry_Absorption( 450.0_fp, 400.0_fp, 500.0_fp, &
                                   Lower, Upper, Actual, Status )
  CALL Require_Status( 'power-law status', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'power-law value', Actual, Expected, &
                      2.0e-6_fp, Failures )

  CALL Interpolate_Dry_Absorption( 450.0_fp, 400.0_fp, 500.0_fp, &
                                   0.0_fp, 0.0_fp, Actual, Status )
  CALL Require_Status( 'two-zero bracket', Status, DRYABS_OK, Failures )

  CALL Interpolate_Dry_Absorption( 450.0_fp, 400.0_fp, 500.0_fp, &
                                   0.0_fp, 1.0_fp, Actual, Status )
  CALL Require_Status( 'mixed-zero bracket', Status, DRYABS_MIXED_ZERO, &
                       Failures )

  CALL Ambient_To_Reference_Mm1( 10.0_fp, 800.0_fp, 280.0_fp, &
                                  1000.0_fp, 250.0_fp, Actual, Status )
  CALL Require_Status( 'reference-volume status', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'reference-volume value', Actual, 14.0_fp, &
                      1.0e-12_fp, Failures )

  Long_Abs = 2.0_fp
  Short_Abs = Long_Abs * ( 405.0_fp / 664.0_fp )**(-3.0_fp)
  CALL Absorption_Angstrom_Exponent( Short_Abs, Long_Abs, 405.0_fp, &
                                     664.0_fp, Actual, Status )
  CALL Require_Status( 'AAE status', Status, DRYABS_OK, Failures )
  CALL Require_Close( 'AAE value', Actual, 3.0_fp, 2.0e-6_fp, Failures )

  Components = (/ 1.0_fp, 0.5_fp, 0.25_fp /)
  CALL Check_Component_Closure( 1.75_fp, Components, 1.0e-12_fp, &
                                1.0e-12_fp, Status )
  CALL Require_Status( 'component closure', Status, DRYABS_OK, Failures )
  CALL Check_Component_Closure( 2.0_fp, Components, 1.0e-12_fp, &
                                1.0e-12_fp, Status )
  CALL Require_Status( 'component nonclosure', Status, &
                       DRYABS_INVALID_INPUT, Failures )

  IF ( Failures /= 0 ) STOP 1
  WRITE( 6, '(a)' ) 'Dry absorption regression: PASS'

CONTAINS

  SUBROUTINE Require_Close( Label, Value, Reference, Tolerance, Failures )

    CHARACTER(LEN=*), INTENT(IN)    :: Label
    REAL(fp),         INTENT(IN)    :: Value, Reference, Tolerance
    INTEGER,          INTENT(INOUT) :: Failures

    IF ( ABS( Value - Reference ) > Tolerance ) THEN
       WRITE( 6, '(a,1x,a,2(1x,es14.6))' ) &
            'FAIL:', TRIM( Label ), Value, Reference
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Require_Close

  SUBROUTINE Require_Status( Label, Value, Reference, Failures )

    CHARACTER(LEN=*), INTENT(IN)    :: Label
    INTEGER,          INTENT(IN)    :: Value, Reference
    INTEGER,          INTENT(INOUT) :: Failures

    IF ( Value /= Reference ) THEN
       WRITE( 6, '(a,1x,a,2(1x,i0))' ) &
            'FAIL:', TRIM( Label ), Value, Reference
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Require_Status

END PROGRAM Dry_Absorption_Regression
