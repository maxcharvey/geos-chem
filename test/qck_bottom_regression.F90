!------------------------------------------------------------------------------
! Default-off compiled regression matrix for Qck positivity corrections.
!------------------------------------------------------------------------------
PROGRAM Qck_Bottom_Regression

  USE Precision_Mod, ONLY : fp, f8
  USE Qck_Positivity_Mod, ONLY : Qck_Bottom_Conservative, Qck_Interior_Correct, &
                                 QCK_BOTTOM_EXACT, QCK_BOTTOM_ROUNDOFF,       &
                                 QCK_BOTTOM_UNFILLABLE, QCK_BOTTOM_INVALID_DONOR, &
                                 QCK_BOTTOM_MICROCLOSURE, Canonical_Diagnostic_Time
  USE Qck_Bottom_Sizing_Mod, ONLY : Qck_Bottom_Sizing_Candidate_Is_Better

  IMPLICIT NONE

  INTEGER :: Failures

  Failures = 0
  CALL Test_Interior_Conservation( Failures )
  CALL Test_Immediate_Donor( Failures )
  CALL Test_Farther_Column_Donor( Failures )
  CALL Test_Roundoff_Closure( Failures )
  CALL Test_Unfillable_Column( Failures )
  CALL Test_Microclosure_Enabled( Failures )
  CALL Test_Microclosure_Over_Limit( Failures )
  CALL Test_Day5_Isop_Ceiling_Boundary( Failures )
  CALL Test_Negative_Donor_Precondition( Failures )
  CALL Test_Canonical_Diagnostic_Time( Failures )
  CALL Test_Sizing_Candidate_Order( Failures )

  IF ( Failures /= 0 ) STOP 1
  WRITE( 6, '(a)' ) 'Qck bottom regression matrix: PASS'

CONTAINS

  SUBROUTINE Test_Sizing_Candidate_Order( Failures )

    INTEGER, INTENT(INOUT) :: Failures

    CALL Require_Logical( 'sizing larger value', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 2, 159, 72, 60, 39, &
              1.0_f8, 20190805, 190000, 1, 1, 72, 1, 1 ), &
         .TRUE., Failures )
    CALL Require_Logical( 'sizing earlier ordinal tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 2, 159, 72, 60, 39, &
              2.0_f8, 20190805, 190000, 3, 1, 1, 1, 1 ), &
         .TRUE., Failures )
    CALL Require_Logical( 'sizing later ordinal loses tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 4, 1, 1, 1, 1, &
              2.0_f8, 20190805, 190000, 3, 159, 72, 60, 39 ), &
         .FALSE., Failures )
    CALL Require_Logical( 'sizing lower species tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 3, 158, 72, 60, 39, &
              2.0_f8, 20190805, 190000, 3, 159, 1, 1, 1 ), &
         .TRUE., Failures )
    CALL Require_Logical( 'sizing lower K tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 3, 159, 71, 60, 39, &
              2.0_f8, 20190805, 190000, 3, 159, 72, 1, 1 ), &
         .TRUE., Failures )
    CALL Require_Logical( 'sizing lower J tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 3, 159, 72, 59, 39, &
              2.0_f8, 20190805, 190000, 3, 159, 72, 60, 1 ), &
         .TRUE., Failures )
    CALL Require_Logical( 'sizing lower I tie', &
         Qck_Bottom_Sizing_Candidate_Is_Better( &
              2.0_f8, 20190805, 190000, 3, 159, 72, 60, 38, &
              2.0_f8, 20190805, 190000, 3, 159, 72, 60, 39 ), &
         .TRUE., Failures )

  END SUBROUTINE Test_Sizing_Candidate_Order

  SUBROUTINE Test_Canonical_Diagnostic_Time( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: This_Date, This_Time

    CALL Canonical_Diagnostic_Time( 20190101, 0, 150, This_Date, This_Time )
    CALL Require_Integer( 'canonical 150-second date', This_Date, 20190101, &
                          Failures )
    CALL Require_Integer( 'canonical 150-second time', This_Time, 230, &
                          Failures )

    CALL Canonical_Diagnostic_Time( 20191231, 235930, 60, This_Date, This_Time )
    CALL Require_Integer( 'canonical year-rollover date', This_Date, 20200101, &
                          Failures )
    CALL Require_Integer( 'canonical year-rollover time', This_Time, 30, &
                          Failures )

    CALL Canonical_Diagnostic_Time( 20200228, 235930, 60, This_Date, This_Time )
    CALL Require_Integer( 'canonical leap-day date', This_Date, 20200229, &
                          Failures )
    CALL Require_Integer( 'canonical leap-day time', This_Time, 30, Failures )

  END SUBROUTINE Test_Canonical_Diagnostic_Time

  SUBROUTINE Test_Interior_Conservation( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    REAL(fp) :: Above, Target, Below, Initial_Total

    Above = 7.0_fp
    Target = -10.0_fp
    Below = 4.0_fp
    Initial_Total = Above + Target + Below
    CALL Qck_Interior_Correct( Above, Target, Below )
    CALL Require_Close( 'interior above', Above, 0.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'interior target', Target, 0.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'interior below', Below, 1.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'interior mass', Above + Target + Below, Initial_Total, &
                        1.0e-12_fp, Failures )

  END SUBROUTINE Test_Interior_Conservation

  SUBROUTINE Test_Immediate_Donor( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Deficit, Available, Withdrawn, Closure, Tolerance
    REAL(fp) :: Initial_Total

    Column = (/ 10.0_fp, -7.0_fp /)
    Initial_Total = SUM( Column )
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count )
    CALL Require_Integer( 'immediate status', Status, QCK_BOTTOM_EXACT, Failures )
    CALL Require_Close( 'immediate donor', Column(1), 3.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'immediate bottom', Column(2), 0.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'immediate mass', SUM( Column ), Initial_Total, &
                        1.0e-12_fp, Failures )
    CALL Require_Integer( 'immediate donor count', Donor_Count, 1, Failures )
    CALL Require_Close( 'immediate closure', Closure, 0.0_fp, 1.0e-12_fp, Failures )

  END SUBROUTINE Test_Immediate_Donor

  SUBROUTINE Test_Farther_Column_Donor( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(3), Deficit, Available, Withdrawn, Closure, Tolerance
    REAL(fp) :: Initial_Total

    Column = (/ 8.0_fp, 7.0_fp, -10.0_fp /)
    Initial_Total = SUM( Column )
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count )
    CALL Require_Integer( 'farther status', Status, QCK_BOTTOM_EXACT, Failures )
    CALL Require_Close( 'farther upper donor', Column(1), 5.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'farther immediate donor', Column(2), 0.0_fp, &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'farther bottom', Column(3), 0.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'farther mass', SUM( Column ), Initial_Total, &
                        1.0e-12_fp, Failures )
    CALL Require_Integer( 'farther donor count', Donor_Count, 2, Failures )

  END SUBROUTINE Test_Farther_Column_Donor

  SUBROUTINE Test_Roundoff_Closure( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Deficit, Available, Withdrawn, Closure, Tolerance

    ! Use a test-local tolerance that is resolvable for both fp=f4 and fp=f8.
    Column = (/ 1.0_fp - 1.0e-5_fp, -1.0_fp /)
    CALL Qck_Bottom_Conservative( Column, 1.0e-4_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count )
    CALL Require_Integer( 'roundoff status', Status, QCK_BOTTOM_ROUNDOFF, Failures )
    CALL Require_Close( 'roundoff donor', Column(1), 0.0_fp, 1.0e-12_fp, Failures )
    CALL Require_Close( 'roundoff bottom', Column(2), 0.0_fp, 1.0e-12_fp, Failures )
    IF ( Closure <= 0.0_fp .or. Closure > Tolerance ) THEN
       WRITE( 6, '(a,2(1x,es14.6))' ) 'FAIL: roundoff closure/tolerance', &
            Closure, Tolerance
       Failures = Failures + 1
    ENDIF
    CALL Require_Integer( 'roundoff donor count', Donor_Count, 1, Failures )

  END SUBROUTINE Test_Roundoff_Closure

  SUBROUTINE Test_Unfillable_Column( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Initial(2), Deficit, Available, Withdrawn, Closure
    REAL(fp) :: Tolerance

    Column = (/ 7.0_fp, -10.0_fp /)
    Initial = Column
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count )
    CALL Require_Integer( 'unfillable status', Status, QCK_BOTTOM_UNFILLABLE, Failures )
    CALL Require_Close( 'unfillable donor unchanged', Column(1), Initial(1), &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'unfillable bottom unchanged', Column(2), Initial(2), &
                        1.0e-12_fp, Failures )
    CALL Require_Integer( 'unfillable donor count', Donor_Count, 0, Failures )

  END SUBROUTINE Test_Unfillable_Column

  SUBROUTINE Test_Microclosure_Enabled( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Deficit, Available, Withdrawn, Closure, Tolerance

    Column = (/ 7.0_fp, -10.0_fp /)
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count, &
         Microclosure_Tolerance=4.0_fp )
    CALL Require_Integer( 'microclosure status', Status, &
                          QCK_BOTTOM_MICROCLOSURE, Failures )
    CALL Require_Close( 'microclosure donor', Column(1), 0.0_fp, &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'microclosure bottom', Column(2), 0.0_fp, &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'microclosure withdrawn', Withdrawn, 7.0_fp, &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'microclosure declared closure', Closure, 3.0_fp, &
                        1.0e-12_fp, Failures )
    CALL Require_Integer( 'microclosure donor count', Donor_Count, 1, Failures )
    IF ( Closure <= Tolerance ) THEN
       WRITE( 6, '(a,2(1x,es14.6))' ) &
            'FAIL: microclosure was ordinary roundoff', Closure, Tolerance
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Test_Microclosure_Enabled

  SUBROUTINE Test_Microclosure_Over_Limit( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Initial(2), Deficit, Available, Withdrawn, Closure
    REAL(fp) :: Tolerance

    Column = (/ 7.0_fp, -10.0_fp /)
    Initial = Column
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count, &
         Microclosure_Tolerance=2.0_fp )
    CALL Require_Integer( 'microclosure over-limit status', Status, &
                          QCK_BOTTOM_UNFILLABLE, Failures )
    CALL Require_Close( 'microclosure over-limit donor unchanged', Column(1), &
                        Initial(1), 1.0e-12_fp, Failures )
    CALL Require_Close( 'microclosure over-limit bottom unchanged', Column(2), &
                        Initial(2), 1.0e-12_fp, Failures )
    CALL Require_Integer( 'microclosure over-limit donor count', Donor_Count, &
                          0, Failures )

  END SUBROUTINE Test_Microclosure_Over_Limit

  SUBROUTINE Test_Day5_Isop_Ceiling_Boundary( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Initial(2), Deficit, Available, Withdrawn, Closure
    REAL(fp) :: Tolerance
    REAL(fp), PARAMETER :: Day5_Available_Hpa = 3.12163e-14_fp
    REAL(fp), PARAMETER :: Day5_Deficit_Hpa = 3.98814e-13_fp
    REAL(fp), PARAMETER :: Event_015kg_Hpa = 2.6949961043547083e-13_fp
    REAL(fp), PARAMETER :: Event_025kg_Hpa = 4.4916601739245140e-13_fp

    ! Reproduce the ISOP QCK_BOTTOM event that stopped the 2019-08-05 run at
    ! 18:50 UTC.  The 0.15 kg ceiling must fail closed without mutation.
    Column = (/ Day5_Available_Hpa, -Day5_Deficit_Hpa /)
    Initial = Column
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count, &
         Microclosure_Tolerance=Event_015kg_Hpa )
    CALL Require_Integer( 'day5 ISOP 0.15kg status', Status, &
                          QCK_BOTTOM_UNFILLABLE, Failures )
    CALL Require_Close( 'day5 ISOP 0.15kg donor unchanged', Column(1), &
                        Initial(1), 1.0e-20_fp, Failures )
    CALL Require_Close( 'day5 ISOP 0.15kg bottom unchanged', Column(2), &
                        Initial(2), 1.0e-20_fp, Failures )

    ! The reviewed 0.25 kg ceiling admits exactly the same measured event as
    ! a declared numerical microclosure.
    Column = Initial
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count, &
         Microclosure_Tolerance=Event_025kg_Hpa )
    CALL Require_Integer( 'day5 ISOP 0.25kg status', Status, &
                          QCK_BOTTOM_MICROCLOSURE, Failures )
    CALL Require_Close( 'day5 ISOP 0.25kg donor', Column(1), 0.0_fp, &
                        1.0e-20_fp, Failures )
    CALL Require_Close( 'day5 ISOP 0.25kg bottom', Column(2), 0.0_fp, &
                        1.0e-20_fp, Failures )
    CALL Require_Close( 'day5 ISOP measured closure', Closure, &
                        Day5_Deficit_Hpa - Day5_Available_Hpa, &
                        1.0e-20_fp, Failures )

  END SUBROUTINE Test_Day5_Isop_Ceiling_Boundary

  SUBROUTINE Test_Negative_Donor_Precondition( Failures )

    INTEGER, INTENT(INOUT) :: Failures
    INTEGER :: Status, Donor_Count
    REAL(fp) :: Column(2), Initial(2), Deficit, Available, Withdrawn, Closure
    REAL(fp) :: Tolerance

    Column = (/ -1.0_fp, -10.0_fp /)
    Initial = Column
    CALL Qck_Bottom_Conservative( Column, 1.0e-9_fp, 1.0e-30_fp, Status, &
         Deficit, Available, Withdrawn, Closure, Tolerance, Donor_Count )
    CALL Require_Integer( 'negative donor status', Status, &
                          QCK_BOTTOM_INVALID_DONOR, Failures )
    CALL Require_Close( 'negative donor unchanged', Column(1), Initial(1), &
                        1.0e-12_fp, Failures )
    CALL Require_Close( 'negative bottom unchanged', Column(2), Initial(2), &
                        1.0e-12_fp, Failures )

  END SUBROUTINE Test_Negative_Donor_Precondition

  SUBROUTINE Require_Close( Label, Actual, Expected, Tolerance, Failures )

    CHARACTER(LEN=*), INTENT(IN) :: Label
    REAL(fp),         INTENT(IN) :: Actual, Expected, Tolerance
    INTEGER,          INTENT(INOUT) :: Failures

    IF ( ABS( Actual - Expected ) > Tolerance ) THEN
       WRITE( 6, '(a,1x,a,2(1x,es14.6))' ) 'FAIL:', TRIM( Label ), Actual, Expected
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Require_Close

  SUBROUTINE Require_Integer( Label, Actual, Expected, Failures )

    CHARACTER(LEN=*), INTENT(IN) :: Label
    INTEGER,          INTENT(IN) :: Actual, Expected
    INTEGER,          INTENT(INOUT) :: Failures

    IF ( Actual /= Expected ) THEN
       WRITE( 6, '(a,1x,a,2(1x,i0))' ) 'FAIL:', TRIM( Label ), Actual, Expected
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Require_Integer

  SUBROUTINE Require_Logical( Label, Actual, Expected, Failures )

    CHARACTER(LEN=*), INTENT(IN) :: Label
    LOGICAL,          INTENT(IN) :: Actual, Expected
    INTEGER,          INTENT(INOUT) :: Failures

    IF ( Actual .neqv. Expected ) THEN
       WRITE( 6, '(a,1x,a)' ) 'FAIL:', TRIM( Label )
       Failures = Failures + 1
    ENDIF

  END SUBROUTINE Require_Logical

END PROGRAM Qck_Bottom_Regression
