!------------------------------------------------------------------------------
! Shared positivity corrections used by the TPCORE Qck routines.
!
! Qck interior behavior is kept byte-for-byte algebraically equivalent to the
! native stencil.  Qck bottom behavior is deliberately stricter: it mutates a
! column only after a nonnegative conservative correction (or an explicitly
! bounded roundoff closure) has been preflighted.
!------------------------------------------------------------------------------
MODULE Qck_Positivity_Mod

  USE Precision_Mod, ONLY : fp

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: QCK_BOTTOM_EXACT          = 0
  INTEGER, PARAMETER, PUBLIC :: QCK_BOTTOM_ROUNDOFF       = 1
  INTEGER, PARAMETER, PUBLIC :: QCK_BOTTOM_UNFILLABLE     = 2
  INTEGER, PARAMETER, PUBLIC :: QCK_BOTTOM_INVALID_DONOR  = 3
  INTEGER, PARAMETER, PUBLIC :: QCK_BOTTOM_MICROCLOSURE   = 4

  REAL(fp), PARAMETER, PUBLIC :: QCK_BOTTOM_RELATIVE_TOLERANCE = 1.0e-9_fp
  REAL(fp), PARAMETER, PUBLIC :: QCK_BOTTOM_ABSOLUTE_TOLERANCE = 1.0e-30_fp

  PUBLIC :: Qck_Interior_Correct
  PUBLIC :: Qck_Bottom_Conservative
  PUBLIC :: Canonical_Diagnostic_Time

CONTAINS

  SUBROUTINE Canonical_Diagnostic_Time( Start_Date, Start_Time, &
                                        Elapsed_Sec, This_Date, This_Time )

    INTEGER, INTENT(IN)  :: Start_Date
    INTEGER, INTENT(IN)  :: Start_Time
    INTEGER, INTENT(IN)  :: Elapsed_Sec
    INTEGER, INTENT(OUT) :: This_Date
    INTEGER, INTENT(OUT) :: This_Time

    INTEGER(KIND=8) :: Day_Offset, Seconds_Of_Day, Total_Seconds
    INTEGER         :: Day, Days_This_Month, Hour, Minute, Month, Second, Year

    Year   = Start_Date / 10000
    Month  = MOD( Start_Date / 100, 100 )
    Day    = MOD( Start_Date, 100 )
    Hour   = Start_Time / 10000
    Minute = MOD( Start_Time / 100, 100 )
    Second = MOD( Start_Time, 100 )

    Total_Seconds = INT( Hour * 3600 + Minute * 60 + Second, KIND=8 ) + &
                    INT( Elapsed_Sec, KIND=8 )
    Day_Offset     = Total_Seconds / 86400_8
    Seconds_Of_Day = MOD( Total_Seconds, 86400_8 )

    DO WHILE ( Day_Offset > 0_8 )
       Days_This_Month = Days_In_Month( Year, Month )
       IF ( Day_Offset <= INT( Days_This_Month - Day, KIND=8 ) ) THEN
          Day = Day + INT( Day_Offset )
          Day_Offset = 0_8
       ELSE
          Day_Offset = Day_Offset - INT( Days_This_Month - Day + 1, KIND=8 )
          Day = 1
          Month = Month + 1
          IF ( Month > 12 ) THEN
             Month = 1
             Year  = Year + 1
          ENDIF
       ENDIF
    ENDDO

    Hour   = INT( Seconds_Of_Day / 3600_8 )
    Minute = INT( MOD( Seconds_Of_Day, 3600_8 ) / 60_8 )
    Second = INT( MOD( Seconds_Of_Day, 60_8 ) )
    This_Date = Year * 10000 + Month * 100 + Day
    This_Time = Hour * 10000 + Minute * 100 + Second

  END SUBROUTINE Canonical_Diagnostic_Time

  PURE INTEGER FUNCTION Days_In_Month( Year, Month ) RESULT( Days )

    INTEGER, INTENT(IN) :: Year
    INTEGER, INTENT(IN) :: Month
    INTEGER, PARAMETER :: Month_Days(12) = &
         (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    LOGICAL :: Leap_Year

    Days = Month_Days(Month)
    Leap_Year = MOD( Year, 4 ) == 0 .and. &
                ( MOD( Year, 100 ) /= 0 .or. MOD( Year, 400 ) == 0 )
    IF ( Month == 2 .and. Leap_Year ) Days = 29

  END FUNCTION Days_In_Month

  SUBROUTINE Qck_Interior_Correct( Above, Target, Below )

    REAL(fp), INTENT(INOUT) :: Above, Target, Below

    REAL(fp) :: Deficit, Withdrawn

    Deficit   = -Target
    Withdrawn = MIN( Deficit, Above )
    Above     = Above - Withdrawn
    Target    = Withdrawn - Deficit
    Below     = Below + Target
    Target    = 0.0_fp

  END SUBROUTINE Qck_Interior_Correct

  SUBROUTINE Qck_Bottom_Conservative( Column, Relative_Tolerance, &
                                       Absolute_Tolerance, Status, Deficit, &
                                       Available, Withdrawn, Closure, &
                                       Tolerance, Donor_Count, &
                                       Microclosure_Tolerance )

    REAL(fp), INTENT(INOUT) :: Column(:)
    REAL(fp), INTENT(IN)    :: Relative_Tolerance, Absolute_Tolerance
    INTEGER,  INTENT(OUT)   :: Status, Donor_Count
    REAL(fp), INTENT(OUT)   :: Deficit, Available, Withdrawn, Closure, Tolerance
    REAL(fp), INTENT(IN), OPTIONAL :: Microclosure_Tolerance

    INTEGER  :: K, N_Levels
    REAL(fp) :: Remaining, Withdrawal, Microclosure_Limit

    Status      = QCK_BOTTOM_EXACT
    Donor_Count = 0
    Deficit     = 0.0_fp
    Available   = 0.0_fp
    Withdrawn   = 0.0_fp
    Closure     = 0.0_fp
    Tolerance   = 0.0_fp
    N_Levels    = SIZE( Column )
    Microclosure_Limit = 0.0_fp
    IF ( PRESENT( Microclosure_Tolerance ) ) THEN
       Microclosure_Limit = Microclosure_Tolerance
    ENDIF

    IF ( N_Levels < 2 .or. Relative_Tolerance < 0.0_fp .or. &
         Absolute_Tolerance < 0.0_fp .or. &
         Microclosure_Limit < 0.0_fp ) THEN
       Status = QCK_BOTTOM_INVALID_DONOR
       RETURN
    ENDIF
    IF ( Column(N_Levels) >= 0.0_fp ) RETURN

    Deficit   = -Column(N_Levels)
    ! Compute the positive donor inventory before mutating anything. The final
    ! tolerance is deliberately scaled by the physical preflight operands.
    DO K = N_Levels - 1, 1, -1
       IF ( Column(K) > 0.0_fp ) Available = Available + Column(K)
    ENDDO
    Tolerance = MAX( Absolute_Tolerance, &
                     Relative_Tolerance * MAX( Deficit, Available ) )

    ! A donor that is negative beyond the final tolerance is a violated
    ! precondition; a tiny negative donor is not eligible to supply mass.
    DO K = N_Levels - 1, 1, -1
       IF ( Column(K) < -Tolerance ) THEN
          Status = QCK_BOTTOM_INVALID_DONOR
          RETURN
       ENDIF
    ENDDO

    ! Simulate the deterministic donor traversal before changing the column.
    ! This protects an unfillable column from partial mutation.
    Remaining = Deficit
    DO K = N_Levels - 1, 1, -1
       IF ( Column(K) <= 0.0_fp ) CYCLE
       IF ( Column(K) >= Remaining ) THEN
          Remaining = 0.0_fp
          EXIT
       ENDIF
       Remaining = Remaining - Column(K)
    ENDDO
    IF ( Remaining > Tolerance .and. &
         ( Microclosure_Limit <= 0.0_fp .or. &
           Remaining > Microclosure_Limit ) ) THEN
       Status = QCK_BOTTOM_UNFILLABLE
       RETURN
    ENDIF

    IF ( Remaining > 0.0_fp ) THEN
       IF ( Remaining <= Tolerance ) THEN
          Status = QCK_BOTTOM_ROUNDOFF
       ELSE
          Status = QCK_BOTTOM_MICROCLOSURE
       ENDIF
       Closure = Remaining
       DO K = N_Levels - 1, 1, -1
          IF ( Column(K) <= 0.0_fp ) CYCLE
          Withdrawal = Column(K)
          Column(K) = 0.0_fp
          Withdrawn = Withdrawn + Withdrawal
          Donor_Count = Donor_Count + 1
       ENDDO
       Column(N_Levels) = 0.0_fp
       RETURN
    ENDIF

    Remaining = Deficit
    DO K = N_Levels - 1, 1, -1
       IF ( Column(K) <= 0.0_fp ) CYCLE
       Withdrawal = MIN( Remaining, Column(K) )
       Column(K) = Column(K) - Withdrawal
       Withdrawn = Withdrawn + Withdrawal
       Donor_Count = Donor_Count + 1
       Remaining = Remaining - Withdrawal
       IF ( Remaining <= 0.0_fp ) EXIT
    ENDDO
    Column(N_Levels) = 0.0_fp

  END SUBROUTINE Qck_Bottom_Conservative

END MODULE Qck_Positivity_Mod
