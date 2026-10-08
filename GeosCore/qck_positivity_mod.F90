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

  REAL(fp), PARAMETER, PUBLIC :: QCK_BOTTOM_RELATIVE_TOLERANCE = 1.0e-9_fp
  REAL(fp), PARAMETER, PUBLIC :: QCK_BOTTOM_ABSOLUTE_TOLERANCE = 1.0e-30_fp

  PUBLIC :: Qck_Interior_Correct
  PUBLIC :: Qck_Bottom_Conservative

CONTAINS

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
                                       Tolerance, Donor_Count )

    REAL(fp), INTENT(INOUT) :: Column(:)
    REAL(fp), INTENT(IN)    :: Relative_Tolerance, Absolute_Tolerance
    INTEGER,  INTENT(OUT)   :: Status, Donor_Count
    REAL(fp), INTENT(OUT)   :: Deficit, Available, Withdrawn, Closure, Tolerance

    INTEGER  :: K, N_Levels
    REAL(fp) :: Remaining, Withdrawal

    Status      = QCK_BOTTOM_EXACT
    Donor_Count = 0
    Deficit     = 0.0_fp
    Available   = 0.0_fp
    Withdrawn   = 0.0_fp
    Closure     = 0.0_fp
    Tolerance   = 0.0_fp
    N_Levels    = SIZE( Column )

    IF ( N_Levels < 2 .or. Relative_Tolerance < 0.0_fp .or. &
         Absolute_Tolerance < 0.0_fp ) THEN
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
    IF ( Remaining > Tolerance ) THEN
       Status = QCK_BOTTOM_UNFILLABLE
       RETURN
    ENDIF

    IF ( Remaining > 0.0_fp ) THEN
       Status  = QCK_BOTTOM_ROUNDOFF
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

