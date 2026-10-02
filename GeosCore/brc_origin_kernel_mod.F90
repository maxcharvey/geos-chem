! Bookkeeping only: apply the actual parent remaining fraction, including its
! cutoff, to diagnostic origin states. Never normalize independent partitions.
MODULE BRC_ORIGIN_KERNEL_MOD
  USE Precision_Mod, ONLY: fp
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_ORIGIN_SPLIT_LOSS
CONTAINS
  PURE SUBROUTINE BRC_ORIGIN_SPLIT_LOSS(Before, After, Parts, Lost, Status)
    REAL(fp), INTENT(IN) :: Before, After
    REAL(fp), INTENT(INOUT) :: Parts(:)
    REAL(fp), INTENT(OUT) :: Lost(SIZE(Parts))
    INTEGER, INTENT(OUT) :: Status
    REAL(fp) :: Remain
    Lost = 0.0_fp
    Status = 0
    IF (.NOT. IEEE_IS_FINITE(Before) .OR. .NOT. IEEE_IS_FINITE(After) .OR. &
        ANY(.NOT. IEEE_IS_FINITE(Parts))) THEN
      Status = 1
      RETURN
    ENDIF
    IF (Before < 0.0_fp .OR. After < 0.0_fp .OR. After > Before .OR. &
        ANY(Parts < 0.0_fp)) THEN
      Status = 2
      RETURN
    ENDIF
    ! A zero carrier cannot supply a conversion. Preserve any origin residual
    ! for the closure audit; deleting it would hide nonadditive transport.
    IF (Before == 0.0_fp) RETURN
    Remain = After / Before
    Lost = Parts - Parts * Remain
    Parts = Parts * Remain
  END SUBROUTINE BRC_ORIGIN_SPLIT_LOSS
END MODULE BRC_ORIGIN_KERNEL_MOD
