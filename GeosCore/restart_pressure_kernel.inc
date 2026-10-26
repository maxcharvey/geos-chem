! Shared production/fixture kernel. Pressure inputs are in the same dry-hPa units.
SUBROUTINE Restart_Pressure_Factors( old_dp, new_dp, factors, valid )
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  REAL(fp), INTENT(IN) :: old_dp(:,:,:), new_dp(:,:,:)
  REAL(fp), INTENT(OUT) :: factors(:,:,:)
  LOGICAL, INTENT(OUT) :: valid
  valid=.FALSE.
  IF ( ANY(SHAPE(old_dp)/=SHAPE(new_dp)) ) RETURN
  IF ( ANY(SHAPE(factors)/=SHAPE(old_dp)) ) RETURN
  IF ( ANY(.NOT.IEEE_IS_FINITE(old_dp)) ) RETURN
  IF ( ANY(.NOT.IEEE_IS_FINITE(new_dp)) ) RETURN
  IF ( ANY(old_dp<=0.0_fp) ) RETURN
  IF ( ANY(new_dp<=0.0_fp) ) RETURN
  factors=1.0_fp
  WHERE (old_dp/=new_dp) factors=old_dp/new_dp
  IF ( ANY(.NOT.IEEE_IS_FINITE(factors)) ) RETURN
  IF ( ANY(factors<=0.0_fp) ) RETURN
  valid=.TRUE.
END SUBROUTINE Restart_Pressure_Factors

SUBROUTINE Restart_Apply_Pressure( conc, old_dp, new_dp, factors )
  IMPLICIT NONE
  REAL(fp), INTENT(INOUT) :: conc(:,:,:)
  REAL(fp), INTENT(IN) :: old_dp(:,:,:),new_dp(:,:,:),factors(:,:,:)
  INTEGER :: i,j,l
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l) COLLAPSE(3)
  DO l=1,SIZE(conc,3)
  DO j=1,SIZE(conc,2)
  DO i=1,SIZE(conc,1)
    ! Equal-pressure cells are literal no-ops, even in a mixed grid.
    IF (old_dp(i,j,l)/=new_dp(i,j,l)) conc(i,j,l)=conc(i,j,l)*factors(i,j,l)
  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
END SUBROUTINE Restart_Apply_Pressure
