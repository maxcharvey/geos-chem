SUBROUTINE Restart_Restore_Surface(joh,jno2,zpj,io3,ino2,valid)
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  REAL(fp),INTENT(IN)::joh(:,:),jno2(:,:)
  REAL(fp),INTENT(INOUT)::zpj(:,:,:,:)
  INTEGER,INTENT(IN)::io3,ino2
  LOGICAL,INTENT(OUT)::valid
  valid=.FALSE.
  IF (SIZE(zpj,1)<1) RETURN
  IF (io3<1.OR.ino2<1.OR.io3==ino2) RETURN
  IF (io3>SIZE(zpj,2).OR.ino2>SIZE(zpj,2)) RETURN
  IF (ANY(SHAPE(joh)/=[SIZE(zpj,3),SIZE(zpj,4)])) RETURN
  IF (ANY(SHAPE(jno2)/=SHAPE(joh))) RETURN
  IF (ANY(.NOT.IEEE_IS_FINITE(joh))) RETURN
  IF (ANY(.NOT.IEEE_IS_FINITE(jno2))) RETURN
  IF (ANY(joh<0.0_fp).OR.ANY(jno2<0.0_fp)) RETURN
  zpj(1,io3,:,:)=joh
  zpj(1,ino2,:,:)=jno2
  valid=.TRUE.
END SUBROUTINE Restart_Restore_Surface
