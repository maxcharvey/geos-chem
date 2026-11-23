MODULE kernels
 USE,INTRINSIC::ISO_FORTRAN_ENV,ONLY:INT64
 IMPLICIT NONE
 INTEGER,PARAMETER::fp=KIND(1.0D0)
CONTAINS
#include "restart_pressure_kernel.inc"
#include "surface_cache_kernel.inc"
END MODULE
PROGRAM kernel_tests
 USE kernels
 USE,INTRINSIC::IEEE_ARITHMETIC
 IMPLICIT NONE
 REAL(fp)::old(2,2,2),new(2,2,2),fac(2,2,2),c(2,2,2),orig(2,2,2),err
 REAL(fp)::small(1,1,1),fsmall(1,1,1),z(2,4,2,2),zold(2,4,2,2),joh(2,2),jno2(2,2),badj(1,1)
 LOGICAL::valid
 INTEGER::i,side,kindbad,n
 n=0
 orig=RESHAPE([1.2345678901234567D-8,1D-300,-0D0,3D-9,1D-6,1D-12,0D0,1D-7],[2,2,2])
 old=10;new=old;c=orig
 CALL Restart_Pressure_Factors(old,new,fac,valid)
 CALL must(valid,'equal pressures valid');CALL Restart_Apply_Pressure(c,old,new,fac)
 CALL must(ALL(TRANSFER(c,[0_INT64],8)==TRANSFER(orig,[0_INT64],8)),'equal exact words')
 old(1,:,:)=10;old(2,:,:)=20;new=15;c=orig
 CALL must(SUM(old)==SUM(new),'global cancellation fixture')
 CALL Restart_Pressure_Factors(old,new,fac,valid);CALL must(valid,'cancelling local differences valid')
 CALL Restart_Apply_Pressure(c,old,new,fac)
 err=MAXVAL(ABS(c*new-orig*old)/MAX(ABS(orig*old),TINY(1D0)))
 CALL must(err<=64*EPSILON(1D0),'local species mass conserved')
 new=old;new(1,1,1)=old(1,1,1)/2;c=orig
 CALL Restart_Pressure_Factors(old,new,fac,valid);CALL Restart_Apply_Pressure(c,old,new,fac)
 CALL must(c(1,1,1)==2*orig(1,1,1),'mixed changed cell')
 c(1,1,1)=orig(1,1,1);CALL must(ALL(TRANSFER(c,[0_INT64],8)==TRANSFER(orig,[0_INT64],8)),'mixed equal cells exact')
 DO side=1,2
 DO kindbad=1,4
  old=10;new=10
  SELECT CASE(kindbad)
  CASE(1);err=0D0
  CASE(2);err=-1D0
  CASE(3);err=IEEE_VALUE(1D0,IEEE_QUIET_NAN)
  CASE(4);err=IEEE_VALUE(1D0,IEEE_POSITIVE_INF)
  END SELECT
  IF(side==1)THEN;old(1,1,1)=err;ELSE;new(1,1,1)=err;ENDIF
  CALL Restart_Pressure_Factors(old,new,fac,valid);CALL must(.NOT.valid,'invalid pressure rejected')
 ENDDO
 ENDDO
 old=1D100;new=1D-300;CALL Restart_Pressure_Factors(old,new,fac,valid);CALL must(.NOT.valid,'overflow ratio rejected')
 old=1D-300;new=1D100;CALL Restart_Pressure_Factors(old,new,fac,valid);CALL must(.NOT.valid,'zero-underflow ratio rejected')
 small=10;CALL Restart_Pressure_Factors(small,new,fsmall,valid);CALL must(.NOT.valid,'shape mismatch rejected')
 CALL Restart_Pressure_Factors(small,small,fsmall,valid);CALL must(valid,'reduced grid accepted')
 z=3D-6;zold=z;joh=1D-6;jno2=2D-5
 CALL Restart_Restore_Surface(joh,jno2,z,2,4,valid);CALL must(valid,'saved daytime cache valid')
 CALL must(ALL(z(1,2,:,:)==joh).AND.ALL(z(1,4,:,:)==jno2),'surface channels restored')
 z(1,2,:,:)=zold(1,2,:,:);z(1,4,:,:)=zold(1,4,:,:)
 CALL must(ALL(TRANSFER(z,[0_INT64],32)==TRANSFER(zold,[0_INT64],32)),'other channels and levels untouched')
 z=zold;CALL Restart_Restore_Surface(joh,jno2,z,0,4,valid);CALL must(.NOT.valid.AND.ALL(z==zold),'invalid mapping no mutation')
 CALL Restart_Restore_Surface(joh,jno2,z,2,2,valid);CALL must(.NOT.valid.AND.ALL(z==zold),'duplicate mapping no mutation')
 badj=0;CALL Restart_Restore_Surface(badj,jno2,z,2,4,valid);CALL must(.NOT.valid.AND.ALL(z==zold),'cache shape no mutation')
 joh(1,1)=IEEE_VALUE(1D0,IEEE_QUIET_NAN);CALL Restart_Restore_Surface(joh,jno2,z,2,4,valid)
 CALL must(.NOT.valid.AND.ALL(z==zold),'invalid cache no mutation')
 joh=0;jno2=0;CALL Restart_Restore_Surface(joh,jno2,z,2,4,valid)
 CALL must(valid.AND.ALL(z(1,2,:,:)==0).AND.ALL(z(1,4,:,:)==0),'legacy cold zero cache explicit')
 WRITE(*,*)'KERNEL_FIXTURES_PASS',n
CONTAINS
 SUBROUTINE must(ok,label)
  LOGICAL,INTENT(IN)::ok
  CHARACTER(*),INTENT(IN)::label
  IF(.NOT.ok)THEN;WRITE(*,*)'FAIL ',label;ERROR STOP 2;ENDIF
  n=n+1;WRITE(*,*)'PASS ',label
 END SUBROUTINE
END PROGRAM
