PROGRAM wrapper_tests
 USE Restart_Validation_Mod
 USE netCDF
 USE,INTRINSIC :: ISO_FORTRAN_ENV, ONLY:REAL32,REAL64
 USE,INTRINSIC :: IEEE_ARITHMETIC
 IMPLICIT NONE
 INTEGER :: k,f,n,rc
 LOGICAL :: valid,expect
 REAL(REAL64) :: lon(2)=[-180D0,0D0],lat(2)=[-90D0,90D0],lev(2)=[0.9D0,0.1D0]
 REAL(REAL32),POINTER :: p2(:,:),p3(:,:,:)
 REAL(REAL32),TARGET :: a2(2,2),a3(2,2,2),bad(1,2)
 CHARACTER(512) :: msg,arg
 n=0
 DO k=0,22
   CALL make_file(k)
   CALL ck(NF90_OPEN('wrapper_fixture.nc',NF90_NOWRITE,f))
   CALL Restart_Validate_File(f,lon,lat,lev,20190601,30000,valid,msg)
   expect=(k==0.OR.k==16.OR.k==17.OR.k==18.OR.k==19)
   CALL must(valid.EQV.expect,'NetCDF wrapper '//TRIM(num(k))//' '//TRIM(msg))
   CALL ck(NF90_CLOSE(f))
 ENDDO
 p2=>NULL();p3=>NULL()
 CALL must(.NOT.Restart_HCO_Valid(p2,.TRUE.,0,0,[2,2]),'null found 2D')
 CALL must(.NOT.Restart_HCO_Valid(p3,.TRUE.,0,0,[2,2,2]),'null found 3D')
 CALL must(Restart_HCO_Valid(p2,.FALSE.,0,0,[2,2]),'legacy absent 2D')
 CALL must(Restart_HCO_Valid(p3,.FALSE.,0,0,[2,2,2]),'legacy absent 3D')
 p2=>a2;p3=>a3;a2=7;a3=8
 CALL must(.NOT.Restart_HCO_Valid(p2,.TRUE.,1,0,[2,2]),'error status 2D')
 CALL must(.NOT.Restart_HCO_Valid(p3,.TRUE.,1,0,[2,2,2]),'error status 3D')
 CALL must(Restart_HCO_Valid(p2,.TRUE.,0,0,[2,2]),'valid 2D')
 CALL must(Restart_HCO_Valid(p3,.TRUE.,0,0,[2,2,2]),'valid 3D')
 p2=>bad
 CALL must(.NOT.Restart_HCO_Valid(p2,.TRUE.,0,0,[2,2]),'wrong shape 2D')
 CALL must(.NOT.Restart_HCO_Valid(p3,.TRUE.,0,0,[2,2,3]),'wrong shape 3D')
 CALL must(ALL(a2==7).AND.ALL(a3==8),'pointer validation does not mutate state')
 CALL GET_COMMAND_ARGUMENT(1,arg)
 IF(LEN_TRIM(arg)>0)THEN
   CALL ck(NF90_OPEN(TRIM(arg),NF90_NOWRITE,f))
   ! Real archived grid is read only to verify the actual metadata contract.
   BLOCK
     REAL(REAL64) :: x(144),y(91),z(47)
     INTEGER :: v
     CALL ck(NF90_INQ_VARID(f,'lon',v));CALL ck(NF90_GET_VAR(f,v,x))
     CALL ck(NF90_INQ_VARID(f,'lat',v));CALL ck(NF90_GET_VAR(f,v,y))
     CALL ck(NF90_INQ_VARID(f,'lev',v));CALL ck(NF90_GET_VAR(f,v,z))
     CALL Restart_Validate_File(f,x,y,z,20190601,30000,valid,msg)
     CALL must(valid,'actual checkpoint metadata '//TRIM(msg))
   END BLOCK
   CALL ck(NF90_CLOSE(f))
 ENDIF
 WRITE(*,*)'WRAPPER_TESTS_PASS',n
CONTAINS
 SUBROUTINE make_file(k)
 INTEGER,INTENT(IN)::k
 INTEGER :: f,ids(4),vids(4),v,i,ty,dl(4),nv
 REAL(REAL64)::x(3),y(2),z(2),t(2),field(3,2,2,2)
 CHARACTER(4)::names(4)=[CHARACTER(4)::'lon','lat','lev','time']
 CHARACTER(80)::units
 x=[lon,180D0];y=lat;z=lev;t=0;field=0
 IF(k==1)x(1)=x(1)+1
 IF(k==2)y=y(2:1:-1)
 IF(k==3)z=z(2:1:-1)
 IF(k==4)x(1)=IEEE_VALUE(1D0,IEEE_QUIET_NAN)
 IF(k==5)t(1)=60
 IF(k==19)t(1)=180
 CALL ck(NF90_CREATE('wrapper_fixture.nc',IOR(NF90_CLOBBER,NF90_NETCDF4),f))
 dl=[2,2,2,1];IF(k==6)dl(1)=3;IF(k==7)dl(4)=2
 DO i=1,4
   CALL ck(NF90_DEF_DIM(f,TRIM(names(i)),dl(i),ids(i)))
   ty=NF90_DOUBLE;IF(k==16)ty=NF90_FLOAT
   CALL ck(NF90_DEF_VAR(f,TRIM(names(i)),ty,[ids(i)],vids(i)))
   SELECT CASE(i)
   CASE(1);units='degrees_east'
   CASE(2);units='degrees_north'
   CASE(3);units='level'
   CASE(4);units='minutes since 2019-06-01 03:00:00'
   END SELECT
   IF(k==8.AND.i==1)units='radians'
   IF(k==9.AND.i==4)units='minutes since 2019-02-30 03:00:00'
   IF(k==19.AND.i==4)units='minutes since 2019-06-01 00:00:00'
   CALL ck(NF90_PUT_ATT(f,vids(i),'units',TRIM(units)))
   IF(i==4)THEN
      units='gregorian';IF(k==10)units='360_day'
      CALL ck(NF90_PUT_ATT(f,vids(i),'calendar',TRIM(units)))
   ENDIF
 ENDDO
 IF(k/=17)THEN
   CALL ck(NF90_DEF_VAR(f,'Met_DELPDRY',NF90_DOUBLE,ids,v))
   units='hPa';IF(k==11)units='Pa'
   CALL ck(NF90_PUT_ATT(f,v,'units',TRIM(units)))
 ENDIF
 IF(k==12)THEN
   CALL ck(NF90_DEF_VAR(f,'SpeciesRst_CO',NF90_DOUBLE,[ids(2),ids(1),ids(3),ids(4)],v))
 ELSEIF(k==13)THEN
   CALL ck(NF90_DEF_VAR(f,'SpeciesRst_CO',NF90_DOUBLE,[ids(1),ids(2),ids(4)],v))
 ELSE
   CALL ck(NF90_DEF_VAR(f,'SpeciesRst_CO',NF90_DOUBLE,ids,v))
 ENDIF
 units='mol mol-1 dry';IF(k==20)units='kg kg-1'
 CALL ck(NF90_PUT_ATT(f,v,'units',TRIM(units)))
 IF(k==21)THEN
   CALL ck(NF90_DEF_VAR(f,'Chem_KPPHvalue',NF90_DOUBLE,[ids(1),ids(2),ids(4)],v))
   CALL ck(NF90_PUT_ATT(f,v,'units','1'))
 ENDIF
 IF(k==22)THEN
   CALL ck(NF90_DEF_VAR(f,'Chem_WetDepNitrogen',NF90_DOUBLE,ids,v))
   CALL ck(NF90_PUT_ATT(f,v,'units','molec cm-2 s-1'))
 ENDIF
 IF(k/=18)THEN
   IF(k==14)THEN
      CALL ck(NF90_DEF_VAR(f,'Chem_JOH',NF90_DOUBLE,ids,v))
   ELSE
      CALL ck(NF90_DEF_VAR(f,'Chem_JOH',NF90_DOUBLE,[ids(1),ids(2),ids(4)],v))
   ENDIF
   units='1';IF(k==15)units='invalid'
   CALL ck(NF90_PUT_ATT(f,v,'units',TRIM(units)))
 ENDIF
 CALL ck(NF90_ENDDEF(f))
 CALL ck(NF90_PUT_VAR(f,vids(1),x(:dl(1))))
 CALL ck(NF90_PUT_VAR(f,vids(2),y))
 CALL ck(NF90_PUT_VAR(f,vids(3),z))
 CALL ck(NF90_PUT_VAR(f,vids(4),t(:dl(4))))
 CALL ck(NF90_CLOSE(f))
 END SUBROUTINE
 SUBROUTINE ck(code)
 INTEGER,INTENT(IN)::code
 IF(code/=NF90_NOERR)THEN
   PRINT *,NF90_STRERROR(code);ERROR STOP 1
 ENDIF
 END SUBROUTINE
 SUBROUTINE must(ok,label)
 LOGICAL,INTENT(IN)::ok
 CHARACTER(*),INTENT(IN)::label
 IF(.NOT.ok)THEN;PRINT *,'FAIL ',label;ERROR STOP 2;ENDIF
 n=n+1;PRINT *,'PASS ',label
 END SUBROUTINE
 FUNCTION num(i)RESULT(s)
 INTEGER,INTENT(IN)::i
 CHARACTER(10)::s
 WRITE(s,'(I0)')i
 END FUNCTION
END PROGRAM
