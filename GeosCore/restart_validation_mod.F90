! Validation shared by Classic restart wrappers and standalone malformed-file tests.
MODULE Restart_Validation_Mod
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: REAL32, REAL64, INT64
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE netCDF
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: fp=REAL64
  PUBLIC :: Restart_Validate_File, Restart_HCO_Valid
  PUBLIC :: Restart_Pressure_Factors, Restart_Apply_Pressure, Restart_Restore_Surface
  INTERFACE Restart_HCO_Valid
     MODULE PROCEDURE HCO_Valid_2D, HCO_Valid_3D
  END INTERFACE
CONTAINS

  LOGICAL FUNCTION HCO_Valid_2D(ptr,found,status,success,expected) RESULT(valid)
    REAL(REAL32), POINTER, INTENT(IN) :: ptr(:,:)
    LOGICAL, INTENT(IN) :: found
    INTEGER, INTENT(IN) :: status,success,expected(2)
    valid=.FALSE.
    IF(status/=success) RETURN
    IF(found) THEN
       IF(.NOT.ASSOCIATED(ptr)) RETURN
       IF(ANY(SHAPE(ptr)/=expected)) RETURN
    ENDIF
    valid=.TRUE.
  END FUNCTION

  LOGICAL FUNCTION HCO_Valid_3D(ptr,found,status,success,expected) RESULT(valid)
    REAL(REAL32), POINTER, INTENT(IN) :: ptr(:,:,:)
    LOGICAL, INTENT(IN) :: found
    INTEGER, INTENT(IN) :: status,success,expected(3)
    valid=.FALSE.
    IF(status/=success) RETURN
    IF(found) THEN
       IF(.NOT.ASSOCIATED(ptr)) RETURN
       IF(ANY(SHAPE(ptr)/=expected)) RETURN
    ENDIF
    valid=.TRUE.
  END FUNCTION

  SUBROUTINE Restart_Validate_File(ncid,lon,lat,lev,date,clock,valid,message)
    INTEGER, INTENT(IN) :: ncid,date,clock
    REAL(fp), INTENT(IN) :: lon(:),lat(:),lev(:)
    LOGICAL, INTENT(OUT) :: valid
    CHARACTER(*), INTENT(OUT) :: message
    INTEGER :: ids(4),lens(4),i,v,nvars,rank,dims(NF90_MAX_VAR_DIMS),kind,stat
    REAL(fp) :: times(1),scale,origin_seconds,expected_seconds
    CHARACTER(4), PARAMETER :: names(4)=[CHARACTER(4)::'lon','lat','lev','time']
    CHARACTER(NF90_MAX_NAME) :: name
    CHARACTER(256) :: units,calendar
    valid=.FALSE.;message='Invalid local REAL8 restart dimensions'
    lens=[SIZE(lon),SIZE(lat),SIZE(lev),1]
    DO i=1,4
       stat=NF90_INQ_DIMID(ncid,TRIM(names(i)),ids(i));IF(stat/=NF90_NOERR) RETURN
       stat=NF90_INQUIRE_DIMENSION(ncid,ids(i),len=rank);IF(stat/=NF90_NOERR) RETURN
       IF(rank/=lens(i)) RETURN
    ENDDO
    CALL Check_Coordinate(ncid,'lon',ids(1),lon,'degrees_east',valid,message)
    IF(.NOT.valid) RETURN
    CALL Check_Coordinate(ncid,'lat',ids(2),lat,'degrees_north',valid,message)
    IF(.NOT.valid) RETURN
    CALL Check_Coordinate(ncid,'lev',ids(3),lev,'level',valid,message)
    IF(.NOT.valid) RETURN
    valid=.FALSE.;message='Invalid local REAL8 restart time coordinate'
    stat=NF90_INQ_VARID(ncid,'time',v);IF(stat/=NF90_NOERR) RETURN
    stat=NF90_INQUIRE_VARIABLE(ncid,v,xtype=kind,ndims=rank,dimids=dims)
    IF(stat/=NF90_NOERR) RETURN
    IF(rank/=1) RETURN
    IF(dims(1)/=ids(4)) RETURN
    IF(kind/=NF90_DOUBLE.AND.kind/=NF90_FLOAT) RETURN
    stat=NF90_GET_VAR(ncid,v,times);IF(stat/=NF90_NOERR) RETURN
    IF(.NOT.IEEE_IS_FINITE(times(1))) RETURN
    calendar='standard'
    stat=NF90_GET_ATT(ncid,v,'calendar',calendar)
    IF(stat/=NF90_NOERR.AND.stat/=NF90_ENOTATT) RETURN
    IF(TRIM(calendar)/='standard'.AND.TRIM(calendar)/='gregorian'.AND. &
       TRIM(calendar)/='proleptic_gregorian') RETURN
    stat=NF90_GET_ATT(ncid,v,'units',units);IF(stat/=NF90_NOERR) RETURN
    CALL Time_Origin(TRIM(units),scale,origin_seconds,valid)
    IF(.NOT.valid) RETURN
    CALL Date_Seconds(date,clock,expected_seconds,valid)
    IF(.NOT.valid) RETURN
    valid=.FALSE.
    IF(.NOT.IEEE_IS_FINITE(times(1)*scale)) RETURN
    IF(ABS(origin_seconds+times(1)*scale-expected_seconds)>1.0e-6_fp) RETURN
    ! Validate the layout of every field the local reader can load, before mutation.
    message='Invalid local REAL8 restart field layout'
    stat=NF90_INQUIRE(ncid,nVariables=nvars);IF(stat/=NF90_NOERR) RETURN
    DO v=1,nvars
       stat=NF90_INQUIRE_VARIABLE(ncid,v,name=name,xtype=kind,ndims=rank,dimids=dims)
       IF(stat/=NF90_NOERR) RETURN
       IF(INDEX(name,'SpeciesRst_')/=1.AND.INDEX(name,'Chem_')/=1.AND.name/='Met_DELPDRY') CYCLE
       message='Invalid local REAL8 restart field layout: '//TRIM(name)
       IF(rank/=3.AND.rank/=4) RETURN
       IF(rank==3) THEN
          IF(ANY(dims(1:3)/=[ids(1),ids(2),ids(4)])) RETURN
       ELSE
          IF(ANY(dims(1:4)/=ids)) RETURN
       ENDIF
       IF(kind/=NF90_DOUBLE.AND.kind/=NF90_FLOAT.AND.kind/=NF90_INT) RETURN
       IF(INDEX(name,'SpeciesRst_')==1.OR.name=='Met_DELPDRY') THEN
          IF(rank/=4) RETURN
       ENDIF
       SELECT CASE(TRIM(name))
       CASE('Chem_KPPHvalue','Chem_H2O2AfterChem','Chem_SO2AfterChem', &
            'Chem_AeroH2OSNA','Chem_ORVCsesq','Chem_StatePSC')
          IF(rank/=4) RETURN
       CASE('Chem_WetDepNitrogen','Chem_DryDepNitrogen','Chem_JOH','Chem_JNO2')
          IF(rank/=3) RETURN
       END SELECT
       IF(INDEX(name,'SpeciesRst_')==1) THEN
          IF(kind/=NF90_DOUBLE.AND.kind/=NF90_FLOAT) RETURN
          units='';stat=NF90_GET_ATT(ncid,v,'units',units);IF(stat/=NF90_NOERR) RETURN
          IF(TRIM(units)/='mol mol-1 dry') RETURN
       ENDIF
       IF(name=='Chem_JOH'.OR.name=='Chem_JNO2') THEN
          IF(rank/=3) RETURN
          units='';stat=NF90_GET_ATT(ncid,v,'units',units);IF(stat/=NF90_NOERR) RETURN
          ! Existing Classic History writes these cache fields with units '1'.
          IF(TRIM(units)/='1'.AND.TRIM(units)/='s-1') RETURN
       ENDIF
       IF(name=='Met_DELPDRY') THEN
          units='';stat=NF90_GET_ATT(ncid,v,'units',units);IF(stat/=NF90_NOERR) RETURN
          IF(TRIM(units)/='hPa') RETURN
       ENDIF
    ENDDO
    valid=.TRUE.;message=''
  END SUBROUTINE

  SUBROUTINE Check_Coordinate(ncid,name,dimid,expected,expected_units,valid,message)
    INTEGER, INTENT(IN) :: ncid,dimid
    CHARACTER(*), INTENT(IN) :: name,expected_units
    REAL(fp), INTENT(IN) :: expected(:)
    LOGICAL, INTENT(OUT) :: valid
    CHARACTER(*), INTENT(OUT) :: message
    INTEGER :: v,stat,rank,dims(NF90_MAX_VAR_DIMS),kind
    REAL(fp) :: values(SIZE(expected)),tol(SIZE(expected))
    CHARACTER(256) :: units
    valid=.FALSE.;message='Invalid local REAL8 restart coordinate: '//name
    stat=NF90_INQ_VARID(ncid,name,v);IF(stat/=NF90_NOERR) RETURN
    stat=NF90_INQUIRE_VARIABLE(ncid,v,xtype=kind,ndims=rank,dimids=dims)
    IF(stat/=NF90_NOERR) RETURN
    IF(rank/=1) RETURN
    IF(dims(1)/=dimid) RETURN
    IF(kind/=NF90_FLOAT.AND.kind/=NF90_DOUBLE) RETURN
    stat=NF90_GET_ATT(ncid,v,'units',units);IF(stat/=NF90_NOERR) RETURN
    IF(TRIM(units)/=expected_units) RETURN
    stat=NF90_GET_VAR(ncid,v,values);IF(stat/=NF90_NOERR) RETURN
    IF(ANY(.NOT.IEEE_IS_FINITE(values)).OR.ANY(.NOT.IEEE_IS_FINITE(expected))) RETURN
    tol=32*EPSILON(1.0_fp)*MAX(1.0_fp,ABS(expected))
    IF(kind==NF90_FLOAT) tol=MAX(1.0e-12_fp,2*REAL(EPSILON(1.0_REAL32),fp)*ABS(expected))
    IF(ANY(ABS(values-expected)>tol)) RETURN
    valid=.TRUE.;message=''
  END SUBROUTINE

  SUBROUTINE Time_Origin(units,scale,seconds,valid)
    CHARACTER(*),INTENT(IN) :: units
    REAL(fp),INTENT(OUT) :: scale,seconds
    LOGICAL,INTENT(OUT) :: valid
    INTEGER :: sep,y,m,d,hh,mm,ss,stat
    CHARACTER(64) :: origin
    valid=.FALSE.;scale=0;seconds=0
    sep=INDEX(units,' since ');IF(sep<2) RETURN
    SELECT CASE(units(:sep-1))
    CASE('seconds');scale=1
    CASE('minutes');scale=60
    CASE('hours');scale=3600
    CASE('days');scale=86400
    CASE DEFAULT;RETURN
    END SELECT
    origin=TRIM(units(sep+7:))
    IF(LEN_TRIM(origin)/=19) RETURN
    IF(origin(5:5)/='-'.OR.origin(8:8)/='-'.OR.origin(11:11)/=' '.OR. &
       origin(14:14)/=':'.OR.origin(17:17)/=':') RETURN
    READ(origin,'(I4,1X,I2,1X,I2,1X,I2,1X,I2,1X,I2)',IOSTAT=stat)y,m,d,hh,mm,ss
    IF(stat/=0) RETURN
    CALL Date_Seconds(y*10000+m*100+d,hh*10000+mm*100+ss,seconds,valid)
  END SUBROUTINE

  SUBROUTINE Date_Seconds(date,clock,seconds,valid)
    INTEGER,INTENT(IN) :: date,clock
    REAL(fp),INTENT(OUT) :: seconds
    LOGICAL,INTENT(OUT) :: valid
    INTEGER :: y,m,d,hh,mm,ss,mdays(12),yp
    INTEGER(INT64) :: days
    valid=.FALSE.;seconds=0
    y=date/10000;m=MOD(date/100,100);d=MOD(date,100)
    hh=clock/10000;mm=MOD(clock/100,100);ss=MOD(clock,100)
    IF(y<1900.OR.y>9999.OR.m<1.OR.m>12) RETURN
    IF(clock<0.OR.hh>23.OR.mm>59.OR.ss>59) RETURN
    mdays=[31,28,31,30,31,30,31,31,30,31,30,31]
    IF(MOD(y,4)==0.AND.(MOD(y,100)/=0.OR.MOD(y,400)==0))mdays(2)=29
    IF(d<1.OR.d>mdays(m)) RETURN
    yp=y-1
    days=365_INT64*yp+yp/4-yp/100+yp/400+SUM(mdays(:m-1))+d-1
    seconds=REAL(days,fp)*86400+hh*3600+mm*60+ss
    valid=.TRUE.
  END SUBROUTINE

#include "restart_pressure_kernel.inc"
#include "surface_cache_kernel.inc"
END MODULE Restart_Validation_Mod
