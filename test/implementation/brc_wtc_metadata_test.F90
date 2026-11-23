MODULE WTC_METADATA_FIXTURE
 IMPLICIT NONE
 TYPE OptInput
 LOGICAL :: LBRC=.true.
 END TYPE
 TYPE Species
 LOGICAL :: Is_Gas=.false.,Is_Aerosol=.false.,Is_Omitted=.false.,Is_Kpp=.false.
 LOGICAL :: Is_Advected=.false.,Is_DryDep=.false.,Is_WetDep=.false.
 LOGICAL :: Is_DryAlt=.false.,Is_HygroGrowth=.false.
 END TYPE
CONTAINS
 SUBROUTINE Validate(Input_Opt,ThisSpc,Name,RC)
 TYPE(OptInput),INTENT(IN) :: Input_Opt
 TYPE(Species),INTENT(INOUT) :: ThisSpc
 CHARACTER(*),INTENT(IN) :: Name
 INTEGER,INTENT(OUT) :: RC
 CHARACTER(255) :: errMsg,thisLoc,dbSpcName,modelSpcName
 RC=0;thisLoc='fixture';dbSpcName=Name;modelSpcName=Name
 INCLUDE 'wtc_metadata_under_test.inc'
 END SUBROUTINE
 SUBROUTINE GC_Error(Message,RC,Location)
 CHARACTER(*),INTENT(IN) :: Message,Location
 INTEGER,INTENT(OUT) :: RC
 RC=1
 END SUBROUTINE
END MODULE
PROGRAM WTC_METADATA_TEST
 USE WTC_METADATA_FIXTURE
 IMPLICIT NONE
 TYPE(OptInput) :: Opt
 TYPE(Species) :: Info
 INTEGER :: RC,Mode
 Info=Species();Info%Is_HygroGrowth=.true.
 CALL Validate(Opt,Info,'WTC',RC)
 IF(RC/=0.OR.Info%Is_Gas.OR.Info%Is_Aerosol) ERROR STOP 'Diagnostic metadata rejected or made physical'
 DO Mode=1,8
    Info=Species();Info%Is_HygroGrowth=.true.
    SELECT CASE(Mode)
    CASE(1);Info%Is_Gas=.true.
    CASE(2);Info%Is_Aerosol=.true.
    CASE(3);Info%Is_Advected=.true.
    CASE(4);Info%Is_Kpp=.true.
    CASE(5);Info%Is_DryDep=.true.
    CASE(6);Info%Is_WetDep=.true.
    CASE(7);Info%Is_DryAlt=.true.
    CASE(8);Info%Is_HygroGrowth=.false.
    END SELECT
    CALL Validate(Opt,Info,'WTC',RC)
    IF(RC==0) ERROR STOP 'Invalid WTC physical metadata accepted'
 ENDDO
 Info=Species();Info%Is_Gas=.true.
 CALL Validate(Opt,Info,'CO',RC)
 IF(RC/=0) ERROR STOP 'Ordinary gas rejected'
 Info=Species();Info%Is_Aerosol=.true.
 CALL Validate(Opt,Info,'OCPI',RC)
 IF(RC/=0) ERROR STOP 'Ordinary aerosol rejected'
 Info=Species()
 CALL Validate(Opt,Info,'UNKNOWN',RC)
 IF(RC==0) ERROR STOP 'Ordinary invalid metadata accepted'
 Info=Species();Info%Is_Kpp=.true.
 CALL Validate(Opt,Info,'KPP_FAMILY',RC)
 IF(RC/=0.OR..NOT.Info%Is_Gas) ERROR STOP 'Ordinary KPP fallback changed'
 WRITE(*,'(a)') 'PASS: actual WTC metadata validation and unchanged ordinary gas/aerosol rules'
END PROGRAM
