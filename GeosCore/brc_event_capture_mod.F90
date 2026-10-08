! Default-off, serial read-only process boundaries for native origin replay.
! These snapshots diagnose gaps; they do not define transfer operators.
MODULE BRC_EVENT_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE Input_Opt_Mod, ONLY: OptInput
  USE PhysConstants, ONLY: AIRMW
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS, GET_TS_DYN, GET_TS_CONV, GET_TS_CHEM, GET_TS_EMIS
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_EVENT_CAPTURE, BRC_EVENT_CAPTURE_STEP
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  LOGICAL, SAVE :: Initialized=.FALSE.
  INTEGER, SAVE :: Ids(7)=0, Step=0, Sequence=0, Limit=2
CONTAINS
  INTEGER FUNCTION BRC_EVENT_CAPTURE_STEP() RESULT(CaptureStep)
    CaptureStep=0
    IF (.NOT. Initialized) RETURN
    IF (LEN_TRIM(Directory)==0 .OR. Step<1 .OR. Step>Limit) RETURN
    CaptureStep=Step
  END FUNCTION BRC_EVENT_CAPTURE_STEP

  SUBROUTINE BRC_EVENT_CAPTURE(Label,Active,Input_Opt,State_Chm,State_Grid,State_Met)
    CHARACTER(LEN=*), INTENT(IN) :: Label
    LOGICAL, INTENT(IN) :: Active
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    CHARACTER(LEN=32) :: Value,StoredLabel
    CHARACTER(LEN=1200) :: Path
    INTEGER :: Status,Length,S,U
    INTEGER(int32) :: Header(29)
    REAL(fp) :: MW(7)
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_PHYSICAL_EVENT_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==-1) ERROR STOP 'BRC event capture directory truncated'
      IF (Status==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('BRC_PHYSICAL_EVENT_CAPTURE_STEPS',Value,Length,Status)
        IF (Status==-1) ERROR STOP 'BRC event capture step limit truncated'
        IF (Status==0 .AND. Length>0) THEN
          IF (VERIFY(TRIM(ADJUSTL(Value)),'0123456789')/=0) ERROR STOP 'Malformed BRC event capture limit'
          READ(Value,*,IOSTAT=Status) Limit
          IF (Status/=0) ERROR STOP 'Invalid BRC event capture limit'
        ENDIF
        IF (Limit<1 .OR. Limit>6) ERROR STOP 'BRC event capture limit must be1..6'
        IF (STORAGE_SIZE(1.0_fp)/=64) ERROR STOP 'BRC event capture currently requires REAL8'
        IF (State_Grid%NestedGrid) ERROR STOP 'BRC event capture requires global domain'
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
        ENDDO
        IF (ANY(Ids<1)) ERROR STOP 'BRC event capture requires seven parents'
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (.NOT. Input_Opt%amIRoot) RETURN
    IF (Label=='step_begin') Step=Step+1
    IF (Step>Limit) RETURN
    IF (Step<1 .OR. LEN_TRIM(Label)>32) ERROR STOP 'Invalid BRC event capture label/order'
    Sequence=Sequence+1
    Header(1:11)=[1_int32,INT(STORAGE_SIZE(1.0_fp)/8,int32),INT(Z'01020304',int32), &
         INT(State_Grid%NX,int32),INT(State_Grid%NY,int32),INT(State_Grid%NZ,int32), &
         INT(Sequence,int32),INT(Step,int32),INT(GET_NYMD(),int32),INT(GET_NHMS(),int32), &
         INT(MERGE(1,0,Active),int32)]
    Header(12:18)=INT(Ids,int32)
    DO S=1,7
      Header(18+S)=INT(State_Chm%Species(Ids(S))%Units,int32)
      MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
      IF (.NOT. ALL(IEEE_IS_FINITE(State_Chm%Species(Ids(S))%Conc))) &
           ERROR STOP 'Nonfinite BRC event parent'
    ENDDO
    Header(26:29)=INT([GET_TS_DYN(),GET_TS_CONV(),GET_TS_CHEM(),GET_TS_EMIS()],int32)
    IF (.NOT. ALL(IEEE_IS_FINITE(State_Met%AD)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%AREA_M2))) ERROR STOP 'Nonfinite BRC event geometry'
    StoredLabel=Label
    WRITE(Path,'(a,"/event_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),ACCESS='stream',FORM='unformatted',STATUS='new', &
         ACTION='write',IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot create immutable BRC event capture'
    WRITE(U,IOSTAT=Status) 'BRCPE001',Header,StoredLabel,Names,MW,REAL(AIRMW,fp), &
         State_Met%AD,State_Met%AREA_M2
    IF (Status/=0) ERROR STOP 'Cannot write BRC event header/geometry'
    DO S=1,7
      WRITE(U,IOSTAT=Status) State_Chm%Species(Ids(S))%Conc
      IF (Status/=0) ERROR STOP 'Cannot write BRC event parent'
    ENDDO
    CLOSE(U,IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot close BRC event capture'
  END SUBROUTINE BRC_EVENT_CAPTURE
END MODULE BRC_EVENT_CAPTURE_MOD
