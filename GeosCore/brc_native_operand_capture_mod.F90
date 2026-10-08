! Serial, default-off native operands. No physical or country-state writes.
MODULE BRC_NATIVE_OPERAND_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE Input_Opt_Mod, ONLY: OptInput
  USE PhysConstants, ONLY: AIRMW
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_NATIVE_OPERAND_CAPTURE
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  LOGICAL, SAVE :: Initialized=.FALSE.
  INTEGER, SAVE :: Ids(7)=0, Sequence=0
CONTAINS
  SUBROUTINE BRC_NATIVE_OPERAND_CAPTURE(Kind,Phase,Caller,Input_Opt,State_Chm, &
       State_Grid,State_Met,Mapping,UpdateMR,Applied,PreviousUnits)
    INTEGER, INTENT(IN) :: Kind,Phase,Mapping(:)
    CHARACTER(LEN=*), INTENT(IN) :: Caller
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    LOGICAL, OPTIONAL, INTENT(IN) :: UpdateMR,Applied
    INTEGER, OPTIONAL, INTENT(IN) :: PreviousUnits
    CHARACTER(LEN=32) :: StoredCaller
    CHARACTER(LEN=1200) :: Path
    INTEGER :: Step,Status,Length,S,U
    INTEGER(int32) :: Header(37)
    REAL(fp) :: MW(7),OldSum,NewSum
    Step=BRC_EVENT_CAPTURE_STEP()
    IF (Step==0 .OR. .NOT. Input_Opt%amIRoot) RETURN
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_NATIVE_OPERAND_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==-1) ERROR STOP 'BRC native operand capture directory truncated'
      IF (Status==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
        IF (STORAGE_SIZE(1.0_fp)/=64 .OR. State_Grid%NestedGrid) &
             ERROR STOP 'BRC native operand capture requires REAL8/global'
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
        ENDDO
        IF (ANY(Ids<1)) ERROR STOP 'BRC native operand capture requires seven parents'
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (Kind<1 .OR. Kind>3 .OR. Phase<1 .OR. Phase>4 .OR. &
        (Kind/=2 .AND. Phase>2) .OR. LEN_TRIM(Caller)>32 .OR. &
        SIZE(Mapping)<1 .OR. ANY(Mapping<1) .OR. ANY(Mapping>State_Chm%nSpecies)) &
         ERROR STOP 'Invalid BRC native operand identity'
    Sequence=Sequence+1
    Header(1:16)=[1_int32,8_int32,INT(Z'01020304',int32), &
         INT(State_Grid%NX,int32),INT(State_Grid%NY,int32),INT(State_Grid%NZ,int32), &
         INT(Sequence,int32),INT(Step,int32),INT(GET_NYMD(),int32),INT(GET_NHMS(),int32), &
         INT(Kind,int32),INT(Phase,int32),-1_int32,-1_int32,-1_int32,INT(SIZE(Mapping),int32)]
    IF (PRESENT(UpdateMR)) Header(13)=INT(MERGE(1,0,UpdateMR),int32)
    IF (PRESENT(Applied)) Header(14)=INT(MERGE(1,0,Applied),int32)
    IF (PRESENT(PreviousUnits)) Header(15)=INT(PreviousUnits,int32)
    Header(17:23)=INT(Ids,int32)
    DO S=1,7
      Header(23+S)=INT(State_Chm%Species(Ids(S))%Units,int32)
      Header(30+S)=INT(MERGE(1,0,ANY(Mapping==Ids(S))),int32)
      MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
      IF (.NOT. ALL(IEEE_IS_FINITE(State_Chm%Species(Ids(S))%Conc))) &
           ERROR STOP 'Nonfinite BRC native operand parent'
    ENDDO
    IF (.NOT. ALL(IEEE_IS_FINITE(State_Met%AD)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%AREA_M2)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%DP_DRY_PREV)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%DELP_DRY))) &
         ERROR STOP 'Nonfinite BRC native operand geometry'
    OldSum=SUM(State_Met%DP_DRY_PREV)
    NewSum=SUM(State_Met%DELP_DRY)
    StoredCaller=Caller
    WRITE(Path,'(a,"/native_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),ACCESS='stream',FORM='unformatted',STATUS='new', &
         ACTION='write',IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot create immutable BRC native operand capture'
    WRITE(U,IOSTAT=Status) 'BRCNO001',Header,StoredCaller,Names,MW,REAL(AIRMW,fp), &
         INT(Mapping,int32),OldSum,NewSum,State_Met%AD,State_Met%AREA_M2, &
         State_Met%DP_DRY_PREV,State_Met%DELP_DRY
    IF (Status/=0) ERROR STOP 'Cannot write BRC native operand header/geometry'
    DO S=1,7
      WRITE(U,IOSTAT=Status) State_Chm%Species(Ids(S))%Conc
      IF (Status/=0) ERROR STOP 'Cannot write BRC native operand parent'
    ENDDO
    CLOSE(U,IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot close BRC native operand capture'
  END SUBROUTINE BRC_NATIVE_OPERAND_CAPTURE
END MODULE BRC_NATIVE_OPERAND_CAPTURE_MOD
