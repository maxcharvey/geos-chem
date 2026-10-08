! Four serial callsite pairs; observer has no physical or country writes.
MODULE BRC_CHEM_UNITS_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE Input_Opt_Mod, ONLY: OptInput
  USE PhysConstants, ONLY: AIRMW,AVO
  USE Time_Mod, ONLY: GET_NYMD,GET_NHMS
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_CHEM_UNITS_CAPTURE
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  CHARACTER(LEN=32), SAVE :: StoredCaller=''
  LOGICAL, SAVE :: Initialized=.FALSE.,PairOpen=.FALSE.
  INTEGER, SAVE :: Ids(7)=0,Sequence=0,BeforeUnits=-1,StoredTarget=-1,StoredMapKind=0
CONTAINS
  SUBROUTINE BRC_CHEM_UNITS_CAPTURE(Caller,Phase,Target,Input_Opt,State_Chm, &
       State_Grid,State_Met,PreviousUnits,Mapping)
    CHARACTER(LEN=*), INTENT(IN) :: Caller
    INTEGER, INTENT(IN) :: Phase,Target
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    INTEGER, OPTIONAL, INTENT(IN) :: PreviousUnits
    INTEGER, OPTIONAL, INTENT(IN) :: Mapping(:)
    INTEGER, ALLOCATABLE :: ObservedMapping(:)
    INTEGER :: Step,Status,Length,S,U,Current,MapKind
    INTEGER(int32) :: H(39)
    REAL(fp) :: MW(7)
    CHARACTER(LEN=32) :: Label
    CHARACTER(LEN=1200) :: Path
    Step=BRC_EVENT_CAPTURE_STEP()
    IF (Step==0 .OR. .NOT. Input_Opt%amIRoot) RETURN
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_CHEM_UNITS_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==-1) ERROR STOP 'BRC chemistry unit directory truncated'
      IF (Status==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
#if defined(TOMAS) || defined(APM) || defined(ADJOINT)
        ERROR STOP 'Unsupported BRC chemistry unit observer build'
#endif
        IF (STORAGE_SIZE(1.0_fp)/=64 .OR. State_Grid%NestedGrid) &
             ERROR STOP 'BRC chemistry units require REAL8/global'
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
        ENDDO
        IF (ANY(Ids<=0)) ERROR STOP 'BRC chemistry units require seven parents'
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    MapKind=1
    IF (PRESENT(Mapping)) THEN
      MapKind=2
      IF (SIZE(Mapping)/=State_Chm%nAdvect .OR. SIZE(Mapping)/=SIZE(State_Chm%Map_Advect)) &
           ERROR STOP 'Unsupported advected chemistry mapping size'
      IF (ANY(Mapping/=State_Chm%Map_Advect)) ERROR STOP 'Unbound advected chemistry mapping'
      ObservedMapping=Mapping
    ELSE
      ObservedMapping=State_Chm%Map_All
    ENDIF
    IF (Phase<1 .OR. Phase>2 .OR. LEN_TRIM(Caller)>32 .OR. SIZE(ObservedMapping)<1 .OR. &
        ANY(ObservedMapping<1) .OR. ANY(ObservedMapping>State_Chm%nSpecies)) &
         ERROR STOP 'Unsupported BRC chemistry units identity/map'
    SELECT CASE(Caller)
    CASE('chem_to_kg','fullchem_to_mnd','fullchem_restore','chem_restore')
      IF (MapKind/=1 .OR. SIZE(ObservedMapping)/=State_Chm%nSpecies) &
           ERROR STOP 'Expected full chemistry mapping'
    CASE('sulfate_early_to_vv','sulfate_early_restore','sulfate_late_to_vv','sulfate_late_restore')
      IF (MapKind/=2) ERROR STOP 'Expected actual sulfate advected mapping'
    CASE DEFAULT
      ERROR STOP 'Unknown BRC chemistry unit caller'
    END SELECT
    DO S=1,SIZE(ObservedMapping)
      IF (COUNT(ObservedMapping==ObservedMapping(S))/=1) ERROR STOP 'Duplicate chemistry unit mapping'
    ENDDO
    IF (.NOT. ALL([(ANY(ObservedMapping==Ids(S)),S=1,7)])) &
         ERROR STOP 'Unmapped chemistry parent'
    IF (.NOT. IEEE_IS_FINITE(REAL(AIRMW,fp)) .OR. AIRMW<=0 .OR. &
        .NOT. IEEE_IS_FINITE(REAL(AVO,fp)) .OR. AVO<=0) &
         ERROR STOP 'Invalid chemistry physical constants'
    Current=State_Chm%Species(ObservedMapping(1))%Units
    IF (.NOT. ANY(Current==[1,2,5,6]) .OR. .NOT. ANY(Target==[1,2,5,6])) &
         ERROR STOP 'Unsupported BRC chemistry unit conversion'
    IF (Phase==1) THEN
      IF (PairOpen .OR. PRESENT(PreviousUnits)) ERROR STOP 'Invalid chemistry unit before capture'
      PairOpen=.TRUE.;BeforeUnits=Current;StoredTarget=Target;StoredCaller=Caller;StoredMapKind=MapKind
    ELSE
      IF (.NOT. PairOpen .OR. StoredCaller/=Caller .OR. StoredTarget/=Target .OR. StoredMapKind/=MapKind .OR. Current/=Target) &
           ERROR STOP 'Unmatched chemistry unit pair'
      PairOpen=.FALSE.
      IF (PRESENT(PreviousUnits)) THEN
        IF (PreviousUnits/=BeforeUnits) ERROR STOP 'Chemistry previous unit binding failed'
      ENDIF
    ENDIF
    Sequence=Sequence+1
    H(1:16)=[2_int32,8_int32,INT(Z'01020304',int32), &
         INT(State_Grid%NX,int32),INT(State_Grid%NY,int32),INT(State_Grid%NZ,int32), &
         INT(Sequence,int32),INT(Step,int32),INT(GET_NYMD(),int32),INT(GET_NHMS(),int32), &
         INT(Phase,int32),INT(Target,int32),INT(BeforeUnits,int32),-1_int32, &
         INT(SIZE(ObservedMapping),int32),INT(State_Chm%nSpecies,int32)]
    IF (PRESENT(PreviousUnits)) H(14)=INT(PreviousUnits,int32)
    H(17:23)=INT(Ids,int32)
    DO S=1,7
      H(23+S)=INT(State_Chm%Species(Ids(S))%Units,int32)
      H(30+S)=INT(MERGE(1,0,ANY(ObservedMapping==Ids(S))),int32)
      MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
      IF (.NOT. IEEE_IS_FINITE(MW(S)) .OR. MW(S)<=0) ERROR STOP 'Invalid chemistry molecular weight'
      IF (H(23+S)/=Current .OR. .NOT. ALL(IEEE_IS_FINITE(State_Chm%Species(Ids(S))%Conc))) &
           ERROR STOP 'Invalid BRC chemistry unit parent'
    ENDDO
    H(38)=INT(MERGE(1,0,BeforeUnits==1 .AND. Target==5),int32)
    IF (.NOT. ALL(IEEE_IS_FINITE(State_Met%AD)) .OR. ANY(State_Met%AD<=0) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%AIRVOL)) .OR. ANY(State_Met%AIRVOL<=0)) &
         ERROR STOP 'Invalid BRC chemistry unit geometry'
    H(39)=INT(MapKind,int32)
    Label=Caller
    WRITE(Path,'(a,"/chemunit_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),ACCESS='stream',FORM='unformatted',STATUS='new', &
         ACTION='write',IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot create immutable chemistry unit capture'
    WRITE(U,IOSTAT=Status) 'BRCCU002',H,Label,Names,MW,REAL(AIRMW,fp),REAL(AVO,fp), &
         INT(ObservedMapping,int32),State_Met%AD,State_Met%AIRVOL
    IF (Status/=0) ERROR STOP 'Cannot write chemistry unit geometry'
    DO S=1,7
      WRITE(U,IOSTAT=Status) State_Chm%Species(Ids(S))%Conc
      IF (Status/=0) ERROR STOP 'Cannot write chemistry unit parent'
    ENDDO
    CLOSE(U,IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot close chemistry unit capture'
  END SUBROUTINE BRC_CHEM_UNITS_CAPTURE
END MODULE BRC_CHEM_UNITS_CAPTURE_MOD
