! Default-off observer: only diagnostic scratch and immutable files are written.
MODULE BRC_CHEM_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE Input_Opt_Mod, ONLY: OptInput
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS, GET_TS_CHEM
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_CHEM_SNAPSHOT, BRC_CHEM_CELL_BEFORE, BRC_CHEM_CELL_AFTER
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  LOGICAL, SAVE :: Initialized=.FALSE., Enabled=.FALSE., CalledStage=.FALSE.
  INTEGER, SAVE :: Sequence=0, CurrentStage=0, Ids(7)=0
  REAL(fp), SAVE :: StoredSmall=0.0_fp
  ! TC0, raw CNEW, final CNEW, K, FREQ, RKT, gain, factor, incoming.
  REAL(fp), ALLOCATABLE, SAVE :: Values(:,:,:,:)
  ! before, after, rate present, RKT present, cutoff evaluated, cutoff true.
  INTEGER(int32), ALLOCATABLE, SAVE :: Masks(:,:,:,:)
CONTAINS
  SUBROUTINE BRC_CHEM_SNAPSHOT(Stage,Phase,Called,SourceId,DestId, &
       Input_Opt,State_Chm,State_Grid,Small,OMOC,C1,C2,C3,C4)
    INTEGER, INTENT(IN) :: Stage,Phase,SourceId,DestId
    LOGICAL, INTENT(IN) :: Called
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    REAL(fp), INTENT(IN) :: Small,OMOC,C1(:,:,:),C2(:,:,:),C3(:,:,:),C4(:,:,:)
    INTEGER :: Step,Status,Length,S,U,F,NX,NY,NZ
    INTEGER(int32) :: Header(37)
    REAL(fp) :: MW(7)
    REAL(fp), ALLOCATABLE :: Missing(:,:,:)
    CHARACTER(LEN=1200) :: Path
    Step=BRC_EVENT_CAPTURE_STEP()
    Enabled=.FALSE.
    IF (Step==0 .OR. .NOT. Input_Opt%amIRoot) RETURN
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_CHEM_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==-1) ERROR STOP 'BRC chemistry capture directory truncated'
      IF (Status==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
        IF (STORAGE_SIZE(1.0_fp)/=64 .OR. State_Grid%NestedGrid) &
             ERROR STOP 'BRC chemistry capture requires REAL8/global'
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
        ENDDO
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (Stage<0 .OR. Stage>5 .OR. Phase<1 .OR. Phase>2 .OR. &
        Small<=0 .OR. OMOC<=0 .OR. .NOT. IEEE_IS_FINITE(Small) .OR. &
        .NOT. IEEE_IS_FINITE(OMOC)) ERROR STOP 'Invalid BRC chemistry capture identity'
    NX=State_Grid%NX; NY=State_Grid%NY; NZ=State_Grid%NZ
    IF (.NOT. ALLOCATED(Values)) THEN
      ALLOCATE(Values(NX,NY,NZ,9),Masks(NX,NY,NZ,6),STAT=Status)
      IF (Status/=0) ERROR STOP 'BRC chemistry capture allocation failed'
    ENDIF
    IF (ANY(SHAPE(Values)/=[NX,NY,NZ,9]) .OR. &
        ANY(SHAPE(C1)/=[NX,NY,NZ]) .OR. ANY(SHAPE(C2)/=[NX,NY,NZ]) .OR. &
        ANY(SHAPE(C3)/=[NX,NY,NZ]) .OR. ANY(SHAPE(C4)/=[NX,NY,NZ])) &
         ERROR STOP 'BRC chemistry capture shape mismatch'
    IF (Phase==1 .OR. Stage==0) THEN
      Values=0.0_fp; Masks=0_int32
      CurrentStage=Stage; CalledStage=Called; StoredSmall=Small
    ELSE
      IF (Stage/=CurrentStage .OR. (Called .NEQV. CalledStage)) &
           ERROR STOP 'BRC chemistry capture unmatched stage'
      IF (Called .AND. (ANY(Masks(:,:,:,1)/=1) .OR. ANY(Masks(:,:,:,2)/=1))) &
           ERROR STOP 'BRC chemistry capture incomplete cell coverage'
    ENDIF
    Enabled=Called .AND. Stage>0 .AND. Phase==1
    Sequence=Sequence+1
    Header(1:16)=[1_int32,8_int32,INT(Z'01020304',int32), &
         INT(NX,int32),INT(NY,int32),INT(NZ,int32),INT(Sequence,int32), &
         INT(Step,int32),INT(GET_NYMD(),int32),INT(GET_NHMS(),int32), &
         INT(Stage,int32),INT(Phase,int32),INT(MERGE(1,0,Called),int32), &
         INT(Input_Opt%BrC_Bleach_Scheme,int32),INT(SourceId,int32),INT(DestId,int32)]
    Header(17:23)=INT(Ids,int32); MW=0.0_fp
    DO S=1,7
      Header(23+S)=-1_int32; Header(30+S)=INT(MERGE(1,0,Ids(S)>0),int32)
      IF (Ids(S)>0) THEN
        Header(23+S)=INT(State_Chm%Species(Ids(S))%Units,int32)
        MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
        IF (Header(23+S)/=1 .OR. .NOT. ALL(IEEE_IS_FINITE(State_Chm%Species(Ids(S))%Conc))) &
             ERROR STOP 'BRC chemistry capture requires finite kg parents'
      ENDIF
    ENDDO
    IF (.NOT. ALL(IEEE_IS_FINITE(C1)) .OR. .NOT. ALL(IEEE_IS_FINITE(C2)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(C3)) .OR. .NOT. ALL(IEEE_IS_FINITE(C4)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(Values))) ERROR STOP 'Nonfinite BRC chemistry captured operand'
    WRITE(Path,'(a,"/chem_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),ACCESS='stream',FORM='unformatted',STATUS='new', &
         ACTION='write',IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot create immutable BRC chemistry capture'
    WRITE(U,IOSTAT=Status) 'BRCCC001',Header,Names,MW,REAL(GET_TS_CHEM(),fp),Small,OMOC
    IF (Status/=0) ERROR STOP 'Cannot write BRC chemistry header'
    IF (ANY(Ids<=0)) THEN
      ALLOCATE(Missing(NX,NY,NZ)); Missing=0.0_fp
    ENDIF
    DO S=1,7
      IF (Ids(S)>0) THEN
        WRITE(U,IOSTAT=Status) State_Chm%Species(Ids(S))%Conc
      ELSE
        WRITE(U,IOSTAT=Status) Missing
      ENDIF
      IF (Status/=0) ERROR STOP 'Cannot write BRC chemistry parent'
    ENDDO
    WRITE(U,IOSTAT=Status) C1,C2,C3,C4
    IF (Status/=0) ERROR STOP 'Cannot write BRC chemistry pending arrays'
    DO F=1,9
      WRITE(U,IOSTAT=Status) Values(:,:,:,F)
      IF (Status/=0) ERROR STOP 'Cannot write BRC chemistry cell operand'
    ENDDO
    DO F=1,6
      WRITE(U,IOSTAT=Status) Masks(:,:,:,F)
      IF (Status/=0) ERROR STOP 'Cannot write BRC chemistry cell mask'
    ENDDO
    CLOSE(U,IOSTAT=Status)
    IF (Status/=0) ERROR STOP 'Cannot close BRC chemistry capture'
  END SUBROUTINE BRC_CHEM_SNAPSHOT

  SUBROUTINE BRC_CHEM_CELL_BEFORE(I,J,L,TC0,Raw,CutoffEvaluated,Rate,Freq,RKT)
    INTEGER, INTENT(IN) :: I,J,L
    REAL(fp), INTENT(IN) :: TC0,Raw
    LOGICAL, INTENT(IN) :: CutoffEvaluated
    REAL(fp), OPTIONAL, INTENT(IN) :: Rate,Freq,RKT
    IF (.NOT. Enabled) RETURN
    IF (Masks(I,J,L,1)/=0) ERROR STOP 'Duplicate BRC chemistry before cell'
    Values(I,J,L,1)=TC0; Values(I,J,L,2)=Raw; Masks(I,J,L,1)=1
    IF (PRESENT(Rate) .NEQV. PRESENT(Freq)) ERROR STOP 'Incomplete BRC chemistry rate operands'
    IF (PRESENT(Rate)) THEN
      Values(I,J,L,4)=Rate; Values(I,J,L,5)=Freq; Masks(I,J,L,3)=1
    ENDIF
    IF (PRESENT(RKT)) THEN
      Values(I,J,L,6)=RKT; Masks(I,J,L,4)=1
    ENDIF
    Masks(I,J,L,5)=INT(MERGE(1,0,CutoffEvaluated),int32)
    IF (CutoffEvaluated) Masks(I,J,L,6)=INT(MERGE(1,0,Raw<StoredSmall),int32)
  END SUBROUTINE BRC_CHEM_CELL_BEFORE

  SUBROUTINE BRC_CHEM_CELL_AFTER(I,J,L,Final,Gain,Factor,Incoming)
    INTEGER, INTENT(IN) :: I,J,L
    REAL(fp), INTENT(IN) :: Final,Gain,Factor,Incoming
    IF (.NOT. Enabled) RETURN
    IF (Masks(I,J,L,1)/=1 .OR. Masks(I,J,L,2)/=0) &
         ERROR STOP 'Invalid BRC chemistry after cell'
    Values(I,J,L,3)=Final; Values(I,J,L,7)=Gain
    Values(I,J,L,8)=Factor; Values(I,J,L,9)=Incoming; Masks(I,J,L,2)=1
  END SUBROUTINE BRC_CHEM_CELL_AFTER
END MODULE BRC_CHEM_CAPTURE_MOD
