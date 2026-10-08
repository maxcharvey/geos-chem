! Diagnostic-only native tendency observer. Default off, immutable streams.
MODULE BRC_TEND_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE Input_Opt_Mod, ONLY: OptInput
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE PhysConstants, ONLY: g0_100, AIRMW
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS, GET_TS_DYN
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_TEND_SNAPSHOT, BRC_TEND_SPECIES, BRC_TEND_COLUMN
  PUBLIC :: BRC_TEND_CELL_BEGIN, BRC_TEND_CELL_END
  PUBLIC :: BRC_TEND_DEP_BASE, BRC_TEND_DEP_INPUT, BRC_TEND_DEP_APPLIED
  PUBLIC :: BRC_TEND_EM_INPUT, BRC_TEND_EM_APPLIED
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  LOGICAL, SAVE :: Initialized=.FALSE., Enabled=.FALSE., Early=.FALSE., HasDT=.FALSE.
  INTEGER, SAVE :: Ids(7)=0, OriginalUnits=0, PreviousUnits=-1
  INTEGER, SAVE :: Sequence=0, LastPhase=-1, StoredStep=-1
  REAL(fp), SAVE :: StoredTS=0.0_fp, StoredDT=0.0_fp, StoredMW(7)=0.0_fp
  INTEGER, ALLOCATABLE, SAVE :: Slots(:), Mapping(:)
  INTEGER(int32), SAVE :: SpeciesFlags(7,5)=0
  INTEGER(int32), ALLOCATABLE, SAVE :: Bounds(:,:,:,:), Masks(:,:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Values(:,:,:,:,:)
CONTAINS
  SUBROUTINE BRC_TEND_SNAPSHOT(Phase,Input_Opt,State_Chm,State_Grid,State_Met, &
       OnlyAbovePBL,DT,NativeTS,Previous)
    INTEGER, INTENT(IN) :: Phase
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    LOGICAL, INTENT(IN) :: OnlyAbovePBL
    REAL(fp), OPTIONAL, INTENT(IN) :: DT,NativeTS
    INTEGER, OPTIONAL, INTENT(IN) :: Previous
    INTEGER :: Step,IOS,Length,NX,NY,NZ,S,N,U,F,LIM
    INTEGER(int32) :: H(49)
    REAL(fp) :: MW(7)
    CHARACTER(LEN=1200) :: Path
    Step=BRC_EVENT_CAPTURE_STEP()
    Enabled=.FALSE.
    IF (Step==0 .OR. .NOT. Input_Opt%amIRoot) RETURN
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_TEND_CAPTURE_DIR',Directory,Length,IOS)
      IF (IOS==-1) ERROR STOP 'BRC tendency directory truncated'
      IF (IOS==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
#if defined(TOMAS) || defined(APM) || defined(ADJOINT) || !defined(MODEL_CLASSIC)
        ERROR STOP 'BRC tendency capture unsupported build'
#endif
        IF (STORAGE_SIZE(1.0_fp)/=64 .OR. State_Grid%NestedGrid) &
             ERROR STOP 'BRC tendency capture requires REAL8 global'
        ALLOCATE(Slots(State_Chm%nSpecies),STAT=IOS)
        IF (IOS/=0) ERROR STOP 'BRC tendency slot allocation failed'
        Slots=0
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
          IF (Ids(S)<=0 .OR. Ids(S)>SIZE(Slots)) ERROR STOP 'Missing BRC tendency parent'
          IF (Slots(Ids(S))/=0) ERROR STOP 'Duplicate BRC tendency identity'
          Slots(Ids(S))=S
        ENDDO
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (Phase<0 .OR. Phase>3) ERROR STOP 'Invalid BRC tendency phase'
    IF (Input_Opt%LTURB .AND. .NOT. Input_Opt%LNLPBL) &
         ERROR STOP 'BRC tendency capture excludes subsequent full PBL mixing'
    NX=State_Grid%NX; NY=State_Grid%NY; NZ=State_Grid%NZ
    IF (State_Chm%nAdvect/=SIZE(State_Chm%Map_Advect) .OR. &
        State_Chm%nAdvect<7 .OR. State_Chm%nAdvect>State_Chm%nSpecies) &
         ERROR STOP 'Invalid BRC tendency mapping count'
    DO N=1,State_Chm%nAdvect
      IF (State_Chm%Map_Advect(N)<1 .OR. &
          State_Chm%Map_Advect(N)>State_Chm%nSpecies) ERROR STOP 'Invalid BRC tendency map id'
      IF (COUNT(State_Chm%Map_Advect==State_Chm%Map_Advect(N))/=1) &
           ERROR STOP 'Duplicate BRC tendency mapping'
    ENDDO
    DO S=1,7
      IF (COUNT(State_Chm%Map_Advect==Ids(S))/=1) ERROR STOP 'BRC tendency parent not advected'
      MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
      IF (.NOT. IEEE_IS_FINITE(MW(S)) .OR. MW(S)<=0.0_fp) &
           ERROR STOP 'Invalid BRC tendency molecular weight'
      IF (.NOT. ALL(IEEE_IS_FINITE(State_Chm%Species(Ids(S))%Conc)) .OR. &
          ANY(State_Chm%Species(Ids(S))%Conc<0.0_fp)) ERROR STOP 'Invalid BRC tendency stock'
    ENDDO
    IF (Phase==0) THEN
      IF (PRESENT(NativeTS) .OR. PRESENT(Previous)) ERROR STOP 'Premature BRC tendency scalar read'
      IF (LastPhase/=-1 .AND. LastPhase/=3) ERROR STOP 'Incomplete previous BRC tendency event'
      IF (StoredStep==Step) ERROR STOP 'Duplicate BRC tendency event'
      StoredStep=Step; StoredMW=MW
      Early=.NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%DoEmissions
      HasDT=PRESENT(DT); StoredDT=0.0_fp; StoredTS=0.0_fp; PreviousUnits=-1
      IF (HasDT) THEN
        StoredDT=DT
        IF (.NOT. IEEE_IS_FINITE(DT)) ERROR STOP 'Nonfinite BRC tendency DT'
        IF (.NOT. Early .AND. DT<=0.0_fp) ERROR STOP 'Unsupported BRC tendency DT'
      ENDIF
      OriginalUnits=State_Chm%Species(State_Chm%Map_Advect(1))%Units
      IF (OriginalUnits/=2 .AND. OriginalUnits/=6) ERROR STOP 'Unsupported tendency input units'
      IF (ALLOCATED(Mapping)) DEALLOCATE(Mapping)
      ALLOCATE(Mapping(State_Chm%nAdvect),STAT=IOS)
      IF (IOS/=0) ERROR STOP 'BRC tendency mapping allocation failed'
      Mapping=State_Chm%Map_Advect
    ELSE
      IF (Step/=StoredStep .OR. .NOT. ALLOCATED(Mapping)) &
           ERROR STOP 'BRC tendency unmatched event'
      IF (SIZE(Mapping)/=State_Chm%nAdvect) ERROR STOP 'BRC tendency mapping count changed'
      IF (ANY(Mapping/=State_Chm%Map_Advect) .OR. ANY(MW/=StoredMW)) &
           ERROR STOP 'BRC tendency identities changed'
      IF (PRESENT(DT) .NEQV. HasDT) ERROR STOP 'BRC tendency DT presence changed'
      IF (PRESENT(DT)) THEN
        IF (DT/=StoredDT) ERROR STOP 'BRC tendency DT changed'
      ENDIF
      IF (Early) THEN
        IF (Phase/=3 .OR. LastPhase/=0 .OR. PRESENT(NativeTS) .OR. PRESENT(Previous)) &
             ERROR STOP 'Invalid BRC tendency early return'
      ELSE
        IF (Phase/=LastPhase+1) ERROR STOP 'BRC tendency phase order'
        IF (Phase==1) THEN
          IF (.NOT. PRESENT(NativeTS) .OR. .NOT. PRESENT(Previous)) &
               ERROR STOP 'Missing actual BRC tendency TS/previous'
          StoredTS=NativeTS; PreviousUnits=Previous
          IF (.NOT. IEEE_IS_FINITE(StoredTS) .OR. StoredTS<=0.0_fp .OR. &
              PreviousUnits/=OriginalUnits) ERROR STOP 'Invalid actual tendency TS/previous'
          IF (HasDT) THEN
            IF (StoredTS/=StoredDT) ERROR STOP 'Native BRC tendency optional timestep mismatch'
          ELSE
            IF (StoredTS/=REAL(GET_TS_DYN(),fp)) ERROR STOP 'Native BRC tendency timestep mismatch'
          ENDIF
          IF (.NOT. ALLOCATED(Values)) THEN
            ALLOCATE(Values(NX,NY,NZ,7,12),Masks(NX,NY,NZ,7,8), &
                 Bounds(NX,NY,7,5),STAT=IOS)
            IF (IOS/=0) ERROR STOP 'BRC tendency scratch allocation failed'
          ENDIF
          IF (ANY(SHAPE(Values)/=[NX,NY,NZ,7,12])) ERROR STOP 'BRC tendency grid changed'
          Values=0.0_fp; Masks=0_int32; Bounds=0_int32; SpeciesFlags=0_int32
          SpeciesFlags(:,2)=-1_int32
        ELSE IF (PRESENT(Previous)) THEN
          IF (Previous/=PreviousUnits) ERROR STOP 'BRC tendency previous units changed'
        ENDIF
        IF (PRESENT(NativeTS)) THEN
          IF (NativeTS/=StoredTS) ERROR STOP 'BRC tendency native TS changed'
        ENDIF
      ENDIF
    ENDIF
    DO S=1,7
      N=State_Chm%Species(Ids(S))%Units
      IF ((Phase==0 .OR. Phase==3) .AND. N/=OriginalUnits) &
           ERROR STOP 'BRC tendency entry/restore units mismatch'
      IF ((Phase==1 .OR. Phase==2) .AND. N/=4) ERROR STOP 'BRC tendency body units mismatch'
    ENDDO
    IF (.NOT. ALL(IEEE_IS_FINITE(State_Met%AD)) .OR. ANY(State_Met%AD<=0.0_fp) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%AIRVOL)) .OR. ANY(State_Met%AIRVOL<=0.0_fp) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%DELP_DRY)) .OR. ANY(State_Met%DELP_DRY<=0.0_fp) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Met%PBL_TOP_L)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(State_Grid%Area_M2)) .OR. ANY(State_Grid%Area_M2<=0.0_fp)) &
         ERROR STOP 'Invalid BRC tendency geometry'
    IF (Phase==2) THEN
      IF (ANY(SpeciesFlags(:,1)/=1) .OR. ANY(Masks(:,:,:,:,1)/=Masks(:,:,:,:,2))) &
           ERROR STOP 'Incomplete BRC tendency body coverage'
      IF (.NOT. ALL(IEEE_IS_FINITE(Values))) ERROR STOP 'Nonfinite BRC tendency operand'
    ENDIF
    H(1:13)=[1_int32,8_int32,INT(Z'01020304',int32),INT(NX,int32),INT(NY,int32), &
         INT(NZ,int32),INT(Sequence+1,int32),INT(Step,int32),INT(GET_NYMD(),int32), &
         INT(GET_NHMS(),int32),INT(Phase,int32),INT(State_Chm%nAdvect,int32), &
         INT(State_Chm%nSpecies,int32)]
#ifdef MODEL_CLASSIC
    LIM=MERGE(1,0,Input_Opt%LIMGRID)
#else
    LIM=0
#endif
    H(14:28)=INT([MERGE(1,0,OnlyAbovePBL),MERGE(1,0,HasDT), &
         MERGE(1,0,.NOT. Early .AND. Phase>0),MERGE(1,0,PreviousUnits/=-1), &
         PreviousUnits,OriginalUnits,MERGE(1,0,Input_Opt%LTURB), &
         MERGE(1,0,Input_Opt%LNLPBL),MERGE(1,0,Input_Opt%DoEmissions), &
         MERGE(1,0,Input_Opt%LDRYD),MERGE(1,0,Input_Opt%PBL_DRYDEP), &
         MERGE(1,0,Input_Opt%LINEAR_CHEM),LIM,MERGE(1,0,Early),MERGE(1,0,Phase==2)],int32)
    H(29:35)=INT(Ids,int32); H(43:49)=1_int32
    DO S=1,7
      H(35+S)=INT(State_Chm%Species(Ids(S))%Units,int32)
    ENDDO
    WRITE(Path,'(a,"/tendency_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence+1,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),ACCESS='stream',FORM='unformatted',STATUS='new', &
         ACTION='write',IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'Cannot create immutable BRC tendency capture'
    WRITE(U,IOSTAT=IOS) 'BRCTD001',H,Names,MW,g0_100,AIRMW,StoredTS,StoredDT, &
         INT(Mapping,int32),REAL(State_Met%AD,fp),REAL(State_Met%AIRVOL,fp), &
         REAL(State_Met%DELP_DRY,fp),REAL(State_Met%PBL_TOP_L,fp), &
         REAL(State_Grid%Area_M2,fp),INT(State_Met%ChemGridLev,int32)
    IF (IOS/=0) ERROR STOP 'Cannot write BRC tendency header/geometry'
    DO S=1,7
      WRITE(U,IOSTAT=IOS) State_Chm%Species(Ids(S))%Conc
      IF (IOS/=0) ERROR STOP 'Cannot write BRC tendency stock'
    ENDDO
    IF (Phase==2) THEN
      WRITE(U,IOSTAT=IOS) SpeciesFlags,Bounds
      IF (IOS/=0) ERROR STOP 'Cannot write BRC tendency flags/bounds'
      DO F=1,12
        WRITE(U,IOSTAT=IOS) Values(:,:,:,:,F)
        IF (IOS/=0) ERROR STOP 'Cannot write BRC tendency native operand'
      ENDDO
      DO F=1,8
        WRITE(U,IOSTAT=IOS) Masks(:,:,:,:,F)
        IF (IOS/=0) ERROR STOP 'Cannot write BRC tendency native mask'
      ENDDO
    ENDIF
    CLOSE(U,IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'Cannot close BRC tendency capture'
    Sequence=Sequence+1; LastPhase=Phase; Enabled=Phase==1
  END SUBROUTINE BRC_TEND_SNAPSHOT

  INTEGER FUNCTION FindSlot(N) RESULT(S)
    INTEGER, INTENT(IN) :: N
    S=0
    IF (.NOT. Enabled) RETURN
    IF (N<1 .OR. N>SIZE(Slots)) ERROR STOP 'Invalid native tendency species id'
    S=Slots(N)
  END FUNCTION FindSlot

  SUBROUTINE BRC_TEND_SPECIES(N,DryDepId,DryDepSpec,EmisSpec,ChemGridOnly)
    INTEGER, INTENT(IN) :: N,DryDepId
    LOGICAL, INTENT(IN) :: DryDepSpec,EmisSpec,ChemGridOnly
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (SpeciesFlags(S,1)/=0) ERROR STOP 'Duplicate native tendency species flags'
    SpeciesFlags(S,:)=INT([1,DryDepId,MERGE(1,0,DryDepSpec), &
         MERGE(1,0,EmisSpec),MERGE(1,0,ChemGridOnly)],int32)
  END SUBROUTINE BRC_TEND_SPECIES

  SUBROUTINE BRC_TEND_COLUMN(N,I,J,PBL_TOP,L1,DRYD_TOP,EMIS_TOP,L2)
    INTEGER, INTENT(IN) :: N,I,J,PBL_TOP,L1,DRYD_TOP,EMIS_TOP,L2
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (SpeciesFlags(S,1)/=1 .OR. ALL(SpeciesFlags(S,3:4)==0)) &
         ERROR STOP 'Invalid native tendency active column'
    IF (ANY(Bounds(I,J,S,:)/=0)) ERROR STOP 'Duplicate native tendency column'
    Bounds(I,J,S,:)=INT([PBL_TOP,L1,DRYD_TOP,EMIS_TOP,L2],int32)
  END SUBROUTINE BRC_TEND_COLUMN

  SUBROUTINE BRC_TEND_CELL_BEGIN(N,I,J,L,Q)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: Q
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (SpeciesFlags(S,1)/=1 .OR. Masks(I,J,L,S,1)/=0) &
         ERROR STOP 'Invalid native tendency cell entry'
    IF (L<Bounds(I,J,S,2) .OR. L>Bounds(I,J,S,5)) &
         ERROR STOP 'Native tendency cell outside evaluated bounds'
    Masks(I,J,L,S,1)=1_int32
    Values(I,J,L,S,1:3)=Q
  END SUBROUTINE BRC_TEND_CELL_BEGIN

  SUBROUTINE BRC_TEND_DEP_BASE(N,I,J,L,FRQ)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: FRQ
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,1)/=1 .OR. Masks(I,J,L,S,3)/=0) &
         ERROR STOP 'Invalid native tendency deposition base'
    Values(I,J,L,S,4)=FRQ
  END SUBROUTINE BRC_TEND_DEP_BASE

  SUBROUTINE BRC_TEND_DEP_INPUT(N,I,J,L,FRQ,PNOXLOSS,Found,TMP)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: FRQ,PNOXLOSS
    LOGICAL, INTENT(IN) :: Found
    REAL(fp), OPTIONAL, INTENT(IN) :: TMP
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,1)/=1 .OR. Masks(I,J,L,S,3)/=0) &
         ERROR STOP 'Invalid native tendency deposition input'
    IF (PRESENT(TMP) .NEQV. Found) ERROR STOP 'Native tendency deposition TMP presence'
    IF (PNOXLOSS/=0.0_fp) ERROR STOP 'Native tendency selected PNOXLOSS unsupported'
    Masks(I,J,L,S,3)=1_int32; Masks(I,J,L,S,4)=INT(MERGE(1,0,Found),int32)
    Values(I,J,L,S,6)=FRQ; Values(I,J,L,S,7)=PNOXLOSS
    IF (PRESENT(TMP)) Values(I,J,L,S,5)=TMP
  END SUBROUTINE BRC_TEND_DEP_INPUT

  SUBROUTINE BRC_TEND_DEP_APPLIED(N,I,J,L,RKT,FRAC,FLUX,Q)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: RKT,FRAC,FLUX,Q
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,3)/=1 .OR. Masks(I,J,L,S,5)/=0) &
         ERROR STOP 'Invalid native tendency deposition application'
    Masks(I,J,L,S,5)=1_int32
    Values(I,J,L,S,2)=Q
    Values(I,J,L,S,8)=RKT; Values(I,J,L,S,9)=FRAC; Values(I,J,L,S,10)=FLUX
  END SUBROUTINE BRC_TEND_DEP_APPLIED

  SUBROUTINE BRC_TEND_EM_INPUT(N,I,J,L,Found,TMP)
    INTEGER, INTENT(IN) :: N,I,J,L
    LOGICAL, INTENT(IN) :: Found
    REAL(fp), OPTIONAL, INTENT(IN) :: TMP
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,1)/=1 .OR. Masks(I,J,L,S,6)/=0) &
         ERROR STOP 'Invalid native tendency emission input'
    IF (PRESENT(TMP) .NEQV. Found) ERROR STOP 'Native tendency emission TMP presence'
    Masks(I,J,L,S,6)=1_int32; Masks(I,J,L,S,7)=INT(MERGE(1,0,Found),int32)
    IF (PRESENT(TMP)) Values(I,J,L,S,11)=TMP
  END SUBROUTINE BRC_TEND_EM_INPUT

  SUBROUTINE BRC_TEND_EM_APPLIED(N,I,J,L,FLUX)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: FLUX
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,7)/=1 .OR. Masks(I,J,L,S,8)/=0) &
         ERROR STOP 'Invalid native tendency emission application'
    Masks(I,J,L,S,8)=1_int32; Values(I,J,L,S,12)=FLUX
  END SUBROUTINE BRC_TEND_EM_APPLIED

  SUBROUTINE BRC_TEND_CELL_END(N,I,J,L,Q)
    INTEGER, INTENT(IN) :: N,I,J,L
    REAL(fp), INTENT(IN) :: Q
    INTEGER :: S
    S=FindSlot(N); IF (S==0) RETURN
    IF (Masks(I,J,L,S,1)/=1 .OR. Masks(I,J,L,S,2)/=0) &
         ERROR STOP 'Invalid native tendency cell exit'
    Masks(I,J,L,S,2)=1_int32; Values(I,J,L,S,3)=Q
  END SUBROUTINE BRC_TEND_CELL_END
END MODULE BRC_TEND_CAPTURE_MOD
