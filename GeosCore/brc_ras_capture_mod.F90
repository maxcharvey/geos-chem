! Default-off, selected-column native RAS observer. No physical writes.
MODULE BRC_RAS_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE Input_Opt_Mod, ONLY: OptInput
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE PhysConstants, ONLY: g0_100, AIRMW
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_RAS_SNAPSHOT, BRC_RAS_COLUMN, BRC_RAS_TRACE, BRC_RAS_ACTIVE
  INTEGER, PARAMETER :: NZ=47, NCOL=2, MAXROWS=32768, MAXVALUES=20
  INTEGER, PARAMETER :: Columns(2,2)=RESHAPE([46,19,131,28],[2,2])
  CHARACTER(LEN=8), PARAMETER :: Names(7)=[CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  TYPE ColumnBuffer
    INTEGER(int32) :: Meta(8)=0
    REAL(fp) :: Scalars(5)=0.0_fp, BMASS(NZ)=0.0_fp, PDOWN(NZ)=0.0_fp
    INTEGER(int32), ALLOCATABLE :: Rows(:,:)
    REAL(fp), ALLOCATABLE :: Q(:,:), V(:,:)
  END TYPE ColumnBuffer
  TYPE(ColumnBuffer), SAVE :: Buffer(NCOL)
  LOGICAL, SAVE :: Initialized=.FALSE., Enabled=.FALSE.
  CHARACTER(LEN=1024), SAVE :: Directory=''
  INTEGER, SAVE :: Ids(7)=0, Positions(7)=0, Sequence=0, LastPhase=-1, StoredStep=0
  INTEGER, ALLOCATABLE, SAVE :: Slots(:)
CONTAINS
  SUBROUTINE BRC_RAS_SNAPSHOT(Phase,Input_Opt,State_Chm,State_Grid,State_Met,FSOL,DT)
    INTEGER, INTENT(IN) :: Phase
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    REAL(fp), INTENT(IN) :: FSOL(:,:,:,:),DT
    INTEGER :: Step,IOS,Length,S,N,C,I,J,U,R
    INTEGER(int32) :: H(48),WetIds(7)
    REAL(fp) :: MW(7),Geometry(NZ,18,NCOL),Stock(NZ,7,NCOL),Area(NCOL)
    CHARACTER(LEN=1200) :: Path
    Step=BRC_EVENT_CAPTURE_STEP()
    IF (Step==0 .OR. .NOT.Input_Opt%amIRoot) THEN
      Enabled=.FALSE.
      RETURN
    ENDIF
    IF (.NOT.Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_RAS_CAPTURE_DIR',Directory,Length,IOS)
      IF (IOS==-1) ERROR STOP 'BRC RAS directory truncated'
      IF (IOS==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
#if defined(LUO_WETDEP) || defined(TOMAS) || defined(APM) || defined(ADJOINT) || !defined(MODEL_CLASSIC)
        ERROR STOP 'BRC RAS observer unsupported build'
#endif
        IF (STORAGE_SIZE(1.0_fp)/=64) ERROR STOP 'BRC RAS requires REAL8'
        ALLOCATE(Slots(State_Chm%nSpecies),STAT=IOS)
        IF (IOS/=0) ERROR STOP 'BRC RAS slot allocation'
        Slots=0
        DO S=1,7
          Ids(S)=Ind_(TRIM(Names(S)))
          IF (Ids(S)<=0 .OR. Ids(S)>SIZE(Slots)) ERROR STOP 'BRC RAS missing parent'
          IF (Slots(Ids(S))/=0) ERROR STOP 'BRC RAS duplicate parent'
          Slots(Ids(S))=S
        ENDDO
        DO C=1,NCOL
          ALLOCATE(Buffer(C)%Rows(8,MAXROWS),Buffer(C)%Q(NZ,MAXROWS), &
               Buffer(C)%V(MAXVALUES,MAXROWS),STAT=IOS)
          IF (IOS/=0) ERROR STOP 'BRC RAS frame allocation'
        ENDDO
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (Phase<0 .OR. Phase>1) ERROR STOP 'BRC RAS phase'
    IF (State_Grid%NestedGrid .OR. State_Grid%NX/=144 .OR. State_Grid%NY/=91 .OR. &
        State_Grid%NZ/=NZ .OR. Input_Opt%Grell_Freitas_Convection .OR. &
        Input_Opt%ITS_A_MERCURY_SIM) ERROR STOP 'BRC RAS runtime mode'
    IF (State_Chm%nAdvect/=SIZE(State_Chm%Map_Advect) .OR. &
        SIZE(FSOL,4)/=State_Chm%nAdvect .OR. ANY(SHAPE(FSOL(:,:,:,1))/=[144,91,NZ])) &
         ERROR STOP 'BRC RAS mapping shape'
    IF (Phase==0) THEN
      IF (LastPhase/=-1 .AND. LastPhase/=1) ERROR STOP 'BRC RAS incomplete previous call'
      StoredStep=Step
      DO C=1,NCOL
        Buffer(C)%Meta=0;Buffer(C)%Scalars=0.0_fp
        Buffer(C)%BMASS=0.0_fp;Buffer(C)%PDOWN=0.0_fp
      ENDDO
    ELSE
      IF (LastPhase/=0 .OR. StoredStep/=Step) ERROR STOP 'BRC RAS unmatched exit'
      IF (ANY([(Buffer(C)%Meta(1)/=1,C=1,NCOL)])) ERROR STOP 'BRC RAS incomplete columns'
    ENDIF
    DO N=1,State_Chm%nAdvect
      I=State_Chm%Map_Advect(N)
      IF (I<=0 .OR. I>SIZE(Slots)) ERROR STOP 'BRC RAS map range'
      IF (COUNT(State_Chm%Map_Advect==I)/=1) ERROR STOP 'BRC RAS map duplicate'
    ENDDO
    DO S=1,7
      IF (COUNT(State_Chm%Map_Advect==Ids(S))/=1) ERROR STOP 'BRC RAS parent map'
      Positions(S)=FINDLOC(State_Chm%Map_Advect,Ids(S),DIM=1)
      IF (State_Chm%Species(Ids(S))%Units/=2) ERROR STOP 'BRC RAS actual units'
      MW(S)=State_Chm%SpcData(Ids(S))%Info%MW_g
      WetIds(S)=INT(State_Chm%SpcData(Ids(S))%Info%WetDepId,int32)
    ENDDO
    IF (.NOT.ALL(IEEE_IS_FINITE(MW)) .OR. ANY(MW<=0.0_fp)) ERROR STOP 'BRC RAS molecular weights'
    DO C=1,NCOL
      I=Columns(1,C);J=Columns(2,C)
      Area(C)=State_Grid%Area_M2(I,J)
      Geometry(:,1,C)=State_Met%AD(I,J,:)
      Geometry(:,2,C)=State_Met%DELP_DRY(I,J,:)
      Geometry(:,3,C)=State_Met%CMFMC(I,J,2:NZ+1)
      Geometry(:,4,C)=State_Met%DTRAIN(I,J,:)
      Geometry(:,5,C)=State_Met%DQRCU(I,J,:)
      Geometry(:,6,C)=State_Met%PFICU(I,J,1:NZ)
      Geometry(:,7,C)=State_Met%PFLCU(I,J,1:NZ)
      Geometry(:,8,C)=State_Met%REEVAPCN(I,J,:)
      Geometry(:,9,C)=State_Met%T(I,J,:)
      Geometry(:,10,C)=State_Met%BXHEIGHT(I,J,:)
      Geometry(:,18,C)=State_Met%AIRVOL(I,J,:)
      DO S=1,7
        Geometry(:,10+S,C)=FSOL(I,J,:,Positions(S))
        Stock(:,S,C)=State_Chm%Species(Ids(S))%Conc(I,J,:)
      ENDDO
    ENDDO
    IF (.NOT.ALL(IEEE_IS_FINITE(Geometry)) .OR. .NOT.ALL(IEEE_IS_FINITE(Stock)) .OR. &
        .NOT.ALL(IEEE_IS_FINITE(Area)) .OR. ANY(Stock<0.0_fp) .OR. ANY(Area<=0.0_fp) .OR. &
        ANY(Geometry(:,1:2,:)<=0.0_fp) .OR. ANY(Geometry(:,18,:)<=0.0_fp) .OR. &
        ANY(Geometry(:,11:17,:)<0.0_fp) .OR. ANY(Geometry(:,11:17,:)>1.0_fp) .OR. &
        .NOT.IEEE_IS_FINITE(DT) .OR. DT<=0.0_fp) ERROR STOP 'BRC RAS nonphysical input'
    H=0
    H(1:16)=INT([1,8,INT(Z'01020304'),144,91,NZ,Sequence+1,Step,GET_NYMD(), &
         GET_NHMS(),Phase,State_Chm%nAdvect,State_Chm%nSpecies,NCOL,MAXROWS,32],int32)
    H(22)=INT(MERGE(1,0,Input_Opt%LIMGRID),int32)
#ifdef DEBUG
    H(24)=1_int32
#endif
    H(25:31)=INT(Ids,int32);H(32:38)=2_int32;H(39:45)=INT(Positions,int32)
    H(46:48)=INT([NZ,MAXVALUES,568],int32)
    WRITE(Path,'(a,"/ras_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence+1,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),STATUS='new',ACCESS='stream',FORM='unformatted', &
         ACTION='write',IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'BRC RAS immutable output open'
    WRITE(U,IOSTAT=IOS) 'BRCRS001',H,Names,MW,g0_100,AIRMW,DT, &
         INT(State_Chm%Map_Advect,int32),INT(Columns,int32),WetIds,Area,Geometry,Stock
    IF (IOS/=0) ERROR STOP 'BRC RAS output input frame'
    DO C=1,NCOL
      WRITE(U,IOSTAT=IOS) Buffer(C)%Meta,Buffer(C)%Scalars,Buffer(C)%BMASS,Buffer(C)%PDOWN
      IF (IOS/=0) ERROR STOP 'BRC RAS output column metadata'
      DO R=1,Buffer(C)%Meta(8)
        WRITE(U,IOSTAT=IOS) Buffer(C)%Rows(:,R),Buffer(C)%Q(:,R),Buffer(C)%V(:,R)
        IF (IOS/=0) ERROR STOP 'BRC RAS output trace frame'
      ENDDO
    ENDDO
    CLOSE(U,IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'BRC RAS close output'
    Sequence=Sequence+1;LastPhase=Phase;Enabled=Phase==0
  END SUBROUTINE BRC_RAS_SNAPSHOT

  INTEGER FUNCTION ColumnSlot(I,J) RESULT(C)
    INTEGER, INTENT(IN) :: I,J
    C=0
    IF (.NOT.Enabled) RETURN
    IF (I==46 .AND. J==19) C=1
    IF (I==131 .AND. J==28) C=2
  END FUNCTION ColumnSlot

  LOGICAL FUNCTION BRC_RAS_ACTIVE(I,J,IC) RESULT(Active)
    INTEGER, INTENT(IN) :: I,J,IC
    Active=.FALSE.
    IF (ColumnSlot(I,J)==0) RETURN
    IF (IC<=0 .OR. IC>SIZE(Slots)) ERROR STOP 'BRC RAS active id'
    Active=Slots(IC)>0
  END FUNCTION BRC_RAS_ACTIVE

  SUBROUTINE BRC_RAS_COLUMN(I,J,CLDBASE,KTOP,NS,NDT,DNS,SDT,MB,TS_DYN,TINYNUM,BMASS,PDOWN)
    INTEGER, INTENT(IN) :: I,J,CLDBASE,KTOP,NS,NDT
    REAL(fp), INTENT(IN) :: DNS,SDT,MB,TS_DYN,TINYNUM,BMASS(:),PDOWN(:)
    INTEGER :: C
    C=ColumnSlot(I,J)
    IF (C==0) RETURN
    IF (Buffer(C)%Meta(1)/=0 .OR. NS<1 .OR. NS>32 .OR. SIZE(BMASS)/=NZ .OR. &
        SIZE(PDOWN)/=NZ .OR. CLDBASE<1 .OR. CLDBASE>NZ .OR. KTOP/=NZ-1) &
         ERROR STOP 'BRC RAS native column bounds'
    Buffer(C)%Meta=INT([1,I,J,CLDBASE,KTOP,NS,NDT,0],int32)
    Buffer(C)%Scalars=[DNS,SDT,MB,TS_DYN,TINYNUM]
    Buffer(C)%BMASS=BMASS;Buffer(C)%PDOWN=PDOWN
    IF (.NOT.ALL(IEEE_IS_FINITE(Buffer(C)%Scalars)) .OR. &
        .NOT.ALL(IEEE_IS_FINITE(BMASS)) .OR. .NOT.ALL(IEEE_IS_FINITE(PDOWN))) &
         ERROR STOP 'BRC RAS native column finite'
  END SUBROUTINE BRC_RAS_COLUMN

  SUBROUTINE BRC_RAS_TRACE(Stage,I,J,IC,NA,NW,ISTEP,K,Q,Values,Aux)
    INTEGER, INTENT(IN) :: Stage,I,J,IC,NA,NW,ISTEP,K,Aux
    REAL(fp), INTENT(IN) :: Q(:),Values(:)
    INTEGER :: C,S,R,NV
    C=ColumnSlot(I,J)
    IF (C==0) RETURN
    IF (IC<=0 .OR. IC>SIZE(Slots)) ERROR STOP 'BRC RAS trace id'
    S=Slots(IC)
    IF (S==0) RETURN
    IF (NA/=Positions(S) .OR. Buffer(C)%Meta(1)/=1 .OR. SIZE(Q)/=NZ) &
         ERROR STOP 'BRC RAS trace binding'
    NV=SIZE(Values)
    IF (NV>MAXVALUES) ERROR STOP 'BRC RAS trace values count'
    R=Buffer(C)%Meta(8)+1
    IF (R>MAXROWS) ERROR STOP 'BRC RAS trace overflow'
    IF (.NOT.ALL(IEEE_IS_FINITE(Q)) .OR. .NOT.ALL(IEEE_IS_FINITE(Values))) &
         ERROR STOP 'BRC RAS trace nonfinite'
    Buffer(C)%Rows(:,R)=INT([Stage,IC,NA,NW,ISTEP,K,NV,Aux],int32)
    Buffer(C)%Q(:,R)=Q;Buffer(C)%V(:,R)=0.0_fp
    IF (NV>0) Buffer(C)%V(1:NV,R)=Values
    Buffer(C)%Meta(8)=R
  END SUBROUTINE BRC_RAS_TRACE
END MODULE BRC_RAS_CAPTURE_MOD
