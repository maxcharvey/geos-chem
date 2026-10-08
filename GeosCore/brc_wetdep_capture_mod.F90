! Default-off native large-scale wetdep observer. No physical state writes.
MODULE BRC_WETDEP_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE Input_Opt_Mod, ONLY: OptInput
  USE State_Chm_Mod, ONLY: ChmState,Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE PhysConstants, ONLY: AIRMW,g0_100
  USE Time_Mod, ONLY: GET_NYMD,GET_NHMS
  USE BRC_EVENT_CAPTURE_MOD, ONLY: BRC_EVENT_CAPTURE_STEP
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_WET_SNAPSHOT,BRC_WET_COLUMN,BRC_WET_ACTIVE,BRC_WET_TRACE,BRC_WET_SAFETY
  INTEGER, PARAMETER :: NZ=47,NC=2,NF=7,NV=28,MAXROWS=8192
  INTEGER, PARAMETER :: Columns(2,2)=RESHAPE([46,19,131,28],[2,2])
  CHARACTER(8), PARAMETER :: Names(NF)=[CHARACTER(8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  TYPE ColumnBuffer
    INTEGER(int32) :: Meta(4)=0
    INTEGER(int32), ALLOCATABLE :: Rows(:,:),Flags(:,:),Negative(:,:)
    REAL(fp), ALLOCATABLE :: Q(:,:,:),Rain(:,:,:),V(:,:)
  END TYPE
  TYPE(ColumnBuffer), SAVE :: Buffer(NC)
  LOGICAL, SAVE :: Initialized=.FALSE.,Enabled=.FALSE.
  CHARACTER(1024), SAVE :: Directory=''
  INTEGER, SAVE :: Ids(NF)=0,Positions(NF)=0,Sequence=0,LastPhase=-1,StoredStep=0
  INTEGER, ALLOCATABLE, SAVE :: Slots(:)
CONTAINS
  SUBROUTINE BRC_WET_SNAPSHOT(Phase,LS,DT,Previous,SO2Id,Input_Opt,State_Chm,State_Grid,State_Met,DSpc)
    INTEGER, INTENT(IN) :: Phase,Previous,SO2Id
    LOGICAL, INTENT(IN) :: LS
    REAL(fp), INTENT(IN) :: DT
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    REAL(fp), OPTIONAL, INTENT(IN) :: DSpc(:,:,:,:)
    INTEGER :: Step,IOS,Length,S,N,C,I,J,U,R
    INTEGER(int32) :: H(64),Props(10,NF),WetIds(NF)
    REAL(fp) :: MW(NF),Eff(3,NF),Geometry(NZ,13,NC),Stock(NZ,NF,NC),Rain(NZ,NF,NC),Area(NC)
    CHARACTER(1200) :: Path
    Step=BRC_EVENT_CAPTURE_STEP()
    IF (Step==0 .OR. .NOT.Input_Opt%amIRoot) THEN
      Enabled=.FALSE.
      RETURN
    END IF
    IF (.NOT.Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_WETDEP_CAPTURE_DIR',Directory,Length,IOS)
      IF (IOS==-1) ERROR STOP 'BRC wet directory truncated'
      IF (IOS==1 .OR. Length==0) Directory=''
      IF (LEN_TRIM(Directory)>0) THEN
#if defined(LUO_WETDEP) || defined(TOMAS) || defined(APM) || defined(ADJOINT) || !defined(MODEL_CLASSIC)
        ERROR STOP 'BRC wet unsupported build'
#endif
        IF (STORAGE_SIZE(1.0_fp)/=64) ERROR STOP 'BRC wet requires REAL8'
        ALLOCATE(Slots(State_Chm%nSpecies),STAT=IOS)
        IF (IOS/=0) ERROR STOP 'BRC wet slot allocation'
        Slots=0
        DO S=1,NF
          Ids(S)=Ind_(TRIM(Names(S)))
          IF (Ids(S)<=0 .OR. Ids(S)>SIZE(Slots)) ERROR STOP 'BRC wet missing parent'
          IF (Slots(Ids(S))/=0) ERROR STOP 'BRC wet duplicate parent'
          Slots(Ids(S))=S
        END DO
        DO C=1,NC
          ALLOCATE(Buffer(C)%Rows(8,MAXROWS),Buffer(C)%Flags(8,MAXROWS), &
               Buffer(C)%Negative(NZ,MAXROWS),Buffer(C)%Q(NZ,NF,MAXROWS), &
               Buffer(C)%Rain(NZ,NF,MAXROWS),Buffer(C)%V(NV,MAXROWS),STAT=IOS)
          IF (IOS/=0) ERROR STOP 'BRC wet frame allocation'
        END DO
      END IF
    END IF
    IF (LEN_TRIM(Directory)==0) RETURN
    IF (State_Grid%NestedGrid .OR. ANY([State_Grid%NX,State_Grid%NY,State_Grid%NZ]/=[144,91,NZ]) &
        .OR. .NOT.LS .OR. .NOT.Input_Opt%ITS_A_FULLCHEM_SIM .OR. Input_Opt%ITS_A_MERCURY_SIM) &
         ERROR STOP 'BRC wet unsupported runtime mode'
    IF (Phase<0 .OR. Phase>3) ERROR STOP 'BRC wet snapshot phase'
    IF (State_Chm%nWetDep/=SIZE(State_Chm%Map_WetDep)) ERROR STOP 'BRC wet mapping size'
    DO N=1,State_Chm%nWetDep
      I=State_Chm%Map_WetDep(N)
      IF (I<=0 .OR. I>SIZE(Slots)) ERROR STOP 'BRC wet mapping range'
      IF (COUNT(State_Chm%Map_WetDep==I)/=1) ERROR STOP 'BRC wet duplicate mapping'
    END DO
    IF (Phase==0) THEN
      IF (LastPhase/=-1 .AND. LastPhase/=3) ERROR STOP 'BRC wet incomplete previous call'
      StoredStep=Step
      DO C=1,NC
        Buffer(C)%Meta=0
      END DO
    ELSE
      IF (LastPhase/=Phase-1 .OR. StoredStep/=Step) ERROR STOP 'BRC wet unmatched snapshot'
    END IF
    IF (PRESENT(DSpc) .NEQV. (Phase>=2)) ERROR STOP 'BRC wet reservoir validity'
    IF (PRESENT(DSpc)) THEN
      IF (ANY(SHAPE(DSpc)/=[State_Chm%nWetDep,NZ,144,91])) ERROR STOP 'BRC wet reservoir shape'
      IF (ANY([(Buffer(C)%Meta(1)/=1,C=1,NC)])) ERROR STOP 'BRC wet incomplete columns'
    END IF
    Props=0;Eff=0.0_fp;Rain=0.0_fp
    DO S=1,NF
      N=Ids(S)
      Positions(S)=FINDLOC(State_Chm%Map_WetDep,N,DIM=1)
      IF (S==1 .AND. Positions(S)/=0) ERROR STOP 'BRC wet FSOAP unexpectedly mapped'
      IF (S>1 .AND. Positions(S)==0) ERROR STOP 'BRC wet missing mapped aerosol'
      IF (N==SO2Id) ERROR STOP 'BRC wet unsupported SO2 coupling'
      Props(:,S)=INT([MERGE(1,0,State_Chm%SpcData(N)%Info%Is_WetDep), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%Is_Gas), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%WD_Is_HNO3), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%WD_Is_H2SO4), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%WD_Is_SO2), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%WD_Is_DSTbin), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%WD_CoarseAer), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%MP_SizeResAer), &
           MERGE(1,0,State_Chm%SpcData(N)%Info%MP_SizeResNum), &
           State_Chm%SpcData(N)%Info%ModelID],int32)
      MW(S)=State_Chm%SpcData(N)%Info%MW_g
      Eff(:,S)=State_Chm%SpcData(N)%Info%WD_RainoutEff
      WetIds(S)=INT(State_Chm%SpcData(N)%Info%WetDepId,int32)
      IF (S>1 .AND. (ANY(Props(2:9,S)/=0) .OR. Props(1,S)/=1)) &
         ERROR STOP 'BRC wet unsupported aerosol properties'
      IF (Phase==0 .OR. Phase==3) THEN
        IF (State_Chm%Species(N)%Units/=2) ERROR STOP 'BRC wet entry/exit units'
      ELSE
        IF (State_Chm%Species(N)%Units/=MERGE(2,4,S==1)) ERROR STOP 'BRC wet body mixed units'
      END IF
    END DO
    Geometry=0.0_fp
    DO C=1,NC
      I=Columns(1,C);J=Columns(2,C)
      Area(C)=State_Grid%Area_M2(I,J)
      Geometry(:,1,C)=State_Met%AD(I,J,:)
      Geometry(:,2,C)=State_Met%AIRVOL(I,J,:)
      Geometry(:,3,C)=State_Met%T(I,J,:)
      Geometry(:,4,C)=State_Met%BXHEIGHT(I,J,:)
      Geometry(:,5,C)=State_Met%QQ(:,I,J)
      Geometry(:,6,C)=State_Met%PDOWN(:,I,J)
      Geometry(:,7,C)=State_Met%REEVAP(:,I,J)
      Geometry(:,8,C)=State_Met%C_H2O(I,J,:)
      Geometry(:,9,C)=State_Met%CLDICE(I,J,:)
      Geometry(:,10,C)=State_Met%CLDLIQ(I,J,:)
      IF (ASSOCIATED(State_Met%CNV_FRC)) Geometry(:,11,C)=State_Met%CNV_FRC(I,J)
      Geometry(:,12,C)=Input_Opt%WETD_CONV_SCAL
      Geometry(:,13,C)=State_Met%DELP_DRY(I,J,:)
      DO S=1,NF
        Stock(:,S,C)=State_Chm%Species(Ids(S))%Conc(I,J,:)
        IF (PRESENT(DSpc)) THEN
          IF (Positions(S)>0) Rain(:,S,C)=DSpc(Positions(S),:,I,J)
        END IF
      END DO
    END DO
    H=0
    H(1:18)=INT([2,8,INT(Z'01020304'),144,91,NZ,Sequence+1,Step,GET_NYMD(), &
         GET_NHMS(),Phase,State_Chm%nWetDep,State_Chm%nSpecies,NC,MAXROWS,NV,MERGE(1,0,LS),Previous],int32)
#ifdef DEBUG
    H(22)=1
#endif
    H(23)=INT(MERGE(1,0,Input_Opt%LIMGRID),int32);H(24)=1
    H(25:31)=INT(Ids,int32)
    H(32:38)=INT([(State_Chm%Species(Ids(S))%Units,S=1,NF)],int32)
    H(39:45)=INT(Positions,int32);H(46:52)=WetIds
    H(53:59)=INT([13,NV,5740,MERGE(1,0,ASSOCIATED(State_Met%CNV_FRC)), &
         MERGE(1,0,PRESENT(DSpc)),SO2Id,NF],int32)
    H(60)=13_int32
    WRITE(Path,'(a,"/wet_",i6.6,"_step",i6.6,".bin")') TRIM(Directory),Sequence+1,Step
    OPEN(NEWUNIT=U,FILE=TRIM(Path),STATUS='new',ACCESS='stream',FORM='unformatted',ACTION='write',IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'BRC wet immutable output open'
    WRITE(U,IOSTAT=IOS) 'BRCWD002',H,Names,MW,AIRMW,g0_100,DT,Input_Opt%WETD_CONV_SCAL, &
         Area,Eff,Props,INT(State_Chm%Map_WetDep,int32),INT(Columns,int32),Geometry,Stock,Rain
    IF (IOS/=0) ERROR STOP 'BRC wet snapshot write'
    DO C=1,NC
      WRITE(U,IOSTAT=IOS) Buffer(C)%Meta
      IF (IOS/=0) ERROR STOP 'BRC wet column metadata write'
      DO R=1,Buffer(C)%Meta(4)
        WRITE(U,IOSTAT=IOS) Buffer(C)%Rows(:,R),Buffer(C)%Flags(:,R),Buffer(C)%Negative(:,R), &
             Buffer(C)%Q(:,:,R),Buffer(C)%Rain(:,:,R),Buffer(C)%V(:,R)
        IF (IOS/=0) ERROR STOP 'BRC wet trace write'
      END DO
    END DO
    CLOSE(U,IOSTAT=IOS)
    IF (IOS/=0) ERROR STOP 'BRC wet output close'
    Sequence=Sequence+1;LastPhase=Phase;Enabled=Phase==1
  END SUBROUTINE
  INTEGER FUNCTION ColumnSlot(I,J) RESULT(C)
    INTEGER, INTENT(IN) :: I,J
    C=0
    IF (.NOT.Enabled) RETURN
    IF (I==46 .AND. J==19) C=1
    IF (I==131 .AND. J==28) C=2
  END FUNCTION
  LOGICAL FUNCTION BRC_WET_ACTIVE(I,J,N) RESULT(Active)
    INTEGER, INTENT(IN) :: I,J,N
    Active=.FALSE.
    IF (ColumnSlot(I,J)==0) RETURN
    IF (N<=0 .OR. N>SIZE(Slots)) ERROR STOP 'BRC wet active species range'
    Active=Slots(N)>0
  END FUNCTION
  SUBROUTINE BRC_WET_COLUMN(I,J,DSpc)
    INTEGER, INTENT(IN) :: I,J
    REAL(fp), INTENT(IN) :: DSpc(:,:,:,:)
    INTEGER :: C
    C=ColumnSlot(I,J)
    IF (C==0) RETURN
    IF (Buffer(C)%Meta(1)/=0 .OR. ANY(DSpc(:,:,I,J)/=0.0_fp)) ERROR STOP 'BRC wet column initialization'
    Buffer(C)%Meta=INT([1,I,J,0],int32)
  END SUBROUTINE
  SUBROUTINE NewRow(C,Kind,Phase,I,J,L,N,NW,R)
    INTEGER, INTENT(IN) :: C,Kind,Phase,I,J,L,N,NW
    INTEGER, INTENT(OUT) :: R
    IF (Buffer(C)%Meta(1)/=1 .OR. L<1 .OR. L>NZ) ERROR STOP 'BRC wet trace bounds'
    R=Buffer(C)%Meta(4)+1
    IF (R>MAXROWS) ERROR STOP 'BRC wet trace capacity'
    Buffer(C)%Meta(4)=R
    Buffer(C)%Rows(:,R)=INT([Kind,Phase,I,J,L,N,NW,R],int32)
    Buffer(C)%Flags(:,R)=0;Buffer(C)%Negative(:,R)=0
    Buffer(C)%Q(:,:,R)=0.0_fp;Buffer(C)%Rain(:,:,R)=0.0_fp;Buffer(C)%V(:,R)=0.0_fp
  END SUBROUTINE
  SUBROUTINE BRC_WET_TRACE(Kind,Phase,I,J,L,N,NW,State_Chm,DSpc,V,Flags)
    INTEGER, INTENT(IN) :: Kind,Phase,I,J,L,N,NW,Flags(8)
    TYPE(ChmState), INTENT(IN) :: State_Chm
    REAL(fp), INTENT(IN) :: DSpc(:,:,:,:),V(NV)
    INTEGER :: C,R,S
    IF (.NOT.BRC_WET_ACTIVE(I,J,N)) RETURN
    C=ColumnSlot(I,J)
    IF (Positions(Slots(N))/=NW .OR. NW<=0) ERROR STOP 'BRC wet trace mapping identity'
    CALL NewRow(C,Kind,Phase,I,J,L,N,NW,R)
    Buffer(C)%Flags(:,R)=INT(Flags,int32);Buffer(C)%V(:,R)=V
    DO S=1,NF
      Buffer(C)%Q(:,S,R)=State_Chm%Species(Ids(S))%Conc(I,J,:)
      IF (Positions(S)>0) Buffer(C)%Rain(:,S,R)=DSpc(Positions(S),:,I,J)
    END DO
    Buffer(C)%Negative(:,R)=INT(MERGE(1,0,Buffer(C)%Q(:,Slots(N),R)<0.0_fp),int32)
  END SUBROUTINE
  SUBROUTINE BRC_WET_SAFETY(I,J,L,N,LS,Spc,DSpc)
    INTEGER, INTENT(IN) :: I,J,L,N
    LOGICAL, INTENT(IN) :: LS
    REAL(fp), INTENT(IN) :: Spc(:),DSpc(:)
    INTEGER :: C,R,S
    IF (.NOT.BRC_WET_ACTIVE(I,J,N)) RETURN
    C=ColumnSlot(I,J);S=Slots(N)
    IF (SIZE(Spc)/=NZ .OR. SIZE(DSpc)/=NZ) ERROR STOP 'BRC wet safety shape'
    CALL NewRow(C,9,0,I,J,L,N,Positions(S),R)
    Buffer(C)%Flags(1:3,R)=INT([MERGE(1,0,LS),1,MERGE(1,0,ANY(Spc<0.0_fp))],int32)
    Buffer(C)%Q(:,S,R)=Spc;Buffer(C)%Rain(:,S,R)=DSpc
    Buffer(C)%Negative(:,R)=INT(MERGE(1,0,Spc<0.0_fp),int32)
  END SUBROUTINE
END MODULE BRC_WETDEP_CAPTURE_MOD
