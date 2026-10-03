! Opt-in READ-ONLY raw parent transport evidence for an offline mechanism test.
! Native fp stream layout is documented in docs/brc/history-and-country-attribution.md.
! No origin partition or parent transport algorithm is changed by this module.
MODULE BRC_TRANSPORT_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Met_Mod, ONLY: MetState
  USE State_Grid_Mod, ONLY: GrdState
  USE UnitConv_Mod, ONLY: KG_SPECIES_PER_KG_DRY_AIR
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_CAPTURE_BEGIN, BRC_CAPTURE_OPEN, BRC_CAPTURE_WRITE, BRC_CAPTURE_CLOSE
  PUBLIC :: BRC_CAPTURE_MIX_PHASE, BRC_CAPTURE_MIX_GUARD, BRC_CAPTURE_MIX_SCALE
  PUBLIC :: BRC_CAPTURE_MIX_DIFF
  CHARACTER(LEN=8), PARAMETER :: Parents(7)=[ CHARACTER(LEN=8) :: &
    'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA' ]
  CHARACTER(LEN=1024), SAVE :: Directory=''
  INTEGER, SAVE :: Ids(7)=0,Step=0,Limit=1
  INTEGER, SAVE :: Tags(7,4)=0,OriginCount=0
  LOGICAL, SAVE :: Initialized=.FALSE.,Capturing=.FALSE.
  LOGICAL, SAVE :: MixingAudit=.FALSE.
  INTEGER, SAVE :: MixIds(35)=0
  INTEGER, ALLOCATABLE, SAVE :: MixUnits(:),MixStages(:)
  INTEGER(int32), ALLOCATABLE, SAVE :: Rollback(:,:,:),SafeScale(:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Threshold(:,:),ScaleNum(:,:,:),ScaleDen(:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: PrePole(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: PreOrigins(:,:,:,:,:)
  INTERFACE BRC_CAPTURE_WRITE
    MODULE PROCEDURE WRITE_2D,WRITE_3D
  END INTERFACE
CONTAINS
  SUBROUTINE BRC_CAPTURE_BEGIN(State_Chm,NX,NY,NZ,JFIRST,JLAST,NG)
    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER, INTENT(IN) :: NX,NY,NZ,JFIRST,JLAST,NG
    INTEGER :: Status,Length,S,O,TagCount
    CHARACTER(LEN=3), PARAMETER :: Origins(4)=['USA','CAN','ROW','UNT']
    CHARACTER(LEN=32) :: Value
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_PARENT_FLUX_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==1 .OR. Length==0) Directory=''
      IF (Status==-1) ERROR STOP 'BRC parent capture directory exceeds1024characters'
      CALL GET_ENVIRONMENT_VARIABLE('BRC_PARENT_FLUX_CAPTURE_MIXING_AUDIT',Value,Length,Status)
      IF (Status==-1) ERROR STOP 'Truncated BRC mixing capture flag'
      IF (Status==0 .AND. Length>0) THEN
        IF (TRIM(Value)/='0' .AND. TRIM(Value)/='1') ERROR STOP 'BRC mixing capture flag must be0or1'
        MixingAudit=TRIM(Value)=='1'
      ENDIF
      IF (MixingAudit .AND. LEN_TRIM(Directory)==0) ERROR STOP 'BRC mixing capture requires capture directory'
      IF (LEN_TRIM(Directory)>0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('BRC_PARENT_FLUX_CAPTURE_STEPS',Value,Length,Status)
        IF (Status==-1) ERROR STOP 'Truncated BRC parent capture step limit'
        IF (Status==0 .AND. Length>0) THEN
          IF (VERIFY(TRIM(ADJUSTL(Value)),'0123456789')/=0) &
            ERROR STOP 'Malformed BRC parent capture step limit'
          READ(Value,*,IOSTAT=Status) Limit
          IF (Status/=0) ERROR STOP 'Invalid BRC parent capture step limit'
        ENDIF
        IF (Limit<1 .OR. Limit>6) ERROR STOP 'BRC parent capture limit must be1..6'
        DO S=1,7
          Ids(S)=Ind_(TRIM(Parents(S)))
          DO O=1,4
            Tags(S,O)=Ind_(TRIM(Parents(S))//'_'//Origins(O))
          ENDDO
        ENDDO
        IF (ANY(Ids<1)) ERROR STOP 'BRC parent capture requires all seven parents'
        ALLOCATE(PrePole(NX,NY,NZ,7))
        TagCount=COUNT(Tags>0)
        IF (TagCount/=0 .AND. TagCount/=28) ERROR STOP 'Incomplete capture origin registration'
        IF (TagCount==28) THEN
          OriginCount=4
          ALLOCATE(PreOrigins(NX,NY,NZ,7,4))
        ENDIF
        IF (MixingAudit) THEN
          IF (TagCount/=28) ERROR STOP 'BRC mixing capture requires28origins'
          DO S=1,7
            MixIds(5*S-4)=Ids(S)
            MixIds(5*S-3:5*S)=Tags(S,:)
          ENDDO
          ALLOCATE(MixUnits(NY),MixStages(NY),Rollback(NX,35,NY),SafeScale(NX,35,NY), &
                   Threshold(35,NY),ScaleNum(NX,35,NY),ScaleDen(NX,35,NY))
          MixUnits=0;MixStages=0
        ENDIF
      ENDIF
    ENDIF
    Step=Step+1
    Capturing=LEN_TRIM(Directory)>0 .AND. Step<=Limit
    IF (.NOT. Capturing) RETURN
    IF (JFIRST/=1 .OR. JLAST/=NY .OR. NG/=0) &
      ERROR STOP 'BRC capture currently requires global grid without halos'
    DO S=1,7
      PrePole(:,:,:,S)=State_Chm%Species(Ids(S))%Conc(:,:,NZ:1:-1)
      DO O=1,OriginCount
        IF (State_Chm%Species(Tags(S,O))%Units/=State_Chm%Species(Ids(S))%Units) &
          ERROR STOP 'Capture parent/origin units differ'
        PreOrigins(:,:,:,S,O)=State_Chm%Species(Tags(S,O))%Conc(:,:,NZ:1:-1)
      ENDDO
    ENDDO
  END SUBROUTINE BRC_CAPTURE_BEGIN

  SUBROUTINE CHECK_IO(Status)
    INTEGER, INTENT(IN) :: Status
    IF (Status/=0) ERROR STOP 'BRC parent capture I/O failed; preserve partial evidence'
  END SUBROUTINE CHECK_IO

  SUBROUTINE BRC_CAPTURE_OPEN(IQ,Units,Dt,Area,Geo,GeoPC,DP1,DP2,Q,CX,CY,WZ, &
                            J1P,J2P,Fill,IORD,JORD,KORD,Cross,Unit)
    INTEGER, INTENT(IN) :: IQ,Units,J1P,J2P,IORD,JORD,KORD
    REAL(fp), INTENT(IN) :: Dt,Area(:),Geo(:),GeoPC,DP1(:,:,:),DP2(:,:,:),Q(:,:,:)
    REAL(fp), INTENT(IN) :: CX(:,:,:),CY(:,:,:),WZ(:,:,:)
    LOGICAL, INTENT(IN) :: Fill,Cross
    INTEGER, INTENT(OUT) :: Unit
    INTEGER :: S,Status,NX,NY,NZ
    INTEGER(int32) :: Header(16)
    CHARACTER(LEN=1200) :: Path
    CHARACTER(LEN=6) :: Counter
    Unit=0
    IF (.NOT. Capturing) RETURN
    S=0
    DO Status=1,7
      IF (Ids(Status)==IQ) S=Status
    ENDDO
    IF (S==0) RETURN
    NX=SIZE(Q,1);NY=SIZE(Q,2);NZ=SIZE(Q,3)
    WRITE(Counter,'(i6.6)') Step
    Path=TRIM(Directory)//'/flux_'//TRIM(Parents(S))//'_step'//Counter//'.bin'
    OPEN(NEWUNIT=Unit,FILE=TRIM(Path),ACCESS='STREAM',FORM='UNFORMATTED', &
         STATUS='NEW',ACTION='WRITE',IOSTAT=Status)
    CALL CHECK_IO(Status)
    Header=[1_int32,INT(NX,int32),INT(NY,int32),INT(NZ,int32),INT(Step,int32), &
      INT(Units,int32),INT(STORAGE_SIZE(Dt)/8,int32),INT(J1P,int32),INT(J2P,int32), &
      INT(MERGE(1,0,Fill),int32),INT(OriginCount,int32),16909060_int32, &
      INT(IORD,int32),INT(JORD,int32),INT(KORD,int32),INT(MERGE(1,0,Cross),int32)]
    WRITE(Unit,IOSTAT=Status) 'BRCFX001',Parents(S),Header,Dt,Area,Geo,GeoPC,PrePole(:,:,:,S)
    CALL CHECK_IO(Status)
    IF (OriginCount>0) THEN
      WRITE(Unit,IOSTAT=Status) PreOrigins(:,:,:,S,:)
      CALL CHECK_IO(Status)
    ENDIF
    WRITE(Unit,IOSTAT=Status) Q,DP1,DP2,CX,CY,WZ
    CALL CHECK_IO(Status)
  END SUBROUTINE BRC_CAPTURE_OPEN

  SUBROUTINE WRITE_2D(Unit,A)
    INTEGER, INTENT(IN) :: Unit
    REAL(fp), INTENT(IN) :: A(:,:)
    INTEGER :: Status
    IF (Unit==0) RETURN
    WRITE(Unit,IOSTAT=Status) A
    CALL CHECK_IO(Status)
  END SUBROUTINE WRITE_2D

  SUBROUTINE WRITE_3D(Unit,A)
    INTEGER, INTENT(IN) :: Unit
    REAL(fp), INTENT(IN) :: A(:,:,:)
    INTEGER :: Status
    IF (Unit==0) RETURN
    WRITE(Unit,IOSTAT=Status) A
    CALL CHECK_IO(Status)
  END SUBROUTINE WRITE_3D

  SUBROUTINE BRC_CAPTURE_CLOSE(Unit,Q,FX,FY,FZ)
    INTEGER, INTENT(IN) :: Unit
    REAL(fp), INTENT(IN) :: Q(:,:,:),FX(:,:,:),FY(:,:,:),FZ(:,:,:)
    INTEGER :: Status
    IF (Unit==0) RETURN
    WRITE(Unit,IOSTAT=Status) Q,FX,FY,FZ
    CALL CHECK_IO(Status)
    CLOSE(Unit,IOSTAT=Status)
    CALL CHECK_IO(Status)
  END SUBROUTINE BRC_CAPTURE_CLOSE
  ! Each latitude owns its stream and scratch, including during outer OpenMP.
  ! All model arguments are read-only. Native kg/kg dry fields may be signed
  ! at intermediate phases. Stream stores full cells, actual guard masks and
  ! native mass-restoration arguments; no inferred partition is written back.
  SUBROUTINE BRC_CAPTURE_MIX_PHASE(Phase,Lat,State_Chm,State_Met,State_Grid,Q,Cflx,Dt,CgFirst,NTop)
    INTEGER, INTENT(IN) :: Phase,Lat,CgFirst,NTop
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(MetState), INTENT(IN) :: State_Met
    TYPE(GrdState), INTENT(IN) :: State_Grid
    REAL(fp), INTENT(IN) :: Q(:,:,:),Cflx(:,:),Dt
    INTEGER :: B,N,NX,NZ,Status
    INTEGER(int32) :: Header(12),PhaseCode
    CHARACTER(LEN=1200) :: Path
    CHARACTER(LEN=6) :: Counter
    CHARACTER(LEN=4) :: Latitude
    IF (.NOT. MixingAudit .OR. .NOT. Capturing) RETURN
    NX=SIZE(Q,1);NZ=SIZE(Q,2)
    IF (Phase/=MixStages(Lat)+1 .OR. Phase<1 .OR. Phase>6) ERROR STOP 'BRC mixing capture phase order'
    IF (Phase==1) THEN
      IF (ANY(MixIds>SIZE(Q,3)) .OR. ANY(MixIds>SIZE(Cflx,2))) &
        ERROR STOP 'BRC mixing capture requires transported origin fields'
      DO B=1,35
        IF (State_Chm%Species(MixIds(B))%Units/=KG_SPECIES_PER_KG_DRY_AIR) &
          ERROR STOP 'BRC mixing capture requires native kg/kg dry units'
      ENDDO
      WRITE(Counter,'(i6.6)') Step;WRITE(Latitude,'(i4.4)') Lat
      Path=TRIM(Directory)//'/mix_step'//Counter//'_lat'//Latitude//'.bin'
      OPEN(NEWUNIT=MixUnits(Lat),FILE=TRIM(Path),ACCESS='STREAM',FORM='UNFORMATTED', &
           STATUS='NEW',ACTION='WRITE',IOSTAT=Status)
      CALL CHECK_IO(Status)
      Header=[1_int32,INT(NX,int32),INT(NZ,int32),INT(Step,int32),INT(Lat,int32), &
              35_int32,INT(STORAGE_SIZE(Dt)/8,int32),INT(KG_SPECIES_PER_KG_DRY_AIR,int32), &
              INT(CgFirst,int32),INT(NTop,int32),16909060_int32,6_int32]
      WRITE(MixUnits(Lat),IOSTAT=Status) 'BRCMX001',Header,Dt, &
        State_Met%AD(:,Lat,NZ:1:-1),State_Grid%AREA_M2(:,Lat)
      CALL CHECK_IO(Status)
      DO B=1,35
        WRITE(MixUnits(Lat),IOSTAT=Status) Cflx(:,MixIds(B))
        CALL CHECK_IO(Status)
      ENDDO
      Rollback(:,:,Lat)=-1;SafeScale(:,:,Lat)=-1
      Threshold(:,Lat)=0;ScaleNum(:,:,Lat)=0;ScaleDen(:,:,Lat)=0
    ENDIF
    PhaseCode=INT(Phase,int32)
    WRITE(MixUnits(Lat),IOSTAT=Status) PhaseCode
    CALL CHECK_IO(Status)
    DO B=1,35
      N=MixIds(B)
      WRITE(MixUnits(Lat),IOSTAT=Status) Q(:,:,N)
      CALL CHECK_IO(Status)
    ENDDO
    MixStages(Lat)=Phase
    IF (Phase==6) THEN
      IF (ANY(Rollback(:,:,Lat)<0) .OR. ANY(SafeScale(:,:,Lat)<0)) ERROR STOP 'Incomplete BRC mixing guard capture'
      WRITE(MixUnits(Lat),IOSTAT=Status) Threshold(:,Lat),Rollback(:,:,Lat), &
        ScaleNum(:,:,Lat),ScaleDen(:,:,Lat),SafeScale(:,:,Lat)
      CALL CHECK_IO(Status)
      CLOSE(MixUnits(Lat),IOSTAT=Status)
      CALL CHECK_IO(Status)
      MixUnits(Lat)=0;MixStages(Lat)=0
    ENDIF
  END SUBROUTINE BRC_CAPTURE_MIX_PHASE

  SUBROUTINE BRC_CAPTURE_MIX_DIFF(Lat,Cc,Ze,Term,Dqbot)
    INTEGER, INTENT(IN) :: Lat
    REAL(fp), INTENT(IN) :: Cc(:,:),Ze(:,:),Term(:,:),Dqbot(:,:)
    INTEGER :: B,Status
    IF (.NOT. MixingAudit .OR. .NOT. Capturing) RETURN
    IF (MixStages(Lat)/=3) ERROR STOP 'BRC mixing coefficients outside rollback phase'
    WRITE(MixUnits(Lat),IOSTAT=Status) 'BRCMD001',Cc,Ze,Term
    CALL CHECK_IO(Status)
    DO B=1,35
      WRITE(MixUnits(Lat),IOSTAT=Status) Dqbot(:,MixIds(B))
      CALL CHECK_IO(Status)
    ENDDO
  END SUBROUTINE BRC_CAPTURE_MIX_DIFF

  SUBROUTINE BRC_CAPTURE_MIX_GUARD(Lat,M,Adjust,Qmin)
    INTEGER, INTENT(IN) :: Lat,M
    LOGICAL, INTENT(IN) :: Adjust(:)
    REAL(fp), INTENT(IN) :: Qmin
    INTEGER :: B
    IF (.NOT. MixingAudit .OR. .NOT. Capturing) RETURN
    DO B=1,35
      IF (MixIds(B)/=M) CYCLE
      IF (MixStages(Lat)/=2) ERROR STOP 'BRC mixing guard outside raw phase'
      Rollback(:,B,Lat)=MERGE(1_int32,0_int32,Adjust)
      Threshold(B,Lat)=Qmin
      RETURN
    ENDDO
  END SUBROUTINE BRC_CAPTURE_MIX_GUARD

  SUBROUTINE BRC_CAPTURE_MIX_SCALE(Lat,M,I,Num,Den,Applied)
    INTEGER, INTENT(IN) :: Lat,M,I
    REAL(fp), INTENT(IN) :: Num,Den
    LOGICAL, INTENT(IN) :: Applied
    INTEGER :: B
    IF (.NOT. MixingAudit .OR. .NOT. Capturing) RETURN
    DO B=1,35
      IF (MixIds(B)/=M) CYCLE
      IF (MixStages(Lat)/=5) ERROR STOP 'BRC mixing scale outside clipped phase'
      ScaleNum(I,B,Lat)=Num;ScaleDen(I,B,Lat)=Den
      SafeScale(I,B,Lat)=MERGE(1_int32,0_int32,Applied)
      RETURN
    ENDDO
  END SUBROUTINE BRC_CAPTURE_MIX_SCALE
END MODULE BRC_TRANSPORT_CAPTURE_MOD
