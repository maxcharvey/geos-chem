! Opt-in read-only native gross surface-budget evidence. No flux is modified.
MODULE BRC_SURFACE_BUDGET_CAPTURE_MOD
  USE Precision_Mod, ONLY: fp
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Met_Mod, ONLY: MetState
  USE State_Grid_Mod, ONLY: GrdState
  USE PhysConstants, ONLY: AIRMW
  USE UnitConv_Mod, ONLY: MOLES_SPECIES_PER_MOLES_DRY_AIR
  USE Time_Mod, ONLY: GET_NYMD, GET_NHMS, GET_TS_CONV, GET_TS_DYN, GET_TS_EMIS
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_CAPTURE_SURFACE_BUDGET
  PUBLIC :: BRC_CAPTURE_MIX_TIME
  CHARACTER(LEN=8), PARAMETER :: Parents(7)=[CHARACTER(LEN=8) :: &
    'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA']
  CHARACTER(LEN=3), PARAMETER :: Origins(4)=['USA','CAN','ROW','UNT']
  CHARACTER(LEN=1024), SAVE :: Directory=''
  CHARACTER(LEN=16), SAVE :: Names(35)
  INTEGER, SAVE :: Ids(35)=0,Advected(35)=0,DryDep(35)=0,Step=0,Limit=2
  LOGICAL, SAVE :: Initialized=.FALSE.
CONTAINS
  SUBROUTINE BRC_CAPTURE_MIX_TIME(Path,MixStep,Lat,Dt,MixIds)
    CHARACTER(LEN=*), INTENT(IN) :: Path
    INTEGER, INTENT(IN) :: MixStep,Lat,MixIds(:)
    REAL(fp), INTENT(IN) :: Dt
    INTEGER :: Unit,Status
    IF (LEN_TRIM(Directory)==0 .OR. MixStep>Limit) RETURN
    IF (SIZE(MixIds)/=35 .OR. ANY(MixIds/=Ids)) &
      ERROR STOP 'BRC gross/mixing capture species maps differ'
    OPEN(NEWUNIT=Unit,FILE=TRIM(Path)//'.time',STATUS='NEW',ACTION='WRITE',IOSTAT=Status)
    CALL CHECK_IO(Status)
    WRITE(Unit,'(a,4(1x,i0),1x,es26.17e3)',IOSTAT=Status) &
      'BRCMT001',MixStep,GET_NYMD(),GET_NHMS(),Lat,Dt
    CALL CHECK_IO(Status)
    WRITE(Unit,'(35(1x,i0))',IOSTAT=Status) MixIds
    CALL CHECK_IO(Status)
    CLOSE(Unit,IOSTAT=Status)
    CALL CHECK_IO(Status)
  END SUBROUTINE BRC_CAPTURE_MIX_TIME

  SUBROUTINE CHECK_IO(Status)
    INTEGER, INTENT(IN) :: Status
    IF (Status/=0) ERROR STOP 'BRC gross surface capture I/O failed; preserve partial file'
  END SUBROUTINE CHECK_IO

  SUBROUTINE BRC_CAPTURE_SURFACE_BUDGET(State_Chm,State_Met,State_Grid,Emissions,Losses)
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(MetState), INTENT(IN) :: State_Met
    TYPE(GrdState), INTENT(IN) :: State_Grid
    REAL(fp), INTENT(IN) :: Emissions(:,:,:),Losses(:,:,:)
    INTEGER :: Status,Length,S,O,N,A,D,MapCount,I,J,Unit,NX,NY
    INTEGER(int32) :: Header(18),Maps(4,35)
    INTEGER(int32), ALLOCATABLE :: TopMix(:,:)
    REAL(fp) :: MW(35),DtConv,DtDyn,DtEmis
    REAL(fp), ALLOCATABLE :: Frequency(:,:)
    CHARACTER(LEN=32) :: Value
    CHARACTER(LEN=1200) :: Path
    CHARACTER(LEN=6) :: Counter,Clock
    CHARACTER(LEN=8) :: Date
    IF (.NOT. Initialized) THEN
      Initialized=.TRUE.
      CALL GET_ENVIRONMENT_VARIABLE('BRC_SURFACE_BUDGET_CAPTURE_DIR',Directory,Length,Status)
      IF (Status==1 .OR. Length==0) Directory=''
      IF (Status==-1) ERROR STOP 'BRC gross surface capture directory truncated'
      IF (LEN_TRIM(Directory)>0) THEN
        CALL GET_ENVIRONMENT_VARIABLE('BRC_SURFACE_BUDGET_CAPTURE_STEPS',Value,Length,Status)
        IF (Status==-1) ERROR STOP 'BRC gross surface capture step flag truncated'
        IF (Status==0 .AND. Length>0) THEN
          IF (VERIFY(TRIM(ADJUSTL(Value)),'0123456789')/=0) &
            ERROR STOP 'Malformed BRC gross surface capture step flag'
          READ(Value,*,IOSTAT=Status) Limit
          IF (Status/=0) ERROR STOP 'Invalid BRC gross surface capture step limit'
        ENDIF
        IF (Limit<1 .OR. Limit>6) ERROR STOP 'BRC gross surface capture limit must be1..6'
        DO S=1,7
          Names(5*S-4)=Parents(S)
          DO O=1,4
            Names(5*S-4+O)=TRIM(Parents(S))//'_'//Origins(O)
          ENDDO
        ENDDO
        DO S=1,35
          Ids(S)=Ind_(TRIM(Names(S)))
          IF (Ids(S)<1 .OR. Ids(S)>State_Chm%nSpecies) &
            ERROR STOP 'BRC gross surface capture missing species'
          IF (COUNT(Ids(:S)==Ids(S))/=1) ERROR STOP 'BRC gross surface capture duplicate species'
          MapCount=0
          DO A=1,State_Chm%nAdvect
            IF (State_Chm%Map_Advect(A)==Ids(S)) THEN
              Advected(S)=A;MapCount=MapCount+1
            ENDIF
          ENDDO
          IF (MapCount/=1) ERROR STOP 'BRC gross surface capture ambiguous advected map'
          MapCount=0
          DO D=1,State_Chm%nDryDep
            IF (State_Chm%Map_DryDep(D)==Ids(S)) THEN
              DryDep(S)=D;MapCount=MapCount+1
            ENDIF
          ENDDO
          IF (MapCount>1) ERROR STOP 'BRC gross surface capture ambiguous drydep map'
          ! No drydep index is valid; map0 explicitly marks the absent frequency.
        ENDDO
      ENDIF
    ENDIF
    IF (LEN_TRIM(Directory)==0) RETURN
    Step=Step+1
    IF (Step>Limit) RETURN
    NX=State_Grid%NX;NY=State_Grid%NY
    IF (State_Grid%NestedGrid .OR. NX/=State_Grid%GlobalNX .OR. NY/=State_Grid%GlobalNY) &
      ERROR STOP 'BRC gross surface capture requires global grid'
    IF (SIZE(Emissions,1)/=NX .OR. SIZE(Emissions,2)/=NY .OR. &
        SIZE(Emissions,3)/=State_Chm%nAdvect .OR. ANY(SHAPE(Losses)/=SHAPE(Emissions))) &
      ERROR STOP 'BRC gross surface capture dimensions mismatch'
    DO S=1,35
      N=Ids(S)
      IF (State_Chm%Species(N)%Units/=MOLES_SPECIES_PER_MOLES_DRY_AIR) &
        ERROR STOP 'BRC gross surface capture requires mol/mol dry'
      Maps(:,S)=[INT(N,int32),INT(Advected(S),int32),INT(DryDep(S),int32), &
                 INT(State_Chm%Species(N)%Units,int32)]
      MW(S)=State_Chm%SpcData(N)%Info%MW_g
    ENDDO
    DtConv=REAL(GET_TS_CONV(),fp);DtDyn=REAL(GET_TS_DYN(),fp);DtEmis=REAL(GET_TS_EMIS(),fp)
    ALLOCATE(TopMix(NX,NY),Frequency(NX,NY))
    DO J=1,NY
    DO I=1,NX
      TopMix(I,J)=INT(MAX(1,FLOOR(State_Met%PBL_TOP_L(I,J))),int32)
    ENDDO
    ENDDO
    WRITE(Counter,'(i6.6)') Step
    WRITE(Date,'(i8.8)') GET_NYMD()
    WRITE(Clock,'(i6.6)') GET_NHMS()
    Path=TRIM(Directory)//'/surface_step'//Counter//'_'//Date//'_'//Clock//'.bin'
    OPEN(NEWUNIT=Unit,FILE=TRIM(Path),ACCESS='STREAM',FORM='UNFORMATTED', &
         STATUS='NEW',ACTION='WRITE',IOSTAT=Status)
    CALL CHECK_IO(Status)
    Header=[1_int32,INT(NX,int32),INT(NY,int32),35_int32,INT(Step,int32), &
      INT(STORAGE_SIZE(DtConv)/8,int32),16909060_int32,INT(GET_NYMD(),int32), &
      INT(GET_NHMS(),int32),INT(MOLES_SPECIES_PER_MOLES_DRY_AIR,int32), &
      1_int32,5_int32,16_int32,7_int32,4_int32,INT(State_Grid%NZ,int32), &
      INT(STORAGE_SIZE(State_Grid%XMid(1,1))/8,int32), &
      INT(STORAGE_SIZE(State_Grid%Area_M2(1,1))/8,int32)]
    WRITE(Unit,IOSTAT=Status) 'BRCED001',Header,Names,Maps,MW,AIRMW,DtConv,DtDyn,DtEmis
    CALL CHECK_IO(Status)
    WRITE(Unit,IOSTAT=Status) State_Grid%XMid,State_Grid%YMid,State_Grid%Area_M2, &
      State_Met%AD(:,:,1),State_Met%PBL_TOP_L,TopMix
    CALL CHECK_IO(Status)
    DO S=1,35
      N=Ids(S);A=Advected(S);D=DryDep(S)
      Frequency=0.0_fp
      IF (D>0) Frequency=State_Chm%DryDepFreq(:,:,D)
      WRITE(Unit,IOSTAT=Status) State_Chm%Species(N)%Conc(:,:,1),Emissions(:,:,A), &
        Losses(:,:,A),State_Chm%SurfaceFlux(:,:,A),Frequency
      CALL CHECK_IO(Status)
    ENDDO
    CLOSE(Unit,IOSTAT=Status)
    CALL CHECK_IO(Status)
    DEALLOCATE(TopMix,Frequency)
  END SUBROUTINE BRC_CAPTURE_SURFACE_BUDGET
END MODULE BRC_SURFACE_BUDGET_CAPTURE_MOD
