!------------------------------------------------------------------------------
! Runtime-gated, process-boundary checkpoints for the plume transport tests.
!
! These diagnostics are disabled unless GC_PLUME_CHECKPOINTS is set to 1,
! true, or yes.  They do not participate in any science calculation.
!------------------------------------------------------------------------------
MODULE Plume_Checkpoint_Mod

  USE ErrCode_Mod,     ONLY : GC_SUCCESS, GC_Error
  USE Input_Opt_Mod,   ONLY : OptInput
  USE Precision_Mod,   ONLY : f8
  USE State_Chm_Mod,   ONLY : ChmState, Ind_
  USE State_Grid_Mod,  ONLY : GrdState
  USE State_Met_Mod,   ONLY : MetState
  USE UnitConv_Mod,    ONLY : KG_SPECIES_PER_KG_DRY_AIR
  USE netcdf

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Write_Plume_Checkpoint

  INTEGER, PARAMETER :: MAX_TAGS = 10
  CHARACTER(LEN=16), PARAMETER :: TAG_NAMES(MAX_TAGS) = (/ &
       'PLUME_SFC       ', 'PLUME_PBL       ', 'PLUME_6535      ', &
       'PLUME_LEV       ', 'PLUME_PROFILE   ', 'PLUME_SFC_PI    ', &
       'PLUME_PBL_PI    ', 'PLUME_6535_PI   ', 'PLUME_LEV_PI    ', &
       'PLUME_PROFILE_PI' /)

  LOGICAL, SAVE :: Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Checkpoints_Enabled = .FALSE.
  INTEGER, SAVE :: N_Active_Tags = 5
  CHARACTER(LEN=512), SAVE :: Checkpoint_Directory = &
       './OutputDir/PlumeCheckpoints'
  CHARACTER(LEN=128), SAVE :: Checkpoint_Run_Id = 'unspecified'

CONTAINS

  SUBROUTINE Initialize_Gate()

    CHARACTER(LEN=512) :: Value
    INTEGER            :: Env_Status

    IF ( Gate_Initialized ) RETURN
    Gate_Initialized = .TRUE.

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_CHECKPOINTS', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status == 0 ) THEN
       SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
          CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
             Checkpoints_Enabled = .TRUE.
          CASE DEFAULT
             Checkpoints_Enabled = .FALSE.
       END SELECT
    ENDIF

    IF ( .not. Checkpoints_Enabled ) RETURN

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_CHECKPOINT_AGED_POOLS', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status == 0 ) THEN
       SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
          CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
             N_Active_Tags = MAX_TAGS
       END SELECT
    ENDIF

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_CHECKPOINT_DIR', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status == 0 .and. LEN_TRIM( Value ) > 0 ) THEN
       Checkpoint_Directory = TRIM( Value )
    ENDIF

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_CHECKPOINT_RUN_ID', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status == 0 .and. LEN_TRIM( Value ) > 0 ) THEN
       Checkpoint_Run_Id = TRIM( Value )
    ENDIF

    WRITE( 6, '(a)' ) 'PLUME CHECKPOINTS ENABLED: ' // &
         TRIM( Checkpoint_Directory )

  END SUBROUTINE Initialize_Gate

  SUBROUTINE Check_Netcdf( Status, Operation, RC )

    INTEGER,          INTENT(IN)    :: Status
    CHARACTER(LEN=*), INTENT(IN)    :: Operation
    INTEGER,          INTENT(INOUT) :: RC

    CHARACTER(LEN=1024) :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc

    IF ( Status == NF90_NOERR ) RETURN
    ThisLoc = ' -> at Write_Plume_Checkpoint (in plume_checkpoint_mod.F90)'
    ErrMsg = TRIM( Operation ) // ': ' // TRIM( NF90_STRERROR( Status ) )
    CALL GC_Error( ErrMsg, RC, ThisLoc )

  END SUBROUTINE Check_Netcdf

  SUBROUTINE Write_Plume_Checkpoint( Checkpoint_Label, Checkpoint_Id, &
                                     State_Chm, State_Grid, State_Met, &
                                     Input_Opt, Model_Date, Model_Time, &
                                     Heartbeat_Index, Elapsed_Seconds, &
                                     State_Time_Offset_Seconds, RC )

    CHARACTER(LEN=*), INTENT(IN)    :: Checkpoint_Label
    INTEGER,          INTENT(IN)    :: Checkpoint_Id
    TYPE(ChmState),   INTENT(IN)    :: State_Chm
    TYPE(GrdState),   INTENT(IN)    :: State_Grid
    TYPE(MetState),   INTENT(IN)    :: State_Met
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt
    INTEGER,          INTENT(IN)    :: Model_Date
    INTEGER,          INTENT(IN)    :: Model_Time
    INTEGER,          INTENT(IN)    :: Heartbeat_Index
    INTEGER,          INTENT(IN)    :: Elapsed_Seconds
    INTEGER,          INTENT(IN)    :: State_Time_Offset_Seconds
    INTEGER,          INTENT(INOUT) :: RC

    CHARACTER(LEN=1024) :: File_Name
    CHARACTER(LEN=32)   :: Timestamp
    CHARACTER(LEN=1024) :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc
    INTEGER             :: Nc_Id, Status
    INTEGER             :: Dim_Lon, Dim_Lat, Dim_Lev, Dim_Tag
    INTEGER             :: Var_Lon, Var_Lat, Var_Lev, Var_Area
    INTEGER             :: Var_Air_Mass, Var_Delp_Dry, Var_Global_Mass
    INTEGER             :: Var_Tag_Id
    INTEGER             :: Var_Mass(MAX_TAGS)
    INTEGER             :: Dims_2d(2), Dims_3d(3)
    INTEGER             :: Model_Id(MAX_TAGS), Tag_Ids(MAX_TAGS)
    INTEGER             :: I, Allocate_Status
    INTEGER, ALLOCATABLE :: Levels(:)
    REAL(f8), ALLOCATABLE :: Longitude(:), Latitude(:)
    REAL(f8), ALLOCATABLE :: Area(:,:), Air_Mass(:,:,:)
    REAL(f8), ALLOCATABLE :: Delp_Dry(:,:,:), Tracer_Mass(:,:,:)
    REAL(f8)              :: Global_Mass(MAX_TAGS)

    RC = GC_SUCCESS
    CALL Initialize_Gate()
    IF ( .not. Checkpoints_Enabled ) RETURN

    ThisLoc = ' -> at Write_Plume_Checkpoint (in plume_checkpoint_mod.F90)'

    IF ( Checkpoint_Id < 0 .or. Checkpoint_Id > 6 ) THEN
       ErrMsg = 'Invalid plume checkpoint identifier'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    DO I = 1, N_Active_Tags
       Model_Id(I) = Ind_( TRIM( TAG_NAMES(I) ), 'S' )
       IF ( Model_Id(I) < 1 ) THEN
          ErrMsg = 'Missing plume species: ' // TRIM( TAG_NAMES(I) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( Ind_( TRIM( TAG_NAMES(I) ), 'A' ) < 1 ) THEN
          ErrMsg = 'Plume species is not advected: ' // TRIM( TAG_NAMES(I) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( State_Chm%Species(Model_Id(I))%Units /= &
            KG_SPECIES_PER_KG_DRY_AIR ) THEN
          ErrMsg = 'Plume checkpoint requires kg/kg dry state units: ' // &
                   TRIM( TAG_NAMES(I) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. ASSOCIATED( State_Chm%Species(Model_Id(I))%Conc ) ) THEN
          ErrMsg = 'Plume species concentration is not associated: ' // &
                   TRIM( TAG_NAMES(I) )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       Tag_Ids(I) = I
    ENDDO

    ALLOCATE( Levels(State_Grid%NZ), Longitude(State_Grid%NX), &
              Latitude(State_Grid%NY), &
              Area(State_Grid%NX,State_Grid%NY), &
              Air_Mass(State_Grid%NX,State_Grid%NY,State_Grid%NZ), &
              Delp_Dry(State_Grid%NX,State_Grid%NY,State_Grid%NZ), &
              Tracer_Mass(State_Grid%NX,State_Grid%NY,State_Grid%NZ), &
              STAT=Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
       ErrMsg = 'Could not allocate plume checkpoint work arrays'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    Levels    = (/ ( I, I=1,State_Grid%NZ ) /)
    Longitude = REAL( State_Grid%XMid(:,1), f8 )
    Latitude  = REAL( State_Grid%YMid(1,:), f8 )
    Area      = REAL( State_Grid%Area_M2, f8 )
    Air_Mass  = REAL( State_Met%AD, f8 )
    Delp_Dry  = REAL( State_Met%DELP_DRY, f8 )

    WRITE( Timestamp, '(i8.8,"_",i6.6,"z")' ) Model_Date, Model_Time
    IF ( Checkpoint_Directory(LEN_TRIM(Checkpoint_Directory): &
                              LEN_TRIM(Checkpoint_Directory)) == '/' ) THEN
       File_Name = TRIM( Checkpoint_Directory ) // &
                   'GEOSChem.PlumeCheckpoint.' // &
                   TRIM( Checkpoint_Label ) // '.' // TRIM( Timestamp ) // '.nc4'
    ELSE
       File_Name = TRIM( Checkpoint_Directory ) // &
                   '/GEOSChem.PlumeCheckpoint.' // &
                   TRIM( Checkpoint_Label ) // '.' // TRIM( Timestamp ) // '.nc4'
    ENDIF

    Status = NF90_CREATE( TRIM( File_Name ), &
                          IOR( NF90_NETCDF4, NF90_NOCLOBBER ), Nc_Id )
    CALL Check_Netcdf( Status, 'Creating ' // TRIM( File_Name ), RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_DIM( Nc_Id, 'lon', State_Grid%NX, Dim_Lon )
    CALL Check_Netcdf( Status, 'Defining lon dimension', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_DEF_DIM( Nc_Id, 'lat', State_Grid%NY, Dim_Lat )
    CALL Check_Netcdf( Status, 'Defining lat dimension', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_DEF_DIM( Nc_Id, 'lev', State_Grid%NZ, Dim_Lev )
    CALL Check_Netcdf( Status, 'Defining lev dimension', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_DEF_DIM( Nc_Id, 'tag', N_Active_Tags, Dim_Tag )
    CALL Check_Netcdf( Status, 'Defining tag dimension', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_VAR( Nc_Id, 'lon', NF90_DOUBLE, (/ Dim_Lon /), Var_Lon )
    CALL Check_Netcdf( Status, 'Defining lon variable', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Lon, 'units', 'degrees_east' )
    CALL Check_Netcdf( Status, 'Defining lon units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_VAR( Nc_Id, 'lat', NF90_DOUBLE, (/ Dim_Lat /), Var_Lat )
    CALL Check_Netcdf( Status, 'Defining lat variable', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Lat, 'units', 'degrees_north' )
    CALL Check_Netcdf( Status, 'Defining lat units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_VAR( Nc_Id, 'lev', NF90_INT, (/ Dim_Lev /), Var_Lev )
    CALL Check_Netcdf( Status, 'Defining lev variable', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Lev, 'long_name', &
                           'one-based GCClassic model level' )
    CALL Check_Netcdf( Status, 'Defining lev metadata', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Dims_2d = (/ Dim_Lon, Dim_Lat /)
    Status = NF90_DEF_VAR( Nc_Id, 'grid_cell_area', NF90_DOUBLE, Dims_2d, &
                           Var_Area, SHUFFLE=.TRUE., DEFLATE_LEVEL=1 )
    CALL Check_Netcdf( Status, 'Defining grid-cell area', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Area, 'units', 'm2' )
    CALL Check_Netcdf( Status, 'Defining grid-cell-area units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Dims_3d = (/ Dim_Lon, Dim_Lat, Dim_Lev /)
    Status = NF90_DEF_VAR( Nc_Id, 'dry_air_mass', NF90_DOUBLE, Dims_3d, &
                           Var_Air_Mass, SHUFFLE=.TRUE., DEFLATE_LEVEL=1 )
    CALL Check_Netcdf( Status, 'Defining dry-air mass', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Air_Mass, 'units', 'kg' )
    CALL Check_Netcdf( Status, 'Defining dry-air-mass units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_VAR( Nc_Id, 'Met_DELPDRY', NF90_DOUBLE, Dims_3d, &
                           Var_Delp_Dry, SHUFFLE=.TRUE., DEFLATE_LEVEL=1 )
    CALL Check_Netcdf( Status, 'Defining dry pressure thickness', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Delp_Dry, 'units', 'hPa' )
    CALL Check_Netcdf( Status, 'Defining DELPDRY units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    DO I = 1, N_Active_Tags
       Status = NF90_DEF_VAR( Nc_Id, 'PlumeMass_' // TRIM( TAG_NAMES(I) ), &
                              NF90_DOUBLE, Dims_3d, Var_Mass(I), &
                              SHUFFLE=.TRUE., DEFLATE_LEVEL=1 )
       CALL Check_Netcdf( Status, 'Defining plume mass variable', RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       Status = NF90_PUT_ATT( Nc_Id, Var_Mass(I), 'units', 'kg' )
       CALL Check_Netcdf( Status, 'Defining plume mass units', RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDDO

    Status = NF90_DEF_VAR( Nc_Id, 'tag_id', NF90_INT, (/ Dim_Tag /), &
                           Var_Tag_Id )
    CALL Check_Netcdf( Status, 'Defining tag identifier', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    ! Reuse this id only until tag_id is written, then define global mass.
    IF ( N_Active_Tags == MAX_TAGS ) THEN
       Status = NF90_PUT_ATT( Nc_Id, Var_Tag_Id, 'tag_names', &
            'PLUME_SFC,PLUME_PBL,PLUME_6535,PLUME_LEV,PLUME_PROFILE,' // &
            'PLUME_SFC_PI,PLUME_PBL_PI,PLUME_6535_PI,PLUME_LEV_PI,' // &
            'PLUME_PROFILE_PI' )
    ELSE
       Status = NF90_PUT_ATT( Nc_Id, Var_Tag_Id, 'tag_names', &
            'PLUME_SFC,PLUME_PBL,PLUME_6535,PLUME_LEV,PLUME_PROFILE' )
    ENDIF
    CALL Check_Netcdf( Status, 'Defining tag names', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_DEF_VAR( Nc_Id, 'global_plume_mass', NF90_DOUBLE, &
                           (/ Dim_Tag /), Var_Global_Mass )
    CALL Check_Netcdf( Status, 'Defining global plume mass', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, Var_Global_Mass, 'units', 'kg' )
    CALL Check_Netcdf( Status, 'Defining global plume mass units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'schema_version', &
                           'plume-checkpoint-v1' )
    CALL Check_Netcdf( Status, 'Writing schema version', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'run_id', &
                           TRIM( Checkpoint_Run_Id ) )
    CALL Check_Netcdf( Status, 'Writing run id', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'checkpoint_label', &
                           TRIM( Checkpoint_Label ) )
    CALL Check_Netcdf( Status, 'Writing checkpoint label', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'checkpoint_id', Checkpoint_Id )
    CALL Check_Netcdf( Status, 'Writing checkpoint id', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'model_date', Model_Date )
    CALL Check_Netcdf( Status, 'Writing model date', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'model_time', Model_Time )
    CALL Check_Netcdf( Status, 'Writing model time', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'heartbeat_index', &
                           Heartbeat_Index )
    CALL Check_Netcdf( Status, 'Writing heartbeat index', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'elapsed_seconds', &
                           Elapsed_Seconds )
    CALL Check_Netcdf( Status, 'Writing elapsed seconds', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'state_time_offset_seconds', &
                           State_Time_Offset_Seconds )
    CALL Check_Netcdf( Status, 'Writing state time offset', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'transport_active', &
                           MERGE( 1, 0, Input_Opt%LTRAN ) )
    CALL Check_Netcdf( Status, 'Writing transport switch', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'transport_timestep_seconds', &
                           Input_Opt%TS_DYN )
    CALL Check_Netcdf( Status, 'Writing transport timestep', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'emissions_timestep_seconds', &
                           Input_Opt%TS_EMIS )
    CALL Check_Netcdf( Status, 'Writing emissions timestep', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'chemistry_timestep_seconds', &
                           Input_Opt%TS_CHEM )
    CALL Check_Netcdf( Status, 'Writing chemistry timestep', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'tracer_state_units', &
                           'kg species kg-1 dry air' )
    CALL Check_Netcdf( Status, 'Writing tracer-state units', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_ATT( Nc_Id, NF90_GLOBAL, 'mass_formula', &
                           'PlumeMass = SpeciesConc * dry_air_mass' )
    CALL Check_Netcdf( Status, 'Writing mass formula', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_ENDDEF( Nc_Id )
    CALL Check_Netcdf( Status, 'Ending checkpoint definition', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_PUT_VAR( Nc_Id, Var_Lon, Longitude )
    CALL Check_Netcdf( Status, 'Writing longitude', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Lat, Latitude )
    CALL Check_Netcdf( Status, 'Writing latitude', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Lev, Levels )
    CALL Check_Netcdf( Status, 'Writing levels', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Area, Area )
    CALL Check_Netcdf( Status, 'Writing grid-cell area', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Air_Mass, Air_Mass )
    CALL Check_Netcdf( Status, 'Writing dry-air mass', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Delp_Dry, Delp_Dry )
    CALL Check_Netcdf( Status, 'Writing dry pressure thickness', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    DO I = 1, N_Active_Tags
       Tracer_Mass = REAL( State_Chm%Species(Model_Id(I))%Conc, f8 ) * Air_Mass
       Global_Mass(I) = SUM( Tracer_Mass )
       Status = NF90_PUT_VAR( Nc_Id, Var_Mass(I), Tracer_Mass )
       CALL Check_Netcdf( Status, 'Writing ' // TRIM( TAG_NAMES(I) ), RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDDO

    Status = NF90_PUT_VAR( Nc_Id, Var_Tag_Id, Tag_Ids(1:N_Active_Tags) )
    CALL Check_Netcdf( Status, 'Writing tag identifiers', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_PUT_VAR( Nc_Id, Var_Global_Mass, &
                           Global_Mass(1:N_Active_Tags) )
    CALL Check_Netcdf( Status, 'Writing global plume mass', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    Status = NF90_SYNC( Nc_Id )
    CALL Check_Netcdf( Status, 'Synchronizing checkpoint file', RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    Status = NF90_CLOSE( Nc_Id )
    CALL Check_Netcdf( Status, 'Closing checkpoint file', RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    DEALLOCATE( Levels, Longitude, Latitude, Area, Air_Mass, Delp_Dry, &
                Tracer_Mass )

  END SUBROUTINE Write_Plume_Checkpoint

END MODULE Plume_Checkpoint_Mod
