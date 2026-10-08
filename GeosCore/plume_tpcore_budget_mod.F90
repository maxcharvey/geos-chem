!------------------------------------------------------------------------------
! Runtime-gated TPCORE mass-boundary ledger for the controlled plume tests.
!
! This module is deliberately diagnostic-only.  It never modifies a tracer
! tendency or state.  The ledger is enabled only with
! GC_PLUME_TPCORE_BUDGETS=1 and uses a deterministic serial traversal for its
! global reductions, so its values do not depend on an OpenMP reduction order.
!------------------------------------------------------------------------------
MODULE Plume_Tpcore_Budget_Mod

  USE ERROR_MOD,         ONLY : Error_Stop
  USE OMP_LIB,           ONLY : OMP_GET_MAX_THREADS
  USE PhysConstants,     ONLY : g0_100
  USE Precision_Mod,     ONLY : fp, f8
  USE Qck_Positivity_Mod, ONLY : QCK_BOTTOM_EXACT, QCK_BOTTOM_MICROCLOSURE, &
                                 QCK_BOTTOM_ROUNDOFF, Canonical_Diagnostic_Time
  USE State_Chm_Mod,     ONLY : ChmState, Ind_
  USE State_Met_Mod,     ONLY : MetState
  USE TIME_MOD,          ONLY : GET_ELAPSED_SEC, GET_NHMSb, GET_NYMDb
  USE UnitConv_Mod,      ONLY : KG_SPECIES_PER_KG_DRY_AIR

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: N_PLUME_TPCORE_TAGS = 5
  CHARACTER(LEN=16), PARAMETER :: TAG_NAMES(N_PLUME_TPCORE_TAGS) = (/ &
       'PLUME_SFC       ', 'PLUME_PBL       ', 'PLUME_6535      ', &
       'PLUME_LEV       ', 'PLUME_PROFILE   ' /)

  TYPE, PUBLIC :: Plume_Tpcore_Metrics
     REAL(f8) :: Global_Mass_Kg          = 0.0_f8
     INTEGER  :: Negative_Cell_Count      = 0
     REAL(f8) :: Minimum_Mixing_Ratio     = 0.0_f8
     REAL(f8) :: Total_Negative_Mass_Kg   = 0.0_f8
  END TYPE Plume_Tpcore_Metrics

  PUBLIC :: Compute_Conc_Metrics
  PUBLIC :: Compute_Dq_Metrics
  PUBLIC :: Tpcore_Budget_Begin
  PUBLIC :: Tpcore_Budget_End
  PUBLIC :: Tpcore_Budget_Enabled
  PUBLIC :: Tpcore_Budget_Tag_Index
  PUBLIC :: Tpcore_Budget_Write_Positivity_Boundaries
  PUBLIC :: Tpcore_Cell_Diagnostics_Enabled
  PUBLIC :: Tpcore_Cell_Diag_Begin
  PUBLIC :: Tpcore_Cell_Diag_Capture_Dq_After_Horizontal
  PUBLIC :: Tpcore_Cell_Diag_Capture_Dq_After_Fzppm
  PUBLIC :: Tpcore_Cell_Diag_Capture_Pre_Floor
  PUBLIC :: Tpcore_Cell_Diag_Capture_Post_Floor
  PUBLIC :: Tpcore_Cell_Diag_Record_Qckxyz
  PUBLIC :: Tpcore_Cell_Diag_Write

  INTEGER, PARAMETER, PUBLIC :: PLUME_QCK_TOP      = 1
  INTEGER, PARAMETER, PUBLIC :: PLUME_QCK_INTERIOR = 2
  INTEGER, PARAMETER, PUBLIC :: PLUME_QCK_BOTTOM   = 3

  LOGICAL, SAVE :: Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Budgets_Enabled  = .FALSE.
  LOGICAL, SAVE :: Call_Active      = .FALSE.
  LOGICAL, SAVE :: File_Open        = .FALSE.
  LOGICAL, SAVE :: Cell_Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Cell_Diagnostics_Enabled_Flag = .FALSE.
  LOGICAL, SAVE :: Cell_File_Open = .FALSE.
  LOGICAL, SAVE :: Donor_Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Donor_Diagnostics_Enabled_Flag = .FALSE.
  LOGICAL, SAVE :: Donor_File_Open = .FALSE.
  INTEGER, SAVE :: Budget_Unit      = -1
  INTEGER, SAVE :: Cell_Unit        = -1
  INTEGER, SAVE :: Donor_Unit       = -1
  INTEGER, SAVE :: Call_Index       = 0
  INTEGER, SAVE :: Tag_Ids(N_PLUME_TPCORE_TAGS) = 0
  INTEGER, SAVE :: Model_Date       = 0
  INTEGER, SAVE :: Model_Time       = 0
  INTEGER, SAVE :: Elapsed_Seconds  = 0
  INTEGER, SAVE :: Heartbeat_Index  = 0
  INTEGER, SAVE :: Omp_Thread_Count = 0
  INTEGER, SAVE :: Qckxyz_Invoked   = 0
  CHARACTER(LEN=128), SAVE :: Budget_Run_Id = ''
  CHARACTER(LEN=128), SAVE :: Budget_Manifest_Id = ''
  REAL(f8), SAVE :: Previous_Mass(N_PLUME_TPCORE_TAGS) = 0.0_f8

  INTEGER, SAVE :: Cell_IM = 0
  INTEGER, SAVE :: Cell_JM = 0
  INTEGER, SAVE :: Cell_KM = 0
  INTEGER, ALLOCATABLE, SAVE :: Cell_Qck_Event_Type(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Dq_After_Horizontal(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Dq_After_Fzppm(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_At_Correction(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Available_Above(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Deficit(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Withdrawn(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Column_Total_Above(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Column_Positive_Above(:,:,:,:)
  INTEGER, ALLOCATABLE, SAVE :: Cell_Qck_Correction_Status(:,:,:,:)
  INTEGER, ALLOCATABLE, SAVE :: Cell_Qck_Donor_Count(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Full_Column_Withdrawn(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Declared_Closure(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Correction_Tolerance(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Microclosure_Event_Max_Kg(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qck_Microclosure_Call_Max_Kg(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qptr_Pre_Floor(:,:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: Cell_Qptr_Post_Floor(:,:,:,:)

CONTAINS

  SUBROUTINE Initialize_Gate()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( Gate_Initialized ) RETURN
    Gate_Initialized = .TRUE.

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_BUDGETS', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 ) RETURN
    SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
       CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
          Budgets_Enabled = .TRUE.
       CASE DEFAULT
          RETURN
    END SELECT

    ThisLoc = ' -> at Initialize_Gate (in plume_tpcore_budget_mod.F90)'

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_BUDGET_RUN_ID', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_BUDGET_RUN_ID is required when budgets are enabled'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    Budget_Run_Id = TRIM( Value )

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_BUDGET_MANIFEST_ID', &
                                   Value, STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_BUDGET_MANIFEST_ID is required when budgets are enabled'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    Budget_Manifest_Id = TRIM( Value )

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_BUDGET_FILE', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_BUDGET_FILE is required when budgets are enabled'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    OPEN( NEWUNIT=Budget_Unit, FILE=TRIM( Value ), STATUS='NEW', &
          ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
    IF ( IO_Status /= 0 ) THEN
       ErrMsg = 'Could not create TPCORE budget ledger: ' // TRIM( Value )
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    File_Open = .TRUE.

    WRITE( Budget_Unit, '(a)' ) &
         'schema_version,run_id,manifest_id,model_date,model_time,' // &
         'elapsed_seconds,heartbeat_index,omp_thread_count,tpcore_call_index,' // &
         'boundary,tag,global_mass_kg,boundary_delta_kg,global_dry_air_mass_kg,' // &
         'negative_cell_count,minimum_mixing_ratio,total_negative_mass_kg,' // &
         'corrected_cell_count,correction_mass_delta_kg,qckxyz_invoked'
    FLUSH( Budget_Unit )
    WRITE( 6, '(a)' ) 'PLUME TPCORE BUDGETS ENABLED: ' // TRIM( Value )

    CALL Initialize_Cell_Diagnostics()
    CALL Initialize_Donor_Diagnostics()

  END SUBROUTINE Initialize_Gate

  SUBROUTINE Initialize_Cell_Diagnostics()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( Cell_Gate_Initialized ) RETURN
    Cell_Gate_Initialized = .TRUE.
    IF ( .not. Budgets_Enabled ) RETURN

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_CELL_DIAGNOSTICS', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 ) RETURN
    SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
       CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
          CONTINUE
       CASE DEFAULT
          RETURN
    END SELECT

    ThisLoc = ' -> at Initialize_Cell_Diagnostics ' // &
              '(in plume_tpcore_budget_mod.F90)'
    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_CELL_DIAG_FILE', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_CELL_DIAG_FILE is required when cell diagnostics are enabled'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    OPEN( NEWUNIT=Cell_Unit, FILE=TRIM( Value ), STATUS='NEW', &
          ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
    IF ( IO_Status /= 0 ) THEN
       ErrMsg = 'Could not create TPCORE cell diagnostic ledger: ' // TRIM( Value )
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    Cell_File_Open = .TRUE.
    Cell_Diagnostics_Enabled_Flag = .TRUE.

    WRITE( Cell_Unit, '(a)' ) &
         'schema_version,run_id,manifest_id,model_date,model_time,' // &
         'elapsed_seconds,heartbeat_index,omp_thread_count,tpcore_call_index,' // &
         'tag,event_type,i_tpcore,j_tpcore,k_tpcore,k_model,' // &
         'dq_after_horizontal_hpa,dq_after_fzppm_hpa,dq_at_correction_hpa,' // &
         'deficit_hpa,available_above_hpa,withdrawn_from_above_hpa,' // &
         'unfilled_deficit_hpa,q_before_kgkg,q_after_kgkg,delp_hpa,area_m2,' // &
         'event_mass_delta_kg,correction_policy,correction_outcome,' // &
         'donor_count,full_column_withdrawn_hpa,' // &
         'declared_closure_hpa,ordinary_roundoff_tolerance_hpa,' // &
         'microclosure_event_max_kg,microclosure_call_max_kg'
    FLUSH( Cell_Unit )
    WRITE( 6, '(a)' ) 'PLUME TPCORE CELL DIAGNOSTICS ENABLED: ' // TRIM( Value )

  END SUBROUTINE Initialize_Cell_Diagnostics

  SUBROUTINE Initialize_Donor_Diagnostics()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( Donor_Gate_Initialized ) RETURN
    Donor_Gate_Initialized = .TRUE.
    IF ( .not. Budgets_Enabled ) RETURN

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_DONOR_DIAGNOSTICS', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 ) RETURN
    SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
       CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
          CONTINUE
       CASE DEFAULT
          RETURN
    END SELECT

    ThisLoc = ' -> at Initialize_Donor_Diagnostics ' // &
              '(in plume_tpcore_budget_mod.F90)'
    IF ( .not. Cell_Diagnostics_Enabled_Flag ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_DONOR_DIAGNOSTICS requires cell diagnostics'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_TPCORE_DONOR_DIAG_FILE', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       ErrMsg = 'GC_PLUME_TPCORE_DONOR_DIAG_FILE is required when donor diagnostics are enabled'
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF

    OPEN( NEWUNIT=Donor_Unit, FILE=TRIM( Value ), STATUS='NEW', &
          ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
    IF ( IO_Status /= 0 ) THEN
       ErrMsg = 'Could not create TPCORE donor diagnostic ledger: ' // TRIM( Value )
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    Donor_File_Open = .TRUE.
    Donor_Diagnostics_Enabled_Flag = .TRUE.

    WRITE( Donor_Unit, '(a)' ) &
         'schema_version,run_id,manifest_id,model_date,model_time,' // &
         'elapsed_seconds,heartbeat_index,omp_thread_count,tpcore_call_index,' // &
         'tag,i_tpcore,j_tpcore,k_tpcore,k_model,deficit_hpa,' // &
         'immediate_available_above_hpa,withdrawn_from_immediate_hpa,' // &
         'column_total_above_hpa,column_positive_above_hpa,column_net_mass_hpa,' // &
         'unfilled_by_immediate_hpa,unfillable_column_mass_hpa,area_m2,' // &
         'event_mass_delta_kg,correction_policy,correction_outcome,' // &
         'donor_count,full_column_withdrawn_hpa,' // &
         'declared_closure_hpa,ordinary_roundoff_tolerance_hpa,' // &
         'microclosure_event_max_kg,microclosure_call_max_kg'
    FLUSH( Donor_Unit )
    WRITE( 6, '(a)' ) 'PLUME TPCORE DONOR DIAGNOSTICS ENABLED: ' // TRIM( Value )

  END SUBROUTINE Initialize_Donor_Diagnostics

  LOGICAL FUNCTION Tpcore_Budget_Enabled()

    CALL Initialize_Gate()
    Tpcore_Budget_Enabled = Budgets_Enabled

  END FUNCTION Tpcore_Budget_Enabled

  LOGICAL FUNCTION Tpcore_Cell_Diagnostics_Enabled()

    CALL Initialize_Gate()
    Tpcore_Cell_Diagnostics_Enabled = Cell_Diagnostics_Enabled_Flag

  END FUNCTION Tpcore_Cell_Diagnostics_Enabled

  SUBROUTINE Tpcore_Cell_Diag_Begin( IM, JM, KM )

    INTEGER, INTENT(IN) :: IM, JM, KM

    REAL(f8) :: Required_MiB

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    IF ( .not. Call_Active ) THEN
       CALL Budget_Error( 'Cell diagnostics began without an active TPCORE call', &
                          'Tpcore_Cell_Diag_Begin' )
    ENDIF
    IF ( IM < 1 .or. JM < 1 .or. KM < 1 ) THEN
       CALL Budget_Error( 'Invalid TPCORE dimensions for cell diagnostics', &
                          'Tpcore_Cell_Diag_Begin' )
    ENDIF

    ! Fifteen REAL(fp) fields plus three integer fields are retained
    ! for five diagnostic tracers.  Keep this opt-in probe bounded so that it
    ! cannot silently allocate an excessive high-resolution scratch volume.
    Required_MiB = REAL( IM, f8 ) * REAL( JM, f8 ) * REAL( KM, f8 ) * &
         REAL( N_PLUME_TPCORE_TAGS, f8 ) * &
         ( 15.0_f8 * REAL( STORAGE_SIZE( 0.0_fp ) / 8, f8 ) + &
           3.0_f8 * REAL( STORAGE_SIZE( 0 ) / 8, f8 ) ) / ( 1024.0_f8**2 )
    IF ( Required_MiB > 512.0_f8 ) THEN
       CALL Budget_Error( 'Cell diagnostic scratch exceeds the 512 MiB guard', &
                          'Tpcore_Cell_Diag_Begin' )
    ENDIF

    IF ( ALLOCATED( Cell_Qck_Event_Type ) ) THEN
       IF ( SIZE( Cell_Qck_Event_Type, 1 ) /= IM .or. &
            SIZE( Cell_Qck_Event_Type, 2 ) /= JM .or. &
            SIZE( Cell_Qck_Event_Type, 3 ) /= KM ) THEN
          DEALLOCATE( Cell_Qck_Event_Type, Cell_Dq_After_Horizontal, &
                      Cell_Dq_After_Fzppm, Cell_Qck_At_Correction, &
                      Cell_Qck_Available_Above, Cell_Qck_Deficit, &
                      Cell_Qck_Withdrawn, Cell_Qck_Column_Total_Above, &
                      Cell_Qck_Column_Positive_Above, &
                      Cell_Qck_Correction_Status, Cell_Qck_Donor_Count, &
                      Cell_Qck_Full_Column_Withdrawn, Cell_Qck_Declared_Closure, &
                      Cell_Qck_Correction_Tolerance, &
                      Cell_Qck_Microclosure_Event_Max_Kg, &
                      Cell_Qck_Microclosure_Call_Max_Kg, Cell_Qptr_Pre_Floor, &
                      Cell_Qptr_Post_Floor )
       ENDIF
    ENDIF
    IF ( .not. ALLOCATED( Cell_Qck_Event_Type ) ) THEN
       ALLOCATE( Cell_Qck_Event_Type( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Dq_After_Horizontal( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Dq_After_Fzppm( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_At_Correction( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Available_Above( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Deficit( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Withdrawn( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Column_Total_Above( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Column_Positive_Above( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Correction_Status( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Donor_Count( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Full_Column_Withdrawn( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Declared_Closure( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Correction_Tolerance( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Microclosure_Event_Max_Kg( IM, JM, KM, &
                                                      N_PLUME_TPCORE_TAGS ), &
                 Cell_Qck_Microclosure_Call_Max_Kg( IM, JM, KM, &
                                                     N_PLUME_TPCORE_TAGS ), &
                 Cell_Qptr_Pre_Floor( IM, JM, KM, N_PLUME_TPCORE_TAGS ), &
                 Cell_Qptr_Post_Floor( IM, JM, KM, N_PLUME_TPCORE_TAGS ) )
    ENDIF

    Cell_IM = IM
    Cell_JM = JM
    Cell_KM = KM
    Cell_Qck_Event_Type      = 0
    Cell_Dq_After_Horizontal = 0.0_fp
    Cell_Dq_After_Fzppm      = 0.0_fp
    Cell_Qck_At_Correction   = 0.0_fp
    Cell_Qck_Available_Above = 0.0_fp
    Cell_Qck_Deficit         = 0.0_fp
    Cell_Qck_Withdrawn       = 0.0_fp
    Cell_Qck_Column_Total_Above    = 0.0_fp
    Cell_Qck_Column_Positive_Above = 0.0_fp
    Cell_Qck_Correction_Status      = QCK_BOTTOM_EXACT
    Cell_Qck_Donor_Count            = 0
    Cell_Qck_Full_Column_Withdrawn  = 0.0_fp
    Cell_Qck_Declared_Closure       = 0.0_fp
    Cell_Qck_Correction_Tolerance   = 0.0_fp
    Cell_Qck_Microclosure_Event_Max_Kg = 0.0_fp
    Cell_Qck_Microclosure_Call_Max_Kg  = 0.0_fp
    Cell_Qptr_Pre_Floor      = 0.0_fp
    Cell_Qptr_Post_Floor     = 0.0_fp

  END SUBROUTINE Tpcore_Cell_Diag_Begin

  SUBROUTINE Tpcore_Cell_Diag_Capture_Dq_After_Horizontal( Tag_Index, Dq )

    INTEGER,  INTENT(IN) :: Tag_Index
    REAL(fp), INTENT(IN) :: Dq(:,:,:)

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    Cell_Dq_After_Horizontal(:,:,:,Tag_Index) = Dq

  END SUBROUTINE Tpcore_Cell_Diag_Capture_Dq_After_Horizontal

  SUBROUTINE Tpcore_Cell_Diag_Capture_Dq_After_Fzppm( Tag_Index, Dq )

    INTEGER,  INTENT(IN) :: Tag_Index
    REAL(fp), INTENT(IN) :: Dq(:,:,:)

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    Cell_Dq_After_Fzppm(:,:,:,Tag_Index) = Dq

  END SUBROUTINE Tpcore_Cell_Diag_Capture_Dq_After_Fzppm

  SUBROUTINE Tpcore_Cell_Diag_Capture_Pre_Floor( Tag_Index, Conc )

    INTEGER,  INTENT(IN) :: Tag_Index
    REAL(fp), INTENT(IN) :: Conc(:,:,:)

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    Cell_Qptr_Pre_Floor(:,:,:,Tag_Index) = Conc

  END SUBROUTINE Tpcore_Cell_Diag_Capture_Pre_Floor

  SUBROUTINE Tpcore_Cell_Diag_Capture_Post_Floor( Tag_Index, Conc )

    INTEGER,  INTENT(IN) :: Tag_Index
    REAL(fp), INTENT(IN) :: Conc(:,:,:)

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    Cell_Qptr_Post_Floor(:,:,:,Tag_Index) = Conc

  END SUBROUTINE Tpcore_Cell_Diag_Capture_Post_Floor

  SUBROUTINE Tpcore_Cell_Diag_Record_Qckxyz( Tag_Index, Event_Type, I, J, K, &
                                              At_Correction, Available_Above, &
                                              Deficit, Withdrawn, &
                                              Column_Total_Above, &
                                              Column_Positive_Above, &
                                              Correction_Status, Donor_Count, &
                                              Full_Column_Withdrawn, &
                                              Declared_Closure, &
                                              Correction_Tolerance, &
                                              Microclosure_Event_Max_Kg, &
                                              Microclosure_Call_Max_Kg )

    INTEGER,  INTENT(IN) :: Tag_Index, Event_Type, I, J, K
    REAL(fp), INTENT(IN) :: At_Correction, Available_Above, Deficit, Withdrawn
    REAL(fp), INTENT(IN) :: Column_Total_Above, Column_Positive_Above
    INTEGER,  INTENT(IN), OPTIONAL :: Correction_Status, Donor_Count
    REAL(fp), INTENT(IN), OPTIONAL :: Full_Column_Withdrawn, Declared_Closure
    REAL(fp), INTENT(IN), OPTIONAL :: Correction_Tolerance
    REAL(fp), INTENT(IN), OPTIONAL :: Microclosure_Event_Max_Kg
    REAL(fp), INTENT(IN), OPTIONAL :: Microclosure_Call_Max_Kg

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    Cell_Qck_Event_Type(I,J,K,Tag_Index)      = Event_Type
    Cell_Qck_At_Correction(I,J,K,Tag_Index)   = At_Correction
    Cell_Qck_Available_Above(I,J,K,Tag_Index) = Available_Above
    Cell_Qck_Deficit(I,J,K,Tag_Index)         = Deficit
    Cell_Qck_Withdrawn(I,J,K,Tag_Index)       = Withdrawn
    Cell_Qck_Column_Total_Above(I,J,K,Tag_Index) = Column_Total_Above
    Cell_Qck_Column_Positive_Above(I,J,K,Tag_Index) = Column_Positive_Above
    Cell_Qck_Correction_Status(I,J,K,Tag_Index) = QCK_BOTTOM_EXACT
    Cell_Qck_Donor_Count(I,J,K,Tag_Index) = 0
    Cell_Qck_Full_Column_Withdrawn(I,J,K,Tag_Index) = 0.0_fp
    Cell_Qck_Declared_Closure(I,J,K,Tag_Index) = 0.0_fp
    Cell_Qck_Correction_Tolerance(I,J,K,Tag_Index) = 0.0_fp
    Cell_Qck_Microclosure_Event_Max_Kg(I,J,K,Tag_Index) = 0.0_fp
    Cell_Qck_Microclosure_Call_Max_Kg(I,J,K,Tag_Index) = 0.0_fp
    IF ( PRESENT( Correction_Status ) ) THEN
       Cell_Qck_Correction_Status(I,J,K,Tag_Index) = Correction_Status
    ENDIF
    IF ( PRESENT( Donor_Count ) ) THEN
       Cell_Qck_Donor_Count(I,J,K,Tag_Index) = Donor_Count
    ENDIF
    IF ( PRESENT( Full_Column_Withdrawn ) ) THEN
       Cell_Qck_Full_Column_Withdrawn(I,J,K,Tag_Index) = Full_Column_Withdrawn
    ENDIF
    IF ( PRESENT( Declared_Closure ) ) THEN
       Cell_Qck_Declared_Closure(I,J,K,Tag_Index) = Declared_Closure
    ENDIF
    IF ( PRESENT( Correction_Tolerance ) ) THEN
       Cell_Qck_Correction_Tolerance(I,J,K,Tag_Index) = Correction_Tolerance
    ENDIF
    IF ( PRESENT( Microclosure_Event_Max_Kg ) ) THEN
       Cell_Qck_Microclosure_Event_Max_Kg(I,J,K,Tag_Index) = &
            Microclosure_Event_Max_Kg
    ENDIF
    IF ( PRESENT( Microclosure_Call_Max_Kg ) ) THEN
       Cell_Qck_Microclosure_Call_Max_Kg(I,J,K,Tag_Index) = &
            Microclosure_Call_Max_Kg
    ENDIF

  END SUBROUTINE Tpcore_Cell_Diag_Record_Qckxyz

  SUBROUTINE Tpcore_Cell_Diag_Write( Delp2, Area_M2 )

    REAL(fp), INTENT(IN) :: Delp2(:,:,:)
    REAL(fp), INTENT(IN) :: Area_M2(:)

    INTEGER  :: I, J, K, Tag_Index, Event_Type
    REAL(f8) :: Event_Mass_Delta, Unfilled_Deficit, Declared_Closure
    REAL(f8) :: Column_Total_Above, Column_Positive_Above
    REAL(f8) :: Column_Net_Mass, Unfillable_Column_Mass
    REAL(f8) :: Q_Before, Q_After

    IF ( .not. Cell_Diagnostics_Enabled_Flag ) RETURN
    IF ( .not. Cell_File_Open .or. .not. ALLOCATED( Cell_Qck_Event_Type ) ) THEN
       CALL Budget_Error( 'Cell diagnostic storage is not initialized', &
                          'Tpcore_Cell_Diag_Write' )
    ENDIF
    IF ( Donor_Diagnostics_Enabled_Flag .and. .not. Donor_File_Open ) THEN
       CALL Budget_Error( 'Donor diagnostic file is not initialized', &
                          'Tpcore_Cell_Diag_Write' )
    ENDIF
    IF ( ANY( SHAPE( Delp2 ) /= (/ Cell_IM, Cell_JM, Cell_KM /) ) .or. &
         SIZE( Area_M2 ) /= Cell_JM ) THEN
       CALL Budget_Error( 'Cell diagnostic geometry differs from its capture', &
                          'Tpcore_Cell_Diag_Write' )
    ENDIF

    DO Tag_Index = 1, N_PLUME_TPCORE_TAGS
    DO K = 1, Cell_KM
    DO J = 1, Cell_JM
    DO I = 1, Cell_IM
       Event_Type = Cell_Qck_Event_Type(I,J,K,Tag_Index)
       IF ( Event_Type == 0 ) CYCLE
       Unfilled_Deficit = 0.0_f8
       Event_Mass_Delta = 0.0_f8
       Declared_Closure = 0.0_f8
       IF ( Event_Type == PLUME_QCK_BOTTOM ) THEN
          Unfilled_Deficit = REAL( Cell_Qck_Deficit(I,J,K,Tag_Index), f8 ) - &
               REAL( Cell_Qck_Withdrawn(I,J,K,Tag_Index), f8 )
          Declared_Closure = REAL( Cell_Qck_Declared_Closure(I,J,K,Tag_Index), f8 )
          Event_Mass_Delta = Declared_Closure * REAL( Area_M2(J), f8 ) * &
               REAL( g0_100, f8 )
       ENDIF
       CALL Write_Cell_Event_Record( Tag_Index, Qck_Event_Type_Name( Event_Type ), &
            I, J, K, Cell_Dq_After_Horizontal(I,J,K,Tag_Index), &
            Cell_Dq_After_Fzppm(I,J,K,Tag_Index), &
            Cell_Qck_At_Correction(I,J,K,Tag_Index), &
            Cell_Qck_Deficit(I,J,K,Tag_Index), &
            Cell_Qck_Available_Above(I,J,K,Tag_Index), &
            Cell_Qck_Withdrawn(I,J,K,Tag_Index), REAL( Unfilled_Deficit, fp ), &
            0.0_fp, 0.0_fp, Delp2(I,J,K), Area_M2(J), Event_Mass_Delta, &
            Qck_Correction_Policy_Name( Event_Type ), &
            Qck_Correction_Outcome_Name( Event_Type, &
                 Cell_Qck_Correction_Status(I,J,K,Tag_Index) ), &
            Cell_Qck_Donor_Count(I,J,K,Tag_Index), &
            Cell_Qck_Full_Column_Withdrawn(I,J,K,Tag_Index), &
            Cell_Qck_Declared_Closure(I,J,K,Tag_Index), &
            Cell_Qck_Correction_Tolerance(I,J,K,Tag_Index), &
            Cell_Qck_Microclosure_Event_Max_Kg(I,J,K,Tag_Index), &
            Cell_Qck_Microclosure_Call_Max_Kg(I,J,K,Tag_Index) )
       IF ( Event_Type == PLUME_QCK_BOTTOM .and. &
            Donor_Diagnostics_Enabled_Flag ) THEN
          Column_Total_Above = REAL( Cell_Qck_Column_Total_Above(I,J,K,Tag_Index), f8 )
          Column_Positive_Above = REAL( &
               Cell_Qck_Column_Positive_Above(I,J,K,Tag_Index), f8 )
          Column_Net_Mass = Column_Total_Above - &
               REAL( Cell_Qck_Deficit(I,J,K,Tag_Index), f8 )
          Unfillable_Column_Mass = MAX( -Column_Net_Mass, 0.0_f8 )
          CALL Write_Donor_Event_Record( Tag_Index, I, J, K, &
               REAL( Cell_Qck_Deficit(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Available_Above(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Withdrawn(I,J,K,Tag_Index), f8 ), &
               Column_Total_Above, Column_Positive_Above, Column_Net_Mass, &
               Unfilled_Deficit, Unfillable_Column_Mass, REAL( Area_M2(J), f8 ), &
               Event_Mass_Delta, Qck_Correction_Policy_Name( Event_Type ), &
               Qck_Correction_Outcome_Name( Event_Type, &
                    Cell_Qck_Correction_Status(I,J,K,Tag_Index) ), &
               Cell_Qck_Donor_Count(I,J,K,Tag_Index), &
               REAL( Cell_Qck_Full_Column_Withdrawn(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Declared_Closure(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Correction_Tolerance(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Microclosure_Event_Max_Kg(I,J,K,Tag_Index), f8 ), &
               REAL( Cell_Qck_Microclosure_Call_Max_Kg(I,J,K,Tag_Index), f8 ) )
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    DO Tag_Index = 1, N_PLUME_TPCORE_TAGS
    DO K = 1, Cell_KM
    DO J = 1, Cell_JM
    DO I = 1, Cell_IM
       Q_Before = REAL( Cell_Qptr_Pre_Floor(I,J,K,Tag_Index), f8 )
       IF ( Q_Before >= 0.0_f8 ) CYCLE
       Q_After = REAL( Cell_Qptr_Post_Floor(I,J,K,Tag_Index), f8 )
       Event_Mass_Delta = ( Q_After - Q_Before ) * REAL( Delp2(I,J,K), f8 ) * &
            REAL( Area_M2(J), f8 ) * REAL( g0_100, f8 )
       CALL Write_Cell_Event_Record( Tag_Index, 'QPTR_NEGATIVE_FLOOR', I, J, K, &
            0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, &
            REAL( Q_Before, fp ), REAL( Q_After, fp ), Delp2(I,J,K), &
            Area_M2(J), Event_Mass_Delta, 'NOT_APPLICABLE', 'NOT_APPLICABLE', &
            0, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp, 0.0_fp )
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    FLUSH( Cell_Unit )
    IF ( Donor_Diagnostics_Enabled_Flag ) FLUSH( Donor_Unit )

  END SUBROUTINE Tpcore_Cell_Diag_Write

  CHARACTER(LEN=16) FUNCTION Qck_Event_Type_Name( Event_Type )

    INTEGER, INTENT(IN) :: Event_Type

    SELECT CASE ( Event_Type )
       CASE ( PLUME_QCK_TOP )
          Qck_Event_Type_Name = 'QCK_TOP'
       CASE ( PLUME_QCK_INTERIOR )
          Qck_Event_Type_Name = 'QCK_INTERIOR'
       CASE ( PLUME_QCK_BOTTOM )
          Qck_Event_Type_Name = 'QCK_BOTTOM'
       CASE DEFAULT
          CALL Budget_Error( 'Unknown Qckxyz event type', 'Qck_Event_Type_Name' )
    END SELECT

  END FUNCTION Qck_Event_Type_Name

  CHARACTER(LEN=48) FUNCTION Qck_Correction_Policy_Name( Event_Type )

    INTEGER, INTENT(IN) :: Event_Type

    SELECT CASE ( Event_Type )
       CASE ( PLUME_QCK_BOTTOM )
          Qck_Correction_Policy_Name = &
               'full_column_nearest_above_microclosure_v2'
       CASE ( PLUME_QCK_TOP, PLUME_QCK_INTERIOR )
          Qck_Correction_Policy_Name = 'native_qck_v1'
       CASE DEFAULT
          Qck_Correction_Policy_Name = 'not_applicable'
    END SELECT

  END FUNCTION Qck_Correction_Policy_Name

  CHARACTER(LEN=24) FUNCTION Qck_Correction_Outcome_Name( Event_Type, &
                                                            Correction_Status )

    INTEGER, INTENT(IN) :: Event_Type, Correction_Status

    IF ( Event_Type /= PLUME_QCK_BOTTOM ) THEN
       Qck_Correction_Outcome_Name = 'not_applicable'
       RETURN
    ENDIF
    SELECT CASE ( Correction_Status )
       CASE ( QCK_BOTTOM_EXACT )
          Qck_Correction_Outcome_Name = 'exact'
       CASE ( QCK_BOTTOM_ROUNDOFF )
          Qck_Correction_Outcome_Name = 'roundoff_closed'
       CASE ( QCK_BOTTOM_MICROCLOSURE )
          Qck_Correction_Outcome_Name = 'microclosure_closed'
       CASE DEFAULT
          CALL Budget_Error( 'Unknown QCK_BOTTOM correction outcome', &
                             'Qck_Correction_Outcome_Name' )
    END SELECT

  END FUNCTION Qck_Correction_Outcome_Name

  SUBROUTINE Write_Cell_Event_Record( Tag_Index, Event_Type, I, J, K, &
                                      Dq_After_Horizontal, Dq_After_Fzppm, &
                                      Dq_At_Correction, Deficit, Available_Above, &
                                      Withdrawn, Unfilled_Deficit, Q_Before, &
                                      Q_After, Delp, Area, Event_Mass_Delta, &
                                      Correction_Policy, Correction_Outcome, &
                                      Donor_Count, Full_Column_Withdrawn, &
                                      Declared_Closure, Correction_Tolerance, &
                                      Microclosure_Event_Max_Kg, &
                                      Microclosure_Call_Max_Kg )

    INTEGER,          INTENT(IN) :: Tag_Index, I, J, K
    CHARACTER(LEN=*), INTENT(IN) :: Event_Type
    REAL(fp),         INTENT(IN) :: Dq_After_Horizontal, Dq_After_Fzppm
    REAL(fp),         INTENT(IN) :: Dq_At_Correction, Deficit, Available_Above
    REAL(fp),         INTENT(IN) :: Withdrawn, Unfilled_Deficit, Q_Before, Q_After
    REAL(fp),         INTENT(IN) :: Delp, Area
    REAL(f8),         INTENT(IN) :: Event_Mass_Delta
    CHARACTER(LEN=*), INTENT(IN) :: Correction_Policy, Correction_Outcome
    INTEGER,          INTENT(IN) :: Donor_Count
    REAL(fp),         INTENT(IN) :: Full_Column_Withdrawn, Declared_Closure
    REAL(fp),         INTENT(IN) :: Correction_Tolerance
    REAL(fp),         INTENT(IN) :: Microclosure_Event_Max_Kg
    REAL(fp),         INTENT(IN) :: Microclosure_Call_Max_Kg

    WRITE( Cell_Unit, &
           '(3(a,","),6(i0,","),2(a,","),4(i0,","),12(es26.17e3,","),' // &
           '2(a,","),i0,",",4(es26.17e3,","),es26.17e3)' ) &
         'plume-tpcore-cell-events-v3', TRIM( Budget_Run_Id ), &
         TRIM( Budget_Manifest_Id ), Model_Date, Model_Time, Elapsed_Seconds, &
         Heartbeat_Index, Omp_Thread_Count, Call_Index, TRIM( TAG_NAMES(Tag_Index) ), &
         TRIM( Event_Type ), I, J, K, Cell_KM + 1 - K, &
         Dq_After_Horizontal, Dq_After_Fzppm, Dq_At_Correction, Deficit, &
         Available_Above, Withdrawn, Unfilled_Deficit, Q_Before, Q_After, &
         Delp, Area, Event_Mass_Delta, TRIM( Correction_Policy ), &
         TRIM( Correction_Outcome ), Donor_Count, Full_Column_Withdrawn, &
         Declared_Closure, Correction_Tolerance, Microclosure_Event_Max_Kg, &
         Microclosure_Call_Max_Kg

  END SUBROUTINE Write_Cell_Event_Record

  SUBROUTINE Write_Donor_Event_Record( Tag_Index, I, J, K, Deficit, &
                                       Immediate_Available, Withdrawn, &
                                       Column_Total_Above, Column_Positive_Above, &
                                       Column_Net_Mass, Unfilled_Immediate, &
                                       Unfillable_Column_Mass, Area, &
                                       Event_Mass_Delta, Correction_Policy, &
                                       Correction_Outcome, Donor_Count, &
                                       Full_Column_Withdrawn, Declared_Closure, &
                                       Correction_Tolerance, &
                                       Microclosure_Event_Max_Kg, &
                                       Microclosure_Call_Max_Kg )

    INTEGER,  INTENT(IN) :: Tag_Index, I, J, K
    REAL(f8), INTENT(IN) :: Deficit, Immediate_Available, Withdrawn
    REAL(f8), INTENT(IN) :: Column_Total_Above, Column_Positive_Above
    REAL(f8), INTENT(IN) :: Column_Net_Mass, Unfilled_Immediate
    REAL(f8), INTENT(IN) :: Unfillable_Column_Mass, Area, Event_Mass_Delta
    CHARACTER(LEN=*), INTENT(IN) :: Correction_Policy, Correction_Outcome
    INTEGER,          INTENT(IN) :: Donor_Count
    REAL(f8),         INTENT(IN) :: Full_Column_Withdrawn, Declared_Closure
    REAL(f8),         INTENT(IN) :: Correction_Tolerance
    REAL(f8),         INTENT(IN) :: Microclosure_Event_Max_Kg
    REAL(f8),         INTENT(IN) :: Microclosure_Call_Max_Kg

    WRITE( Donor_Unit, &
           '(3(a,","),6(i0,","),a,",",4(i0,","),10(es26.17e3,","),' // &
           '2(a,","),i0,",",4(es26.17e3,","),es26.17e3)' ) &
         'plume-tpcore-donor-events-v3', TRIM( Budget_Run_Id ), &
         TRIM( Budget_Manifest_Id ), Model_Date, Model_Time, Elapsed_Seconds, &
         Heartbeat_Index, Omp_Thread_Count, Call_Index, TRIM( TAG_NAMES(Tag_Index) ), &
         I, J, K, Cell_KM + 1 - K, Deficit, Immediate_Available, Withdrawn, &
         Column_Total_Above, Column_Positive_Above, Column_Net_Mass, &
         Unfilled_Immediate, Unfillable_Column_Mass, Area, Event_Mass_Delta, &
         TRIM( Correction_Policy ), TRIM( Correction_Outcome ), Donor_Count, &
         Full_Column_Withdrawn, Declared_Closure, Correction_Tolerance, &
         Microclosure_Event_Max_Kg, Microclosure_Call_Max_Kg

  END SUBROUTINE Write_Donor_Event_Record

  SUBROUTINE Tpcore_Budget_Begin( State_Chm, Ak, Bk, Ps1, Area_M2, &
                                  N_Advect, Transport_Timestep_S, &
                                  Qckxyz_Was_Invoked )

    TYPE(ChmState), INTENT(IN) :: State_Chm
    REAL(fp),       INTENT(IN) :: Ak(:), Bk(:)
    REAL(fp),       INTENT(IN) :: Ps1(:,:)
    REAL(fp),       INTENT(IN) :: Area_M2(:)
    INTEGER,        INTENT(IN) :: N_Advect
    REAL(fp),       INTENT(IN) :: Transport_Timestep_S
    LOGICAL,        INTENT(IN) :: Qckxyz_Was_Invoked

    INTEGER  :: I
    REAL(f8) :: Global_Dry_Air_Mass
    TYPE(Plume_Tpcore_Metrics) :: Entry_Metrics

    IF ( .not. Tpcore_Budget_Enabled() ) RETURN
    IF ( Call_Active ) THEN
       CALL Budget_Error( 'TPCORE budget call began before the prior call ended', &
                          'Tpcore_Budget_Begin' )
    ENDIF
    IF ( SIZE( Ps1, 2 ) /= SIZE( Area_M2 ) ) THEN
       CALL Budget_Error( 'TPCORE budget latitude extent does not match area', &
                          'Tpcore_Budget_Begin' )
    ENDIF
    IF ( SIZE( Ak ) /= SIZE( Bk ) .or. SIZE( Ak ) < 2 ) THEN
       CALL Budget_Error( 'TPCORE budget hybrid-coordinate arrays are invalid', &
                          'Tpcore_Budget_Begin' )
    ENDIF
    IF ( Transport_Timestep_S <= 0.0_fp ) THEN
       CALL Budget_Error( 'TPCORE budget received a non-positive timestep', &
                          'Tpcore_Budget_Begin' )
    ENDIF

    DO I = 1, N_PLUME_TPCORE_TAGS
       Tag_Ids(I) = Ind_( TRIM( TAG_NAMES(I) ), 'S' )
       IF ( Tag_Ids(I) < 1 .or. Tag_Ids(I) > N_Advect ) THEN
          CALL Budget_Error( 'Missing or non-advected plume tag: ' // &
                             TRIM( TAG_NAMES(I) ), 'Tpcore_Budget_Begin' )
       ENDIF
       IF ( State_Chm%Species(Tag_Ids(I))%Units /= &
            KG_SPECIES_PER_KG_DRY_AIR ) THEN
          CALL Budget_Error( 'Plume tag is not in kg/kg dry-air units: ' // &
                             TRIM( TAG_NAMES(I) ), 'Tpcore_Budget_Begin' )
       ENDIF
       IF ( .not. ASSOCIATED( State_Chm%Species(Tag_Ids(I))%Conc ) ) THEN
          CALL Budget_Error( 'Plume tag concentration is not associated: ' // &
                             TRIM( TAG_NAMES(I) ), 'Tpcore_Budget_Begin' )
       ENDIF
    ENDDO

    Call_Index       = Call_Index + 1
    Call_Active      = .TRUE.
    Elapsed_Seconds  = GET_ELAPSED_SEC()
    CALL Canonical_Diagnostic_Time( GET_NYMDb(), GET_NHMSb(), &
                                    Elapsed_Seconds, Model_Date, Model_Time )
    Heartbeat_Index  = INT( REAL( Elapsed_Seconds, f8 ) / &
                             REAL( Transport_Timestep_S, f8 ) ) + 1
    Omp_Thread_Count = OMP_GET_MAX_THREADS()
    Qckxyz_Invoked   = 0
    IF ( Qckxyz_Was_Invoked ) Qckxyz_Invoked = 1
    Previous_Mass    = 0.0_f8
    Global_Dry_Air_Mass = Dry_Air_Mass_From_Pressure( Ak, Bk, Ps1, Area_M2 )

    DO I = 1, N_PLUME_TPCORE_TAGS
       Entry_Metrics%Global_Mass_Kg = Conc_Mass_From_Pressure( &
            State_Chm%Species(Tag_Ids(I))%Conc, Ak, Bk, Ps1, Area_M2 )
       Entry_Metrics%Negative_Cell_Count = 0
       Entry_Metrics%Minimum_Mixing_Ratio = 0.0_f8
       Entry_Metrics%Total_Negative_Mass_Kg = 0.0_f8
       CALL Write_Boundary_Record( 'TPCORE_ENTRY', I, Entry_Metrics, &
                                   Global_Dry_Air_Mass, 0, 0.0_f8 )
    ENDDO

  END SUBROUTINE Tpcore_Budget_Begin

  INTEGER FUNCTION Tpcore_Budget_Tag_Index( Species_Id )

    INTEGER, INTENT(IN) :: Species_Id
    INTEGER :: I

    Tpcore_Budget_Tag_Index = 0
    IF ( .not. Budgets_Enabled ) RETURN
    DO I = 1, N_PLUME_TPCORE_TAGS
       IF ( Tag_Ids(I) == Species_Id ) THEN
          Tpcore_Budget_Tag_Index = I
          RETURN
       ENDIF
    ENDDO

  END FUNCTION Tpcore_Budget_Tag_Index

  SUBROUTINE Compute_Dq_Metrics( Dq, Delp2, Area_M2, Metrics )

    REAL(fp),                  INTENT(IN)  :: Dq(:,:,:)
    REAL(fp),                  INTENT(IN)  :: Delp2(:,:,:)
    REAL(fp),                  INTENT(IN)  :: Area_M2(:)
    TYPE(Plume_Tpcore_Metrics), INTENT(OUT) :: Metrics

    INTEGER  :: I, J, K
    REAL(f8) :: Dry_Air_Mass, Mixing_Ratio, Tracer_Mass

    IF ( ANY( SHAPE( Dq ) /= SHAPE( Delp2 ) ) ) THEN
       CALL Budget_Error( 'DQ and DELP2 shapes differ', 'Compute_Dq_Metrics' )
    ENDIF
    IF ( SIZE( Dq, 2 ) /= SIZE( Area_M2 ) ) THEN
       CALL Budget_Error( 'DQ latitude extent does not match area', &
                          'Compute_Dq_Metrics' )
    ENDIF

    Metrics%Global_Mass_Kg        = 0.0_f8
    Metrics%Negative_Cell_Count   = 0
    Metrics%Minimum_Mixing_Ratio  = HUGE( 0.0_f8 )
    Metrics%Total_Negative_Mass_Kg = 0.0_f8

    DO K = 1, SIZE( Dq, 3 )
    DO J = 1, SIZE( Dq, 2 )
       Dry_Air_Mass = REAL( Area_M2(J), f8 ) * REAL( g0_100, f8 )
       DO I = 1, SIZE( Dq, 1 )
          Tracer_Mass = REAL( Dq(I,J,K), f8 ) * Dry_Air_Mass
          IF ( Delp2(I,J,K) <= 0.0_fp ) THEN
             CALL Budget_Error( 'Non-positive DELP2 in Qckxyz metric', &
                                'Compute_Dq_Metrics' )
          ENDIF
          Mixing_Ratio = REAL( Dq(I,J,K), f8 ) / REAL( Delp2(I,J,K), f8 )
          Metrics%Global_Mass_Kg = Metrics%Global_Mass_Kg + Tracer_Mass
          Metrics%Minimum_Mixing_Ratio = MIN( Metrics%Minimum_Mixing_Ratio, &
                                               Mixing_Ratio )
          IF ( Tracer_Mass < 0.0_f8 ) THEN
             Metrics%Negative_Cell_Count = Metrics%Negative_Cell_Count + 1
             Metrics%Total_Negative_Mass_Kg = &
                  Metrics%Total_Negative_Mass_Kg - Tracer_Mass
          ENDIF
       ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE Compute_Dq_Metrics

  SUBROUTINE Compute_Conc_Metrics( Conc, Delp2, Area_M2, Metrics )

    REAL(fp),                   INTENT(IN)  :: Conc(:,:,:)
    REAL(fp),                   INTENT(IN)  :: Delp2(:,:,:)
    REAL(fp),                   INTENT(IN)  :: Area_M2(:)
    TYPE(Plume_Tpcore_Metrics), INTENT(OUT) :: Metrics

    INTEGER  :: I, J, K
    REAL(f8) :: Dry_Air_Mass, Mixing_Ratio, Tracer_Mass

    IF ( ANY( SHAPE( Conc ) /= SHAPE( Delp2 ) ) ) THEN
       CALL Budget_Error( 'Concentration and DELP2 shapes differ', &
                          'Compute_Conc_Metrics' )
    ENDIF
    IF ( SIZE( Conc, 2 ) /= SIZE( Area_M2 ) ) THEN
       CALL Budget_Error( 'Concentration latitude extent does not match area', &
                          'Compute_Conc_Metrics' )
    ENDIF

    Metrics%Global_Mass_Kg         = 0.0_f8
    Metrics%Negative_Cell_Count    = 0
    Metrics%Minimum_Mixing_Ratio   = HUGE( 0.0_f8 )
    Metrics%Total_Negative_Mass_Kg = 0.0_f8

    DO K = 1, SIZE( Conc, 3 )
    DO J = 1, SIZE( Conc, 2 )
       Dry_Air_Mass = REAL( Area_M2(J), f8 ) * REAL( g0_100, f8 )
       DO I = 1, SIZE( Conc, 1 )
          IF ( Delp2(I,J,K) <= 0.0_fp ) THEN
             CALL Budget_Error( 'Non-positive DELP2 in concentration metric', &
                                'Compute_Conc_Metrics' )
          ENDIF
          Mixing_Ratio = REAL( Conc(I,J,K), f8 )
          Tracer_Mass = Mixing_Ratio * REAL( Delp2(I,J,K), f8 ) * &
                        Dry_Air_Mass
          Metrics%Global_Mass_Kg = Metrics%Global_Mass_Kg + Tracer_Mass
          Metrics%Minimum_Mixing_Ratio = MIN( Metrics%Minimum_Mixing_Ratio, &
                                               Mixing_Ratio )
          IF ( Mixing_Ratio < 0.0_f8 ) THEN
             Metrics%Negative_Cell_Count = Metrics%Negative_Cell_Count + 1
             Metrics%Total_Negative_Mass_Kg = &
                  Metrics%Total_Negative_Mass_Kg - Tracer_Mass
          ENDIF
       ENDDO
    ENDDO
    ENDDO

  END SUBROUTINE Compute_Conc_Metrics

  SUBROUTINE Tpcore_Budget_Write_Positivity_Boundaries( Qck_Pre_Metrics, &
                                                         Qck_Post_Metrics, &
                                                         Qck_Corrected_Cells, &
                                                         Pre_Floor_Metrics, &
                                                         Post_Floor_Metrics, &
                                                         Delp2, Area_M2 )

    TYPE(Plume_Tpcore_Metrics), INTENT(IN) :: &
         Qck_Pre_Metrics(N_PLUME_TPCORE_TAGS)
    TYPE(Plume_Tpcore_Metrics), INTENT(IN) :: &
         Qck_Post_Metrics(N_PLUME_TPCORE_TAGS)
    INTEGER, INTENT(IN) :: Qck_Corrected_Cells(N_PLUME_TPCORE_TAGS)
    TYPE(Plume_Tpcore_Metrics), INTENT(IN) :: &
         Pre_Floor_Metrics(N_PLUME_TPCORE_TAGS)
    TYPE(Plume_Tpcore_Metrics), INTENT(IN) :: &
         Post_Floor_Metrics(N_PLUME_TPCORE_TAGS)
    REAL(fp),                    INTENT(IN) :: Delp2(:,:,:)
    REAL(fp),                    INTENT(IN) :: Area_M2(:)

    INTEGER  :: I
    REAL(f8) :: Global_Dry_Air_Mass, Correction_Mass_Delta

    IF ( .not. Budgets_Enabled ) RETURN
    IF ( .not. Call_Active ) THEN
       CALL Budget_Error( 'Positivity boundary written without active TPCORE call', &
                          'Tpcore_Budget_Write_Positivity_Boundaries' )
    ENDIF
    Global_Dry_Air_Mass = Dry_Air_Mass_From_Delp( Delp2, Area_M2 )
    DO I = 1, N_PLUME_TPCORE_TAGS
       CALL Write_Boundary_Record( 'PRE_QCKXYZ', I, Qck_Pre_Metrics(I), &
                                   Global_Dry_Air_Mass, 0, 0.0_f8 )
       Correction_Mass_Delta = Qck_Post_Metrics(I)%Global_Mass_Kg - &
                               Qck_Pre_Metrics(I)%Global_Mass_Kg
       CALL Write_Boundary_Record( 'POST_QCKXYZ', I, Qck_Post_Metrics(I), &
                                   Global_Dry_Air_Mass, Qck_Corrected_Cells(I), &
                                   Correction_Mass_Delta )
       CALL Write_Boundary_Record( 'PRE_QPTR_NEGATIVE_FLOOR', I, &
                                   Pre_Floor_Metrics(I), &
                                   Global_Dry_Air_Mass, 0, 0.0_f8 )
       Correction_Mass_Delta = Post_Floor_Metrics(I)%Global_Mass_Kg - &
                               Pre_Floor_Metrics(I)%Global_Mass_Kg
       CALL Write_Boundary_Record( 'POST_QPTR_NEGATIVE_FLOOR', I, &
                                   Post_Floor_Metrics(I), &
                                   Global_Dry_Air_Mass, &
                                   Pre_Floor_Metrics(I)%Negative_Cell_Count, &
                                   Correction_Mass_Delta )
    ENDDO

  END SUBROUTINE Tpcore_Budget_Write_Positivity_Boundaries

  SUBROUTINE Tpcore_Budget_End( State_Chm, State_Met )

    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(MetState), INTENT(IN) :: State_Met

    INTEGER  :: I, J, K, L
    REAL(f8) :: Global_Dry_Air_Mass
    TYPE(Plume_Tpcore_Metrics) :: Exit_Metrics

    IF ( .not. Budgets_Enabled ) RETURN
    IF ( .not. Call_Active ) THEN
       CALL Budget_Error( 'TPCORE budget ended without an active call', &
                          'Tpcore_Budget_End' )
    ENDIF
    IF ( .not. ASSOCIATED( State_Met%AD ) ) THEN
       CALL Budget_Error( 'State_Met%AD is not associated at TPCORE exit', &
                          'Tpcore_Budget_End' )
    ENDIF

    Global_Dry_Air_Mass = 0.0_f8
    DO L = 1, SIZE( State_Met%AD, 3 )
    DO J = 1, SIZE( State_Met%AD, 2 )
    DO I = 1, SIZE( State_Met%AD, 1 )
       Global_Dry_Air_Mass = Global_Dry_Air_Mass + REAL( State_Met%AD(I,J,L), f8 )
    ENDDO
    ENDDO
    ENDDO

    DO I = 1, N_PLUME_TPCORE_TAGS
       Exit_Metrics%Global_Mass_Kg = 0.0_f8
       DO L = 1, SIZE( State_Met%AD, 3 )
       DO J = 1, SIZE( State_Met%AD, 2 )
       DO K = 1, SIZE( State_Met%AD, 1 )
          Exit_Metrics%Global_Mass_Kg = Exit_Metrics%Global_Mass_Kg + &
               REAL( State_Chm%Species(Tag_Ids(I))%Conc(K,J,L), f8 ) * &
               REAL( State_Met%AD(K,J,L), f8 )
       ENDDO
       ENDDO
       ENDDO
       Exit_Metrics%Negative_Cell_Count = 0
       Exit_Metrics%Minimum_Mixing_Ratio = 0.0_f8
       Exit_Metrics%Total_Negative_Mass_Kg = 0.0_f8
       CALL Write_Boundary_Record( 'TPCORE_EXIT_POST_RESET', I, &
                                   Exit_Metrics, Global_Dry_Air_Mass, &
                                   0, 0.0_f8 )
    ENDDO
    Call_Active = .FALSE.

  END SUBROUTINE Tpcore_Budget_End

  SUBROUTINE Write_Boundary_Record( Boundary, Tag_Index, Metrics, &
                                    Global_Dry_Air_Mass, Corrected_Cell_Count, &
                                    Correction_Mass_Delta )

    CHARACTER(LEN=*),             INTENT(IN) :: Boundary
    INTEGER,                      INTENT(IN) :: Tag_Index
    TYPE(Plume_Tpcore_Metrics),   INTENT(IN) :: Metrics
    REAL(f8),                     INTENT(IN) :: Global_Dry_Air_Mass
    INTEGER,                      INTENT(IN) :: Corrected_Cell_Count
    REAL(f8),                     INTENT(IN) :: Correction_Mass_Delta

    REAL(f8) :: Boundary_Mass_Delta

    IF ( .not. File_Open ) THEN
       CALL Budget_Error( 'TPCORE budget file is not open', 'Write_Boundary_Record' )
    ENDIF
    IF ( Tag_Index < 1 .or. Tag_Index > N_PLUME_TPCORE_TAGS ) THEN
       CALL Budget_Error( 'Invalid plume-tag index', 'Write_Boundary_Record' )
    ENDIF

    Boundary_Mass_Delta = Metrics%Global_Mass_Kg - Previous_Mass(Tag_Index)
    IF ( TRIM( Boundary ) == 'TPCORE_ENTRY' ) Boundary_Mass_Delta = 0.0_f8
    WRITE( Budget_Unit, &
           '(a,",",a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",' // &
           'a,",",a,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",i0,",",' // &
           'es26.17e3,",",es26.17e3,",",i0,",",es26.17e3,",",i0)' ) &
         'plume-tpcore-budget-v2', TRIM( Budget_Run_Id ), &
         TRIM( Budget_Manifest_Id ), Model_Date, Model_Time, Elapsed_Seconds, &
         Heartbeat_Index, Omp_Thread_Count, Call_Index, TRIM( Boundary ), &
         TRIM( TAG_NAMES(Tag_Index) ), Metrics%Global_Mass_Kg, &
         Boundary_Mass_Delta, Global_Dry_Air_Mass, Metrics%Negative_Cell_Count, &
         Metrics%Minimum_Mixing_Ratio, Metrics%Total_Negative_Mass_Kg, &
         Corrected_Cell_Count, Correction_Mass_Delta, Qckxyz_Invoked
    FLUSH( Budget_Unit )
    Previous_Mass(Tag_Index) = Metrics%Global_Mass_Kg

  END SUBROUTINE Write_Boundary_Record

  REAL(f8) FUNCTION Dry_Air_Mass_From_Delp( Delp, Area_M2 )

    REAL(fp), INTENT(IN) :: Delp(:,:,:)
    REAL(fp), INTENT(IN) :: Area_M2(:)

    INTEGER :: I, J, K

    IF ( SIZE( Delp, 2 ) /= SIZE( Area_M2 ) ) THEN
       CALL Budget_Error( 'DELP latitude extent does not match area', &
                          'Dry_Air_Mass_From_Delp' )
    ENDIF
    Dry_Air_Mass_From_Delp = 0.0_f8
    DO K = 1, SIZE( Delp, 3 )
    DO J = 1, SIZE( Delp, 2 )
    DO I = 1, SIZE( Delp, 1 )
       Dry_Air_Mass_From_Delp = Dry_Air_Mass_From_Delp + &
            REAL( Delp(I,J,K), f8 ) * REAL( Area_M2(J), f8 ) * &
            REAL( g0_100, f8 )
    ENDDO
    ENDDO
    ENDDO

  END FUNCTION Dry_Air_Mass_From_Delp

  REAL(f8) FUNCTION Dry_Air_Mass_From_Pressure( Ak, Bk, Ps1, Area_M2 )

    REAL(fp), INTENT(IN) :: Ak(:), Bk(:)
    REAL(fp), INTENT(IN) :: Ps1(:,:)
    REAL(fp), INTENT(IN) :: Area_M2(:)

    INTEGER :: I, J, K
    REAL(f8) :: Delp

    IF ( SIZE( Ps1, 2 ) /= SIZE( Area_M2 ) .or. &
         SIZE( Ak ) /= SIZE( Bk ) .or. SIZE( Ak ) < 2 ) THEN
       CALL Budget_Error( 'Invalid entry-pressure geometry', &
                          'Dry_Air_Mass_From_Pressure' )
    ENDIF
    Dry_Air_Mass_From_Pressure = 0.0_f8
    DO K = 1, SIZE( Ak ) - 1
       DO J = 1, SIZE( Ps1, 2 )
       DO I = 1, SIZE( Ps1, 1 )
          Delp = REAL( Ak(K+1) - Ak(K), f8 ) + &
                 REAL( Bk(K+1) - Bk(K), f8 ) * REAL( Ps1(I,J), f8 )
          Dry_Air_Mass_From_Pressure = Dry_Air_Mass_From_Pressure + &
               Delp * REAL( Area_M2(J), f8 ) * REAL( g0_100, f8 )
       ENDDO
       ENDDO
    ENDDO

  END FUNCTION Dry_Air_Mass_From_Pressure

  REAL(f8) FUNCTION Conc_Mass_From_Pressure( Conc, Ak, Bk, Ps1, Area_M2 )

    REAL(fp), INTENT(IN) :: Conc(:,:,:)
    REAL(fp), INTENT(IN) :: Ak(:), Bk(:)
    REAL(fp), INTENT(IN) :: Ps1(:,:)
    REAL(fp), INTENT(IN) :: Area_M2(:)

    INTEGER :: I, J, K, K_Conc
    REAL(f8) :: Delp

    IF ( SIZE( Conc, 1 ) /= SIZE( Ps1, 1 ) .or. &
         SIZE( Conc, 2 ) /= SIZE( Ps1, 2 ) .or. &
         SIZE( Conc, 3 ) /= SIZE( Ak ) - 1 .or. &
         SIZE( Ak ) /= SIZE( Bk ) .or. SIZE( Ps1, 2 ) /= SIZE( Area_M2 ) ) THEN
       CALL Budget_Error( 'Invalid concentration/entry-pressure geometry', &
                          'Conc_Mass_From_Pressure' )
    ENDIF
    Conc_Mass_From_Pressure = 0.0_f8
    DO K = 1, SIZE( Ak ) - 1
       K_Conc = SIZE( Ak ) - K
       DO J = 1, SIZE( Ps1, 2 )
       DO I = 1, SIZE( Ps1, 1 )
          Delp = REAL( Ak(K+1) - Ak(K), f8 ) + &
                 REAL( Bk(K+1) - Bk(K), f8 ) * REAL( Ps1(I,J), f8 )
          Conc_Mass_From_Pressure = Conc_Mass_From_Pressure + &
               REAL( Conc(I,J,K_Conc), f8 ) * Delp * &
               REAL( Area_M2(J), f8 ) * REAL( g0_100, f8 )
       ENDDO
       ENDDO
    ENDDO

  END FUNCTION Conc_Mass_From_Pressure

  SUBROUTINE Budget_Error( Message, Routine )

    CHARACTER(LEN=*), INTENT(IN) :: Message, Routine
    CHARACTER(LEN=255) :: ThisLoc

    ThisLoc = ' -> at ' // TRIM( Routine ) // &
              ' (in plume_tpcore_budget_mod.F90)'
    CALL Error_Stop( Message, ThisLoc )

  END SUBROUTINE Budget_Error

END MODULE Plume_Tpcore_Budget_Mod
