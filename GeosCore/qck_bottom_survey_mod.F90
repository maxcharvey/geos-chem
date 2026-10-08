!------------------------------------------------------------------------------
! Runtime-gated survey ledger for native TPCORE QCK_BOTTOM behavior.
!
! This module is diagnostic-only.  When GC_QCK_BOTTOM_SURVEY=1, Qckxyz uses
! the historical native bottom correction while separately preflighting the
! conservative full-column policy.  Every negative bottom event is written in
! a deterministic spatial order within each species.  The default-off path is
! unchanged.
!------------------------------------------------------------------------------
MODULE Qck_Bottom_Survey_Mod

  USE ERROR_MOD,          ONLY : Error_Stop
  USE OMP_LIB,            ONLY : OMP_GET_MAX_THREADS
  USE PhysConstants,      ONLY : g0_100
  USE Precision_Mod,      ONLY : fp, f8
  USE Qck_Positivity_Mod, ONLY : QCK_BOTTOM_EXACT, QCK_BOTTOM_ROUNDOFF, &
                                  QCK_BOTTOM_UNFILLABLE,                &
                                  QCK_BOTTOM_INVALID_DONOR
  USE TIME_MOD,           ONLY : GET_ELAPSED_SEC, GET_NHMS, GET_NYMD

  IMPLICIT NONE
  PRIVATE

  LOGICAL, SAVE :: Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Survey_Enabled_Flag = .FALSE.
  LOGICAL, SAVE :: File_Open = .FALSE.
  INTEGER, SAVE :: Survey_Unit = -1
  CHARACTER(LEN=128), SAVE :: Survey_Run_Id = ''

  PUBLIC :: Qck_Bottom_Survey_Enabled
  PUBLIC :: Qck_Bottom_Survey_Write

CONTAINS

  LOGICAL FUNCTION Qck_Bottom_Survey_Enabled()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( .not. Gate_Initialized ) THEN
       Gate_Initialized = .TRUE.
       Value = ''
       CALL GET_ENVIRONMENT_VARIABLE( 'GC_QCK_BOTTOM_SURVEY', Value, &
                                      STATUS=Env_Status )
       IF ( Env_Status == 0 ) THEN
          SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
             CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
                Survey_Enabled_Flag = .TRUE.
             CASE DEFAULT
                Survey_Enabled_Flag = .FALSE.
          END SELECT
       ENDIF

       IF ( Survey_Enabled_Flag ) THEN
          ThisLoc = ' -> at Qck_Bottom_Survey_Enabled ' // &
                    '(in qck_bottom_survey_mod.F90)'
          Value = ''
          CALL GET_ENVIRONMENT_VARIABLE( 'GC_QCK_BOTTOM_SURVEY_RUN_ID', &
                                         Value, STATUS=Env_Status )
          IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
             ErrMsg = 'GC_QCK_BOTTOM_SURVEY_RUN_ID is required when survey is enabled'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
          Survey_Run_Id = TRIM( Value )

          Value = ''
          CALL GET_ENVIRONMENT_VARIABLE( 'GC_QCK_BOTTOM_SURVEY_FILE', Value, &
                                         STATUS=Env_Status )
          IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
             ErrMsg = 'GC_QCK_BOTTOM_SURVEY_FILE is required when survey is enabled'
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
          OPEN( NEWUNIT=Survey_Unit, FILE=TRIM( Value ), STATUS='NEW', &
                ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
          IF ( IO_Status /= 0 ) THEN
             ErrMsg = 'Could not create QCK_BOTTOM survey ledger: ' // TRIM( Value )
             CALL Error_Stop( ErrMsg, ThisLoc )
          ENDIF
          File_Open = .TRUE.
          WRITE( Survey_Unit, '(a)' ) &
               'schema_version,run_id,model_date,model_time,elapsed_seconds,' // &
               'omp_thread_count,species_index,i_tpcore,j_tpcore,k_tpcore,' // &
               'feasibility_status,deficit_hpa,full_column_available_hpa,' // &
               'shortfall_hpa,event_tolerance_hpa,immediate_donor_hpa,' // &
               'native_withdrawn_hpa,native_mass_created_hpa,' // &
               'deficit_kg,full_column_available_kg,shortfall_kg,' // &
               'native_mass_created_kg,area_m2,behavior'
          FLUSH( Survey_Unit )
          WRITE( 6, '(a)' ) 'QCK_BOTTOM SURVEY ENABLED: ' // TRIM( Value )
       ENDIF
    ENDIF

    Qck_Bottom_Survey_Enabled = Survey_Enabled_Flag

  END FUNCTION Qck_Bottom_Survey_Enabled

  SUBROUTINE Qck_Bottom_Survey_Write( Species_Index, K_Bottom, I_First, &
                                       J_First, Area_M2, Status, Deficit, &
                                       Available, Shortfall, Tolerance, &
                                       Immediate_Donor )

    INTEGER,  INTENT(IN) :: Species_Index, K_Bottom, I_First, J_First
    REAL(fp), INTENT(IN) :: Area_M2(:)
    INTEGER,  INTENT(IN) :: Status(:,:)
    REAL(fp), INTENT(IN) :: Deficit(:,:), Available(:,:), Shortfall(:,:)
    REAL(fp), INTENT(IN) :: Tolerance(:,:)
    REAL(fp), INTENT(IN) :: Immediate_Donor(:,:)

    INTEGER           :: I_Local, J_Local, I_Global, J_Global
    INTEGER           :: Model_Date, Model_Time, Elapsed_Seconds
    INTEGER           :: Omp_Thread_Count
    REAL(f8)          :: Area, Conversion, Immediate_Withdrawn
    REAL(f8)          :: Deficit_Hpa, Available_Hpa, Shortfall_Hpa
    REAL(f8)          :: Tolerance_Hpa, Immediate_Hpa, Native_Created_Hpa
    CHARACTER(LEN=24) :: Status_Name

    IF ( .not. Survey_Enabled_Flag ) RETURN
    IF ( .not. File_Open ) THEN
       CALL Survey_Error( 'Survey ledger is not open', &
                          'Qck_Bottom_Survey_Write' )
    ENDIF
    IF ( SIZE( Status, 1 ) /= SIZE( Deficit, 1 ) .or. &
         SIZE( Status, 2 ) /= SIZE( Deficit, 2 ) .or. &
         ANY( SHAPE( Status ) /= SHAPE( Available ) ) .or. &
         ANY( SHAPE( Status ) /= SHAPE( Shortfall ) ) .or. &
         ANY( SHAPE( Status ) /= SHAPE( Tolerance ) ) .or. &
         ANY( SHAPE( Status ) /= SHAPE( Immediate_Donor ) ) ) THEN
       CALL Survey_Error( 'Survey field shapes differ', &
                          'Qck_Bottom_Survey_Write' )
    ENDIF

    Model_Date      = GET_NYMD()
    Model_Time      = GET_NHMS()
    Elapsed_Seconds = GET_ELAPSED_SEC()
    Omp_Thread_Count = OMP_GET_MAX_THREADS()

    !$OMP CRITICAL (QCK_BOTTOM_SURVEY_FILE)
    DO J_Local = 1, SIZE( Status, 2 )
       J_Global = J_First + J_Local - 1
       IF ( J_Global < 1 .or. J_Global > SIZE( Area_M2 ) ) THEN
          CALL Survey_Error( 'Survey latitude index is outside area array', &
                             'Qck_Bottom_Survey_Write' )
       ENDIF
       Area = REAL( Area_M2(J_Global), f8 )
       Conversion = Area * REAL( g0_100, f8 )
       DO I_Local = 1, SIZE( Status, 1 )
          IF ( Status(I_Local,J_Local) < 0 ) CYCLE
          I_Global = I_First + I_Local - 1
          Deficit_Hpa   = REAL( Deficit(I_Local,J_Local), f8 )
          Available_Hpa = REAL( Available(I_Local,J_Local), f8 )
          Tolerance_Hpa = REAL( Tolerance(I_Local,J_Local), f8 )
          Immediate_Hpa = REAL( Immediate_Donor(I_Local,J_Local), f8 )
          Immediate_Withdrawn = MIN( Deficit_Hpa, Immediate_Hpa )
          Shortfall_Hpa = REAL( Shortfall(I_Local,J_Local), f8 )
          Native_Created_Hpa = Deficit_Hpa - Immediate_Withdrawn
          Status_Name = Survey_Status_Name( Status(I_Local,J_Local) )
          WRITE( Survey_Unit, &
               '(a,",",a,8(",",i0),",",a,12(",",es26.17e3),",",a)' ) &
               'qck-bottom-survey-v1', TRIM( Survey_Run_Id ), Model_Date, &
               Model_Time, Elapsed_Seconds, Omp_Thread_Count, Species_Index, &
               I_Global, J_Global, K_Bottom, TRIM( Status_Name ), &
               Deficit_Hpa, Available_Hpa, Shortfall_Hpa, Tolerance_Hpa, &
               Immediate_Hpa, Immediate_Withdrawn, Native_Created_Hpa, &
               Deficit_Hpa * Conversion, Available_Hpa * Conversion, &
               Shortfall_Hpa * Conversion, Native_Created_Hpa * Conversion, &
               Area, 'native_immediate_donor_v1'
       ENDDO
    ENDDO
    FLUSH( Survey_Unit )
    !$OMP END CRITICAL (QCK_BOTTOM_SURVEY_FILE)

  END SUBROUTINE Qck_Bottom_Survey_Write

  CHARACTER(LEN=24) FUNCTION Survey_Status_Name( Status )

    INTEGER, INTENT(IN) :: Status

    SELECT CASE ( Status )
       CASE ( QCK_BOTTOM_EXACT )
          Survey_Status_Name = 'exact'
       CASE ( QCK_BOTTOM_ROUNDOFF )
          Survey_Status_Name = 'roundoff_feasible'
       CASE ( QCK_BOTTOM_UNFILLABLE )
          Survey_Status_Name = 'materially_unfillable'
       CASE ( QCK_BOTTOM_INVALID_DONOR )
          Survey_Status_Name = 'invalid_donor'
       CASE DEFAULT
          Survey_Status_Name = 'unknown'
    END SELECT

  END FUNCTION Survey_Status_Name

  SUBROUTINE Survey_Error( Message, Routine )

    CHARACTER(LEN=*), INTENT(IN) :: Message, Routine
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ErrMsg = TRIM( Message )
    ThisLoc = ' -> at ' // TRIM( Routine ) // ' (in qck_bottom_survey_mod.F90)'
    CALL Error_Stop( ErrMsg, ThisLoc )

  END SUBROUTINE Survey_Error

END MODULE Qck_Bottom_Survey_Mod
