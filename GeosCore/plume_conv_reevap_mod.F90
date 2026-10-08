!------------------------------------------------------------------------------
! Runtime-gated measurement of signed RAS surface below-cloud re-evaporation
! omitted by the native F(K,NA)>0 DIAG38 gate. Diagnostic only: this module
! never modifies tracer state or a native diagnostic.
!------------------------------------------------------------------------------
MODULE Plume_Conv_Reevap_Mod

  USE ERROR_MOD,       ONLY : Error_Stop
  USE, INTRINSIC       :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE
  USE OMP_LIB,         ONLY : OMP_GET_MAX_THREADS
  USE Precision_Mod,   ONLY : fp, f8
  USE State_Chm_Mod,   ONLY : ChmState, Ind_
  USE TIME_MOD,        ONLY : GET_ELAPSED_SEC, GET_NHMS, GET_NYMD

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: MAX_TAGS = 10
  CHARACTER(LEN=16), PARAMETER :: TAG_NAMES(MAX_TAGS) = (/ &
       'PLUME_SFC       ', 'PLUME_PBL       ', 'PLUME_6535      ', &
       'PLUME_LEV       ', 'PLUME_PROFILE   ', 'PLUME_SFC_PI    ', &
       'PLUME_PBL_PI    ', 'PLUME_6535_PI   ', 'PLUME_LEV_PI    ', &
       'PLUME_PROFILE_PI' /)

  LOGICAL, SAVE :: Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Enabled          = .FALSE.
  LOGICAL, SAVE :: File_Open        = .FALSE.
  LOGICAL, SAVE :: Call_Active      = .FALSE.
  INTEGER, SAVE :: Ledger_Unit      = -1
  INTEGER, SAVE :: Call_Index       = 0
  INTEGER, SAVE :: Model_Date       = 0
  INTEGER, SAVE :: Model_Time       = 0
  INTEGER, SAVE :: Elapsed_Seconds  = 0
  INTEGER, SAVE :: Omp_Thread_Count = 0
  INTEGER, SAVE :: NX_Saved         = 0
  INTEGER, SAVE :: NY_Saved         = 0
  INTEGER, SAVE :: N_Active_Tags    = 5
  INTEGER, SAVE :: Tag_Ids(MAX_TAGS)  = 0
  INTEGER, SAVE :: Advect_Ids(MAX_TAGS) = 0
  INTEGER, SAVE :: Wetdep_Ids(MAX_TAGS) = 0
  CHARACTER(LEN=128), SAVE :: Run_Id      = ''
  CHARACTER(LEN=128), SAVE :: Manifest_Id = ''

  LOGICAL,  ALLOCATABLE, SAVE :: Seen(:,:,:)
  INTEGER,  ALLOCATABLE, SAVE :: Event_Count(:,:,:)
  REAL(fp), ALLOCATABLE, SAVE :: F_Scav(:,:,:), Area(:,:,:)
  REAL(f8), ALLOCATABLE, SAVE :: Gained_Sum(:,:,:), Gross_Wash_Sum(:,:,:)
  REAL(f8), ALLOCATABLE, SAVE :: Signed_Wetloss_Sum(:,:,:)
  REAL(f8), ALLOCATABLE, SAVE :: Signed_Mass_Sum(:,:,:)
  REAL(f8), ALLOCATABLE, SAVE :: Realized_Mass_Sum(:,:,:)

  PUBLIC :: Conv_Reevap_Begin
  PUBLIC :: Conv_Reevap_Record_Surface_Omission
  PUBLIC :: Conv_Reevap_Write

CONTAINS

  SUBROUTINE Initialize_Gate()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( Gate_Initialized ) RETURN
    Gate_Initialized = .TRUE.

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_RAS_SURFACE_REEVAP_LEDGER', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 ) RETURN
    SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
       CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
          Enabled = .TRUE.
       CASE DEFAULT
          RETURN
    END SELECT

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_RAS_SURFACE_REEVAP_INCLUDE_AGED', &
                                   Value, STATUS=Env_Status )
    IF ( Env_Status == 0 ) THEN
       SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
          CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
             N_Active_Tags = MAX_TAGS
       END SELECT
    ENDIF

    ThisLoc = ' -> at Initialize_Gate (in plume_conv_reevap_mod.F90)'
    CALL Read_Required_Environment( 'GC_RAS_SURFACE_REEVAP_LEDGER_RUN_ID', &
                                    Run_Id, ThisLoc )
    CALL Read_Required_Environment( &
         'GC_RAS_SURFACE_REEVAP_LEDGER_MANIFEST_ID', Manifest_Id, ThisLoc )
    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_RAS_SURFACE_REEVAP_LEDGER_FILE', &
                                   Value, STATUS=Env_Status )
    IF ( Env_Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       CALL Error_Stop( 'GC_RAS_SURFACE_REEVAP_LEDGER_FILE is required', &
                        ThisLoc )
    ENDIF
    OPEN( NEWUNIT=Ledger_Unit, FILE=TRIM( Value ), STATUS='NEW', &
          ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
    IF ( IO_Status /= 0 ) THEN
       ErrMsg = 'Could not create RAS surface re-evaporation ledger: ' // &
                TRIM( Value )
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    File_Open = .TRUE.
    WRITE( Ledger_Unit, '(a)' ) &
         'schema_version,run_id,manifest_id,model_date,model_time,' // &
         'elapsed_seconds,convection_call_index,omp_thread_count,tag,' // &
         'model_species_id,advect_id,wetdep_id,i_gc,j_gc,k_gc,event_count,' // &
         'f_scavenging,area_m2,gained_kg_m2,gross_wash_kg_m2,' // &
         'signed_wetloss_kg_m2,signed_wetloss_kg,' // &
         'realized_signed_loss_kg,state_gain_kg'
    FLUSH( Ledger_Unit )
    WRITE( 6, '(a)' ) 'RAS SURFACE REEVAP LEDGER ENABLED: ' // TRIM( Value )

  END SUBROUTINE Initialize_Gate

  SUBROUTINE Read_Required_Environment( Name, Result, ThisLoc )

    CHARACTER(LEN=*), INTENT(IN)  :: Name, ThisLoc
    CHARACTER(LEN=*), INTENT(OUT) :: Result
    CHARACTER(LEN=1024) :: Value
    INTEGER :: Status

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( Name, Value, STATUS=Status )
    IF ( Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
       CALL Error_Stop( TRIM( Name ) // ' is required when enabled', ThisLoc )
    ENDIF
    Result = TRIM( Value )

  END SUBROUTINE Read_Required_Environment

  SUBROUTINE Conv_Reevap_Begin( State_Chm, NX, NY, Do_Wetloss, Is_Ras )

    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER,        INTENT(IN) :: NX, NY
    LOGICAL,        INTENT(IN) :: Do_Wetloss, Is_Ras

    INTEGER :: A, T

    CALL Initialize_Gate()
    IF ( .not. Enabled ) RETURN
    IF ( Call_Active ) THEN
       CALL Ledger_Error( 'Begin called before prior write', &
                          'Conv_Reevap_Begin' )
    ENDIF
    IF ( .not. Is_Ras ) THEN
       CALL Ledger_Error( 'Surface re-evaporation ledger supports RAS only', &
                          'Conv_Reevap_Begin' )
    ENDIF
    IF ( .not. Do_Wetloss ) THEN
       CALL Ledger_Error( 'WetLossConv must be active for this ledger', &
                          'Conv_Reevap_Begin' )
    ENDIF
    IF ( NX <= 0 .or. NY <= 0 ) THEN
       CALL Ledger_Error( 'Invalid grid dimensions', 'Conv_Reevap_Begin' )
    ENDIF

    DO T = 1, N_Active_Tags
       Tag_Ids(T) = Ind_( TRIM( TAG_NAMES(T) ), 'S' )
       IF ( Tag_Ids(T) < 1 .or. Tag_Ids(T) > SIZE( State_Chm%Species ) ) THEN
          CALL Ledger_Error( 'Missing plume tag: ' // TRIM( TAG_NAMES(T) ), &
                             'Conv_Reevap_Begin' )
       ENDIF
       Advect_Ids(T) = 0
       DO A = 1, State_Chm%nAdvect
          IF ( State_Chm%Map_Advect(A) == Tag_Ids(T) ) Advect_Ids(T) = A
       ENDDO
       IF ( Advect_Ids(T) == 0 ) THEN
          CALL Ledger_Error( 'Plume tag is not advected: ' // &
                             TRIM( TAG_NAMES(T) ), 'Conv_Reevap_Begin' )
       ENDIF
       Wetdep_Ids(T) = State_Chm%SpcData(Tag_Ids(T))%Info%WetDepId
       IF ( Wetdep_Ids(T) <= 0 ) THEN
          CALL Ledger_Error( 'Plume tag is not wet-deposited: ' // &
                             TRIM( TAG_NAMES(T) ), 'Conv_Reevap_Begin' )
       ENDIF
    ENDDO

    CALL Ensure_Storage( NX, NY )
    Seen               = .FALSE.
    Event_Count        = 0
    F_Scav             = 0.0_fp
    Area               = 0.0_fp
    Gained_Sum         = 0.0_f8
    Gross_Wash_Sum     = 0.0_f8
    Signed_Wetloss_Sum = 0.0_f8
    Signed_Mass_Sum    = 0.0_f8
    Realized_Mass_Sum  = 0.0_f8

    Call_Index       = Call_Index + 1
    Model_Date       = GET_NYMD()
    Model_Time       = GET_NHMS()
    Elapsed_Seconds  = GET_ELAPSED_SEC()
    Omp_Thread_Count = OMP_GET_MAX_THREADS()
    Call_Active      = .TRUE.

  END SUBROUTINE Conv_Reevap_Begin

  SUBROUTINE Ensure_Storage( NX, NY )

    INTEGER, INTENT(IN) :: NX, NY

    IF ( ALLOCATED( Seen ) ) THEN
       IF ( NX /= NX_Saved .or. NY /= NY_Saved ) THEN
          CALL Ledger_Error( 'Grid shape changed between convection calls', &
                             'Ensure_Storage' )
       ENDIF
       RETURN
    ENDIF
    NX_Saved = NX
    NY_Saved = NY
    ALLOCATE( Seen(NX,NY,N_Active_Tags), &
              Event_Count(NX,NY,N_Active_Tags) )
    ALLOCATE( F_Scav(NX,NY,N_Active_Tags), &
              Area(NX,NY,N_Active_Tags) )
    ALLOCATE( Gained_Sum(NX,NY,N_Active_Tags), &
              Gross_Wash_Sum(NX,NY,N_Active_Tags) )
    ALLOCATE( Signed_Wetloss_Sum(NX,NY,N_Active_Tags) )
    ALLOCATE( Signed_Mass_Sum(NX,NY,N_Active_Tags) )
    ALLOCATE( Realized_Mass_Sum(NX,NY,N_Active_Tags) )

  END SUBROUTINE Ensure_Storage

  SUBROUTINE Conv_Reevap_Record_Surface_Omission( &
       Species_Id, Wetdep_Id, I, J, K, Scavenging_Fraction, Area_M2, &
       Bmass, Q_Before, Q_After, Gained, Washfrac, Wetloss )

    INTEGER,  INTENT(IN) :: Species_Id, Wetdep_Id, I, J, K
    REAL(fp), INTENT(IN) :: Scavenging_Fraction, Area_M2, Bmass
    REAL(fp), INTENT(IN) :: Q_Before, Q_After, Gained, Washfrac, Wetloss

    INTEGER  :: T
    REAL(f8) :: Area8, Realized

    IF ( .not. Enabled ) RETURN
    IF ( .not. Call_Active ) THEN
       CALL Ledger_Error( 'Record called outside active call', &
                          'Conv_Reevap_Record_Surface_Omission' )
    ENDIF
    T = Tag_Index( Species_Id )
    IF ( T == 0 ) RETURN
    IF ( K /= 1 ) THEN
       CALL Ledger_Error( 'Recorder received non-surface plume event', &
                          'Conv_Reevap_Record_Surface_Omission' )
    ENDIF
    IF ( Wetdep_Id /= Wetdep_Ids(T) .or. Scavenging_Fraction > 0.0_fp ) THEN
       CALL Ledger_Error( 'Invalid omitted-diagnostic predicate', &
                          'Conv_Reevap_Record_Surface_Omission' )
    ENDIF
    IF ( .not. All_Finite( (/ Scavenging_Fraction, Area_M2, Bmass, &
                              Q_Before, Q_After, Gained, Washfrac, Wetloss /) ) ) THEN
       CALL Ledger_Error( 'Non-finite surface re-evaporation value', &
                          'Conv_Reevap_Record_Surface_Omission' )
    ENDIF
    IF ( Wetloss == 0.0_fp ) RETURN

    IF ( Seen(I,J,T) ) THEN
       IF ( F_Scav(I,J,T) /= Scavenging_Fraction .or. &
            Area(I,J,T) /= Area_M2 ) THEN
          CALL Ledger_Error( 'Non-invariant cell metadata', &
                             'Conv_Reevap_Record_Surface_Omission' )
       ENDIF
    ELSE
       Seen(I,J,T)   = .TRUE.
       F_Scav(I,J,T) = Scavenging_Fraction
       Area(I,J,T)   = Area_M2
    ENDIF

    Area8 = REAL( Area_M2, f8 )
    Realized = REAL( Q_Before - Q_After, f8 ) * &
               REAL( Bmass, f8 ) * Area8
    Event_Count(I,J,T)        = Event_Count(I,J,T) + 1
    Gained_Sum(I,J,T)         = Gained_Sum(I,J,T) + REAL( Gained, f8 )
    Gross_Wash_Sum(I,J,T)     = Gross_Wash_Sum(I,J,T) + &
                                REAL( Wetloss + Gained, f8 )
    Signed_Wetloss_Sum(I,J,T) = Signed_Wetloss_Sum(I,J,T) + &
                                REAL( Wetloss, f8 )
    Signed_Mass_Sum(I,J,T)    = Signed_Mass_Sum(I,J,T) + &
                                REAL( Wetloss, f8 ) * Area8
    Realized_Mass_Sum(I,J,T)  = Realized_Mass_Sum(I,J,T) + Realized

  END SUBROUTINE Conv_Reevap_Record_Surface_Omission

  SUBROUTINE Conv_Reevap_Write()

    INTEGER :: I, J, T

    IF ( .not. Enabled ) RETURN
    IF ( .not. Call_Active ) THEN
       CALL Ledger_Error( 'Write called without active begin', &
                          'Conv_Reevap_Write' )
    ENDIF
    DO T = 1, N_Active_Tags
    DO J = 1, NY_Saved
    DO I = 1, NX_Saved
       IF ( Seen(I,J,T) ) CALL Write_Record( I, J, T )
    ENDDO
    ENDDO
    ENDDO
    FLUSH( Ledger_Unit )
    Call_Active = .FALSE.

  END SUBROUTINE Conv_Reevap_Write

  SUBROUTINE Write_Record( I, J, T )

    INTEGER, INTENT(IN) :: I, J, T

    WRITE( Ledger_Unit, &
           '(a,",",a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,",",' // &
           'a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",' // &
           'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",' // &
           'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3)' ) &
         'ras-surface-reevap-ledger-v1', TRIM( Run_Id ), &
         TRIM( Manifest_Id ), Model_Date, Model_Time, Elapsed_Seconds, &
         Call_Index, Omp_Thread_Count, TRIM( TAG_NAMES(T) ), Tag_Ids(T), &
         Advect_Ids(T), Wetdep_Ids(T), I, J, 1, Event_Count(I,J,T), &
         F_Scav(I,J,T), Area(I,J,T), Gained_Sum(I,J,T), &
         Gross_Wash_Sum(I,J,T), Signed_Wetloss_Sum(I,J,T), &
         Signed_Mass_Sum(I,J,T), Realized_Mass_Sum(I,J,T), &
         -Realized_Mass_Sum(I,J,T)

  END SUBROUTINE Write_Record

  INTEGER FUNCTION Tag_Index( Species_Id )

    INTEGER, INTENT(IN) :: Species_Id
    INTEGER :: T

    Tag_Index = 0
    DO T = 1, N_Active_Tags
       IF ( Tag_Ids(T) == Species_Id ) THEN
          Tag_Index = T
          RETURN
       ENDIF
    ENDDO

  END FUNCTION Tag_Index

  LOGICAL FUNCTION All_Finite( Values )

    REAL(fp), INTENT(IN) :: Values(:)
    All_Finite = ALL( IEEE_IS_FINITE( Values ) )

  END FUNCTION All_Finite

  SUBROUTINE Ledger_Error( Message, Routine )

    CHARACTER(LEN=*), INTENT(IN) :: Message, Routine
    CHARACTER(LEN=255) :: ThisLoc

    ThisLoc = ' -> at ' // TRIM( Routine ) // &
              ' (in plume_conv_reevap_mod.F90)'
    CALL Error_Stop( Message, ThisLoc )

  END SUBROUTINE Ledger_Error

END MODULE Plume_Conv_Reevap_Mod
