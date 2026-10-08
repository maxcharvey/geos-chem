!------------------------------------------------------------------------------
! Runtime-gated fixed-lifetime transfer between matched plume tracer pools.
! The caller invokes this while Do_Chemistry stores species as kg per grid box.
! Disabled mode is a no-op and does not allocate or open a file.
!------------------------------------------------------------------------------
MODULE Plume_Aging_Mod

  USE ERROR_MOD,       ONLY : Error_Stop
  USE ErrCode_Mod,     ONLY : GC_SUCCESS
  USE, INTRINSIC       :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE
  USE OMP_LIB,         ONLY : OMP_GET_MAX_THREADS
  USE Precision_Mod,   ONLY : fp, f8
  USE State_Chm_Mod,   ONLY : ChmState, Ind_
  USE State_Grid_Mod,  ONLY : GrdState
  USE TIME_MOD,        ONLY : GET_ELAPSED_SEC, GET_NHMS, GET_NYMD, GET_TS_CHEM
  USE UnitConv_Mod,    ONLY : KG_SPECIES

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: N_PAIRS = 5
  CHARACTER(LEN=16), PARAMETER :: PARENT_NAMES(N_PAIRS) = (/ &
       'PLUME_SFC       ', 'PLUME_PBL       ', 'PLUME_6535      ', &
       'PLUME_LEV       ', 'PLUME_PROFILE   ' /)
  CHARACTER(LEN=16), PARAMETER :: PRODUCT_NAMES(N_PAIRS) = (/ &
       'PLUME_SFC_PI    ', 'PLUME_PBL_PI    ', 'PLUME_6535_PI   ', &
       'PLUME_LEV_PI    ', 'PLUME_PROFILE_PI' /)

  LOGICAL, SAVE :: Gate_Initialized = .FALSE.
  LOGICAL, SAVE :: Enabled          = .FALSE.
  INTEGER, SAVE :: Ledger_Unit      = -1
  INTEGER, SAVE :: Call_Index       = 0
  REAL(f8), SAVE :: Lifetime_S      = 0.0_f8
  CHARACTER(LEN=128), SAVE :: Run_Id      = ''
  CHARACTER(LEN=128), SAVE :: Manifest_Id = ''

  PUBLIC :: Apply_Plume_Aging

CONTAINS

  SUBROUTINE Initialize_Gate()

    CHARACTER(LEN=1024) :: Value
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: Env_Status, IO_Status

    IF ( Gate_Initialized ) RETURN
    Gate_Initialized = .TRUE.

    Value = ''
    CALL GET_ENVIRONMENT_VARIABLE( 'GC_PLUME_AGING', Value, &
                                   STATUS=Env_Status )
    IF ( Env_Status /= 0 ) RETURN
    SELECT CASE ( TRIM( ADJUSTL( Value ) ) )
       CASE ( '1', 'true', 'TRUE', 'yes', 'YES' )
          Enabled = .TRUE.
       CASE DEFAULT
          RETURN
    END SELECT

    ThisLoc = ' -> at Initialize_Gate (in plume_aging_mod.F90)'
    CALL Read_Required_Environment( 'GC_PLUME_AGING_RUN_ID', Run_Id, ThisLoc )
    CALL Read_Required_Environment( 'GC_PLUME_AGING_MANIFEST_ID', &
                                    Manifest_Id, ThisLoc )
    CALL Read_Required_Environment( 'GC_PLUME_AGING_LIFETIME_S', &
                                    Value, ThisLoc )
    READ( Value, *, IOSTAT=IO_Status ) Lifetime_S
    IF ( IO_Status /= 0 .or. .not. IEEE_IS_FINITE( Lifetime_S ) .or. &
         Lifetime_S <= 0.0_f8 ) THEN
       CALL Error_Stop( 'GC_PLUME_AGING_LIFETIME_S must be finite and positive', &
                        ThisLoc )
    ENDIF

    CALL Read_Required_Environment( 'GC_PLUME_AGING_FILE', Value, ThisLoc )
    OPEN( NEWUNIT=Ledger_Unit, FILE=TRIM( Value ), STATUS='NEW', &
          ACTION='WRITE', FORM='FORMATTED', IOSTAT=IO_Status )
    IF ( IO_Status /= 0 ) THEN
       ErrMsg = 'Could not create plume aging ledger: ' // TRIM( Value )
       CALL Error_Stop( ErrMsg, ThisLoc )
    ENDIF
    WRITE( Ledger_Unit, '(a)' ) &
         'schema_version,run_id,manifest_id,model_date,model_time,' // &
         'elapsed_seconds,chemistry_call_index,omp_thread_count,pair_index,' // &
         'parent,product,parent_species_id,product_species_id,' // &
         'chemistry_timestep_s,lifetime_s,decay_factor,parent_before_kg,' // &
         'parent_after_kg,product_before_kg,product_after_kg,' // &
         'analytic_transfer_kg,realized_parent_loss_kg,' // &
         'realized_product_gain_kg,pair_closure_residual_kg'
    FLUSH( Ledger_Unit )
    WRITE( 6, '(a)' ) 'PLUME AGING ENABLED: ' // TRIM( Value )

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

  SUBROUTINE Apply_Plume_Aging( State_Chm, State_Grid, RC )

    TYPE(ChmState),  INTENT(INOUT) :: State_Chm
    TYPE(GrdState),  INTENT(IN)    :: State_Grid
    INTEGER,         INTENT(INOUT) :: RC

    INTEGER  :: A, I, J, K, P, Parent_Id(N_PAIRS), Product_Id(N_PAIRS)
    INTEGER  :: Model_Date, Model_Time, Elapsed_Seconds, Ts_Chem
    LOGICAL  :: Parent_Advected, Product_Advected
    REAL(fp) :: Parent_Old, Parent_New, Transfer
    REAL(f8) :: Decay, Parent_Before, Parent_After, Product_Before
    REAL(f8) :: Product_After, Analytic_Transfer, Parent_Loss, Product_Gain
    CHARACTER(LEN=255) :: ThisLoc

    RC = GC_SUCCESS
    CALL Initialize_Gate()
    IF ( .not. Enabled ) RETURN
    ThisLoc = ' -> at Apply_Plume_Aging (in plume_aging_mod.F90)'

    Ts_Chem = GET_TS_CHEM()
    IF ( Ts_Chem <= 0 ) THEN
       CALL Error_Stop( 'Chemistry timestep must be positive', ThisLoc )
    ENDIF
    Decay = EXP( -REAL( Ts_Chem, f8 ) / Lifetime_S )
    IF ( .not. IEEE_IS_FINITE( Decay ) .or. Decay <= 0.0_f8 .or. &
         Decay > 1.0_f8 ) THEN
       CALL Error_Stop( 'Invalid plume aging decay factor', ThisLoc )
    ENDIF

    DO P = 1, N_PAIRS
       Parent_Id(P)  = Ind_( TRIM( PARENT_NAMES(P) ), 'S' )
       Product_Id(P) = Ind_( TRIM( PRODUCT_NAMES(P) ), 'S' )
       IF ( Parent_Id(P) < 1 .or. Product_Id(P) < 1 ) THEN
          CALL Error_Stop( 'Missing plume aging pair: ' // &
               TRIM( PARENT_NAMES(P) ) // ' -> ' // &
               TRIM( PRODUCT_NAMES(P) ), ThisLoc )
       ENDIF
       IF ( State_Chm%Species(Parent_Id(P))%Units /= KG_SPECIES .or. &
            State_Chm%Species(Product_Id(P))%Units /= KG_SPECIES ) THEN
          CALL Error_Stop( 'Plume aging requires species state in kg', ThisLoc )
       ENDIF
       IF ( .not. ASSOCIATED( State_Chm%Species(Parent_Id(P))%Conc ) .or. &
            .not. ASSOCIATED( State_Chm%Species(Product_Id(P))%Conc ) ) THEN
          CALL Error_Stop( 'Plume aging concentration is not associated', &
                           ThisLoc )
       ENDIF
       IF ( State_Chm%SpcData(Parent_Id(P))%Info%MW_g /= &
            State_Chm%SpcData(Product_Id(P))%Info%MW_g ) THEN
          CALL Error_Stop( 'Plume aging pair molecular weights differ', ThisLoc )
       ENDIF
       Parent_Advected  = .FALSE.
       Product_Advected = .FALSE.
       DO A = 1, State_Chm%nAdvect
          IF ( State_Chm%Map_Advect(A) == Parent_Id(P) ) &
             Parent_Advected = .TRUE.
          IF ( State_Chm%Map_Advect(A) == Product_Id(P) ) &
             Product_Advected = .TRUE.
       ENDDO
       IF ( .not. Parent_Advected .or. .not. Product_Advected ) THEN
          CALL Error_Stop( 'Plume aging pair is not advected', ThisLoc )
       ENDIF
    ENDDO

    Call_Index      = Call_Index + 1
    Model_Date      = GET_NYMD()
    Model_Time      = GET_NHMS()
    Elapsed_Seconds = GET_ELAPSED_SEC()

    DO P = 1, N_PAIRS
       Parent_Before  = 0.0_f8
       Product_Before = 0.0_f8
       DO K = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Parent_Old = State_Chm%Species(Parent_Id(P))%Conc(I,J,K)
          Parent_New = State_Chm%Species(Product_Id(P))%Conc(I,J,K)
          IF ( .not. IEEE_IS_FINITE( Parent_Old ) .or. &
               .not. IEEE_IS_FINITE( Parent_New ) .or. &
               Parent_Old < 0.0_fp .or. Parent_New < 0.0_fp ) THEN
             CALL Error_Stop( 'Non-finite or negative plume aging state', &
                              ThisLoc )
          ENDIF
          Parent_Before  = Parent_Before  + REAL( Parent_Old, f8 )
          Product_Before = Product_Before + REAL( Parent_New, f8 )
       ENDDO
       ENDDO
       ENDDO

       !$OMP PARALLEL DO DEFAULT( SHARED ) PRIVATE( I, J, K, Parent_Old, &
       !$OMP Parent_New, Transfer ) COLLAPSE( 3 )
       DO K = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Parent_Old = State_Chm%Species(Parent_Id(P))%Conc(I,J,K)
          Parent_New = Parent_Old * REAL( Decay, fp )
          Transfer   = Parent_Old - Parent_New
          State_Chm%Species(Parent_Id(P))%Conc(I,J,K)  = Parent_New
          State_Chm%Species(Product_Id(P))%Conc(I,J,K) = &
               State_Chm%Species(Product_Id(P))%Conc(I,J,K) + Transfer
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       Parent_After = 0.0_f8
       Product_After = 0.0_f8
       DO K = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Parent_After = Parent_After + REAL( &
               State_Chm%Species(Parent_Id(P))%Conc(I,J,K), f8 )
          Product_After = Product_After + REAL( &
               State_Chm%Species(Product_Id(P))%Conc(I,J,K), f8 )
       ENDDO
       ENDDO
       ENDDO

       Analytic_Transfer = Parent_Before * ( 1.0_f8 - Decay )
       Parent_Loss       = Parent_Before - Parent_After
       Product_Gain      = Product_After - Product_Before
       WRITE( Ledger_Unit, &
            '(a,",",a,",",a,",",i0,",",i0,",",i0,",",i0,",",i0,",",' // &
            'i0,",",a,",",a,",",i0,",",i0,",",i0,",",es26.17e3,",",' // &
            'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",' // &
            'es26.17e3,",",es26.17e3,",",es26.17e3,",",es26.17e3,",",' // &
            'es26.17e3)' ) &
            'plume-aging-ledger-v1', TRIM( Run_Id ), TRIM( Manifest_Id ), &
            Model_Date, Model_Time, Elapsed_Seconds, Call_Index, &
            OMP_GET_MAX_THREADS(), P, TRIM( PARENT_NAMES(P) ), &
            TRIM( PRODUCT_NAMES(P) ), Parent_Id(P), Product_Id(P), Ts_Chem, &
            Lifetime_S, Decay, Parent_Before, Parent_After, Product_Before, &
            Product_After, Analytic_Transfer, Parent_Loss, Product_Gain, &
            Parent_Loss - Product_Gain
    ENDDO
    FLUSH( Ledger_Unit )

  END SUBROUTINE Apply_Plume_Aging

END MODULE Plume_Aging_Mod
