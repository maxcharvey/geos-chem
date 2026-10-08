!------------------------------------------------------------------------------
! Opt-in, output-only diagnostics for sizing the bounded QCK_BOTTOM policy.
!
! This module never receives or mutates tracer state.  It records scalar
! candidates produced by Qckxyz and emits one deterministic end-of-run summary.
!------------------------------------------------------------------------------
MODULE Qck_Bottom_Sizing_Mod

  USE Precision_Mod, ONLY : fp, f8

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: SPECIES_NAME_LEN = 64

  LOGICAL, SAVE :: Initialized = .FALSE.
  LOGICAL, SAVE :: Enabled = .FALSE.
  INTEGER, SAVE :: Tpcore_Invocation_Count = 0
  REAL(f8), SAVE :: Diagnostic_Cumulative_Closure_Kg = 0.0_f8

  REAL(f8), SAVE :: Max_Attempted_Event_Kg = 0.0_f8
  INTEGER,  SAVE :: Max_Attempted_Date = -1
  INTEGER,  SAVE :: Max_Attempted_Time = -1
  INTEGER,  SAVE :: Max_Attempted_Ordinal = -1
  INTEGER,  SAVE :: Max_Attempted_Species_Index = -1
  INTEGER,  SAVE :: Max_Attempted_I = -1
  INTEGER,  SAVE :: Max_Attempted_J = -1
  INTEGER,  SAVE :: Max_Attempted_K = -1
  INTEGER,  SAVE :: Max_Attempted_Status = -1
  CHARACTER(LEN=SPECIES_NAME_LEN), SAVE :: Max_Attempted_Species_Name = 'NONE'

  REAL(f8), SAVE :: Max_Accepted_Event_Kg = 0.0_f8
  INTEGER,  SAVE :: Max_Accepted_Date = -1
  INTEGER,  SAVE :: Max_Accepted_Time = -1
  INTEGER,  SAVE :: Max_Accepted_Ordinal = -1
  INTEGER,  SAVE :: Max_Accepted_Species_Index = -1
  INTEGER,  SAVE :: Max_Accepted_I = -1
  INTEGER,  SAVE :: Max_Accepted_J = -1
  INTEGER,  SAVE :: Max_Accepted_K = -1
  CHARACTER(LEN=SPECIES_NAME_LEN), SAVE :: Max_Accepted_Species_Name = 'NONE'

  REAL(f8), SAVE :: Max_Species_Call_Kg = 0.0_f8
  INTEGER,  SAVE :: Max_Call_Date = -1
  INTEGER,  SAVE :: Max_Call_Time = -1
  INTEGER,  SAVE :: Max_Call_Ordinal = -1
  INTEGER,  SAVE :: Max_Call_Species_Index = -1
  CHARACTER(LEN=SPECIES_NAME_LEN), SAVE :: Max_Call_Species_Name = 'NONE'

  PUBLIC :: Qck_Bottom_Sizing_Candidate_Is_Better
  PUBLIC :: Qck_Bottom_Sizing_Enabled
  PUBLIC :: Qck_Bottom_Sizing_Next_Tpcore_Ordinal
  PUBLIC :: Qck_Bottom_Sizing_Record_Call
  PUBLIC :: Qck_Bottom_Sizing_Reset
  PUBLIC :: Qck_Bottom_Sizing_Write_Summary

CONTAINS

  SUBROUTINE Qck_Bottom_Sizing_Reset

    Initialized = .FALSE.
    Enabled = .FALSE.
    Tpcore_Invocation_Count = 0
    Diagnostic_Cumulative_Closure_Kg = 0.0_f8

    Max_Attempted_Event_Kg = 0.0_f8
    Max_Attempted_Date = -1
    Max_Attempted_Time = -1
    Max_Attempted_Ordinal = -1
    Max_Attempted_Species_Index = -1
    Max_Attempted_I = -1
    Max_Attempted_J = -1
    Max_Attempted_K = -1
    Max_Attempted_Status = -1
    Max_Attempted_Species_Name = 'NONE'

    Max_Accepted_Event_Kg = 0.0_f8
    Max_Accepted_Date = -1
    Max_Accepted_Time = -1
    Max_Accepted_Ordinal = -1
    Max_Accepted_Species_Index = -1
    Max_Accepted_I = -1
    Max_Accepted_J = -1
    Max_Accepted_K = -1
    Max_Accepted_Species_Name = 'NONE'

    Max_Species_Call_Kg = 0.0_f8
    Max_Call_Date = -1
    Max_Call_Time = -1
    Max_Call_Ordinal = -1
    Max_Call_Species_Index = -1
    Max_Call_Species_Name = 'NONE'

  END SUBROUTINE Qck_Bottom_Sizing_Reset

  LOGICAL FUNCTION Qck_Bottom_Sizing_Enabled() RESULT( Is_Enabled )

    CHARACTER(LEN=1024) :: Value
    INTEGER             :: Status

    IF ( .not. Initialized ) THEN
!$OMP CRITICAL( Qck_Bottom_Sizing_Initialization )
       IF ( .not. Initialized ) THEN
          Value = ''
          CALL GET_ENVIRONMENT_VARIABLE( &
               'GC_QCK_BOTTOM_SIZING_DIAGNOSTICS', Value, STATUS=Status )
          IF ( Status /= 0 .or. LEN_TRIM( Value ) == 0 ) THEN
             Enabled = .FALSE.
          ELSEIF ( TRIM( ADJUSTL( Value ) ) == '1' ) THEN
             Enabled = .TRUE.
          ELSE
             WRITE( 6, '(a)' ) &
                  'ERROR: GC_QCK_BOTTOM_SIZING_DIAGNOSTICS must be unset or 1'
             ERROR STOP 1
          ENDIF
          Initialized = .TRUE.
       ENDIF
!$OMP END CRITICAL( Qck_Bottom_Sizing_Initialization )
    ENDIF

    Is_Enabled = Enabled

  END FUNCTION Qck_Bottom_Sizing_Enabled

  INTEGER FUNCTION Qck_Bottom_Sizing_Next_Tpcore_Ordinal() RESULT( Ordinal )

    IF ( .not. Enabled ) THEN
       Ordinal = -1
       RETURN
    ENDIF

    Tpcore_Invocation_Count = Tpcore_Invocation_Count + 1
    Ordinal = Tpcore_Invocation_Count

  END FUNCTION Qck_Bottom_Sizing_Next_Tpcore_Ordinal

  PURE LOGICAL FUNCTION Qck_Bottom_Sizing_Candidate_Is_Better( &
       Candidate_Value, Candidate_Date, Candidate_Time, Candidate_Ordinal, &
       Candidate_Species, Candidate_K, Candidate_J, Candidate_I, &
       Current_Value, Current_Date, Current_Time, Current_Ordinal, &
       Current_Species, Current_K, Current_J, Current_I ) RESULT( Is_Better )

    REAL(f8), INTENT(IN) :: Candidate_Value, Current_Value
    INTEGER,  INTENT(IN) :: Candidate_Date, Candidate_Time, Candidate_Ordinal
    INTEGER,  INTENT(IN) :: Candidate_Species, Candidate_K, Candidate_J
    INTEGER,  INTENT(IN) :: Candidate_I
    INTEGER,  INTENT(IN) :: Current_Date, Current_Time, Current_Ordinal
    INTEGER,  INTENT(IN) :: Current_Species, Current_K, Current_J, Current_I

    Is_Better = .FALSE.
    IF ( Candidate_Value > Current_Value ) THEN
       Is_Better = .TRUE.
       RETURN
    ENDIF
    IF ( Candidate_Value < Current_Value .or. Candidate_Value <= 0.0_f8 ) RETURN

    IF ( Candidate_Date /= Current_Date ) THEN
       Is_Better = Candidate_Date < Current_Date
    ELSEIF ( Candidate_Time /= Current_Time ) THEN
       Is_Better = Candidate_Time < Current_Time
    ELSEIF ( Candidate_Ordinal /= Current_Ordinal ) THEN
       Is_Better = Candidate_Ordinal < Current_Ordinal
    ELSEIF ( Candidate_Species /= Current_Species ) THEN
       Is_Better = Candidate_Species < Current_Species
    ELSEIF ( Candidate_K /= Current_K ) THEN
       Is_Better = Candidate_K < Current_K
    ELSEIF ( Candidate_J /= Current_J ) THEN
       Is_Better = Candidate_J < Current_J
    ELSE
       Is_Better = Candidate_I < Current_I
    ENDIF

  END FUNCTION Qck_Bottom_Sizing_Candidate_Is_Better

  SUBROUTINE Qck_Bottom_Sizing_Record_Call( &
       Attempted_Event_Kg, Attempted_I, Attempted_J, Attempted_K, &
       Attempted_Status, Accepted_Event_Kg, Accepted_I, Accepted_J, &
       Accepted_K, Species_Call_Kg, Species_Index, Species_Name, Model_Date, &
       Model_Time, Tpcore_Ordinal )

    REAL(f8), INTENT(IN) :: Attempted_Event_Kg, Accepted_Event_Kg
    REAL(f8), INTENT(IN) :: Species_Call_Kg
    INTEGER,  INTENT(IN) :: Attempted_I, Attempted_J, Attempted_K
    INTEGER,  INTENT(IN) :: Attempted_Status
    INTEGER,  INTENT(IN) :: Accepted_I, Accepted_J, Accepted_K
    INTEGER,  INTENT(IN) :: Species_Index, Model_Date, Model_Time
    INTEGER,  INTENT(IN) :: Tpcore_Ordinal
    CHARACTER(LEN=*), INTENT(IN) :: Species_Name

    IF ( .not. Enabled ) RETURN

!$OMP CRITICAL( Qck_Bottom_Sizing_Record )
    Diagnostic_Cumulative_Closure_Kg = &
         Diagnostic_Cumulative_Closure_Kg + Species_Call_Kg

    IF ( Qck_Bottom_Sizing_Candidate_Is_Better( &
         Attempted_Event_Kg, Model_Date, Model_Time, Tpcore_Ordinal, &
         Species_Index, Attempted_K, Attempted_J, Attempted_I, &
         Max_Attempted_Event_Kg, Max_Attempted_Date, Max_Attempted_Time, &
         Max_Attempted_Ordinal, Max_Attempted_Species_Index, &
         Max_Attempted_K, Max_Attempted_J, Max_Attempted_I ) ) THEN
       Max_Attempted_Event_Kg = Attempted_Event_Kg
       Max_Attempted_Date = Model_Date
       Max_Attempted_Time = Model_Time
       Max_Attempted_Ordinal = Tpcore_Ordinal
       Max_Attempted_Species_Index = Species_Index
       Max_Attempted_Species_Name = TRIM( ADJUSTL( Species_Name ) )
       Max_Attempted_I = Attempted_I
       Max_Attempted_J = Attempted_J
       Max_Attempted_K = Attempted_K
       Max_Attempted_Status = Attempted_Status
    ENDIF

    IF ( Qck_Bottom_Sizing_Candidate_Is_Better( &
         Accepted_Event_Kg, Model_Date, Model_Time, Tpcore_Ordinal, &
         Species_Index, Accepted_K, Accepted_J, Accepted_I, &
         Max_Accepted_Event_Kg, Max_Accepted_Date, Max_Accepted_Time, &
         Max_Accepted_Ordinal, Max_Accepted_Species_Index, &
         Max_Accepted_K, Max_Accepted_J, Max_Accepted_I ) ) THEN
       Max_Accepted_Event_Kg = Accepted_Event_Kg
       Max_Accepted_Date = Model_Date
       Max_Accepted_Time = Model_Time
       Max_Accepted_Ordinal = Tpcore_Ordinal
       Max_Accepted_Species_Index = Species_Index
       Max_Accepted_Species_Name = TRIM( ADJUSTL( Species_Name ) )
       Max_Accepted_I = Accepted_I
       Max_Accepted_J = Accepted_J
       Max_Accepted_K = Accepted_K
    ENDIF

    IF ( Qck_Bottom_Sizing_Candidate_Is_Better( &
         Species_Call_Kg, Model_Date, Model_Time, Tpcore_Ordinal, &
         Species_Index, -1, -1, -1, Max_Species_Call_Kg, Max_Call_Date, &
         Max_Call_Time, Max_Call_Ordinal, Max_Call_Species_Index, &
         -1, -1, -1 ) ) THEN
       Max_Species_Call_Kg = Species_Call_Kg
       Max_Call_Date = Model_Date
       Max_Call_Time = Model_Time
       Max_Call_Ordinal = Tpcore_Ordinal
       Max_Call_Species_Index = Species_Index
       Max_Call_Species_Name = TRIM( ADJUSTL( Species_Name ) )
    ENDIF
!$OMP END CRITICAL( Qck_Bottom_Sizing_Record )

  END SUBROUTINE Qck_Bottom_Sizing_Record_Call

  SUBROUTINE Qck_Bottom_Sizing_Write_Summary( Cumulative_Closure_Kg, &
                                               Event_Max_Kg, Call_Max_Kg, &
                                               Run_Max_Kg )

    REAL(f8), INTENT(IN) :: Cumulative_Closure_Kg
    REAL(fp), INTENT(IN) :: Event_Max_Kg, Call_Max_Kg, Run_Max_Kg

    IF ( .not. Enabled ) RETURN

    WRITE( 6, '(a)' ) &
         'QCK_BOTTOM_SIZING_SUMMARY schema_version=1 enabled=1'
    WRITE( 6, '(a,es22.14,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)' ) &
         'QCK_BOTTOM_SIZING_MAX_ATTEMPTED event_kg=', &
         Max_Attempted_Event_Kg, ' species_index=', &
         Max_Attempted_Species_Index, ' species_name=', &
         TRIM( Max_Attempted_Species_Name ), ' nymd=', Max_Attempted_Date, &
         ' nhms=', Max_Attempted_Time, ' ordinal=', Max_Attempted_Ordinal, &
         ' i=', Max_Attempted_I, ' j=', Max_Attempted_J, ' k=', Max_Attempted_K, &
         ' status=', Max_Attempted_Status
    WRITE( 6, '(a,es22.14,a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)' ) &
         'QCK_BOTTOM_SIZING_MAX_ACCEPTED event_kg=', &
         Max_Accepted_Event_Kg, ' species_index=', &
         Max_Accepted_Species_Index, ' species_name=', &
         TRIM( Max_Accepted_Species_Name ), ' nymd=', Max_Accepted_Date, &
         ' nhms=', Max_Accepted_Time, ' ordinal=', Max_Accepted_Ordinal, &
         ' i=', Max_Accepted_I, ' j=', Max_Accepted_J, ' k=', Max_Accepted_K
    WRITE( 6, '(a,es22.14,a,i0,a,a,a,i0,a,i0,a,i0)' ) &
         'QCK_BOTTOM_SIZING_MAX_CALL closure_kg=', Max_Species_Call_Kg, &
         ' species_index=', Max_Call_Species_Index, ' species_name=', &
         TRIM( Max_Call_Species_Name ), ' nymd=', Max_Call_Date, &
         ' nhms=', Max_Call_Time, ' ordinal=', Max_Call_Ordinal
    WRITE( 6, '(a,es22.14,a,es22.14,a,es22.14,a,es22.14,a,es22.14)' ) &
         'QCK_BOTTOM_SIZING_TOTAL closure_kg=', &
         Diagnostic_Cumulative_Closure_Kg, ' bounded_closure_kg=', &
         Cumulative_Closure_Kg, ' event_max_kg=', Event_Max_Kg, &
         ' call_max_kg=', Call_Max_Kg, ' run_max_kg=', Run_Max_Kg

  END SUBROUTINE Qck_Bottom_Sizing_Write_Summary

END MODULE Qck_Bottom_Sizing_Mod
