! Derived white-carbon output. WTC is metadata only; PWTC/SWTC own all mass.
! Derive at sampling in the components' current units, never from a stale mirror.
MODULE BRC_WTC_Output_Mod
  USE Precision_Mod, ONLY: fp, f4
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  TYPE, PUBLIC :: BRC_WTC_Map
     INTEGER :: Total = 0
     INTEGER :: Primary = 0
     INTEGER :: Secondary = 0
  END TYPE
  PUBLIC :: BRC_WTC_Output_Ids, BRC_Copy_Species_Output
  PUBLIC :: BRC_Copy_Species_Slab, BRC_Species_Output_Cell
  PUBLIC :: BRC_WTC_Missing_Restart_Error, BRC_WTC_Total_Matches
  PUBLIC :: BRC_WTC_Components_Valid
CONTAINS
  FUNCTION BRC_WTC_Output_Ids() RESULT(Map)
    TYPE(BRC_WTC_Map) :: Map
    Map%Total = Ind_('WTC')
    IF ( Map%Total <= 0 ) RETURN
    Map%Primary = Ind_('PWTC')
    Map%Secondary = Ind_('SWTC')
    IF ( Map%Primary <= 0 .OR. Map%Secondary <= 0 ) &
         ERROR STOP 'Derived WTC output requires PWTC and SWTC'
  END FUNCTION

  PURE SUBROUTINE BRC_Copy_Species_Output(State, N, Values, Map)
    TYPE(ChmState), INTENT(IN) :: State
    INTEGER, INTENT(IN) :: N
    REAL(fp), INTENT(OUT) :: Values(:,:,:)
    TYPE(BRC_WTC_Map), INTENT(IN) :: Map
    IF ( Map%Total > 0 .AND. N == Map%Total ) THEN
       Values = State%Species(Map%Primary)%Conc + State%Species(Map%Secondary)%Conc
    ELSE
       Values = State%Species(N)%Conc
    ENDIF
  END SUBROUTINE

  PURE SUBROUTINE BRC_Copy_Species_Slab(State, N, I, Values, Map)
    TYPE(ChmState), INTENT(IN) :: State
    INTEGER, INTENT(IN) :: N, I
    REAL(fp), INTENT(OUT) :: Values(:,:)
    TYPE(BRC_WTC_Map), INTENT(IN) :: Map
    IF ( Map%Total > 0 .AND. N == Map%Total ) THEN
       Values = State%Species(Map%Primary)%Conc(I,:,:) + State%Species(Map%Secondary)%Conc(I,:,:)
    ELSE
       Values = State%Species(N)%Conc(I,:,:)
    ENDIF
  END SUBROUTINE

  PURE FUNCTION BRC_Species_Output_Cell(State, N, I, J, L, Map) RESULT(Value)
    TYPE(ChmState), INTENT(IN) :: State
    INTEGER, INTENT(IN) :: N, I, J, L
    TYPE(BRC_WTC_Map), INTENT(IN) :: Map
    REAL(fp) :: Value
    IF ( Map%Total > 0 .AND. N == Map%Total ) THEN
       Value = State%Species(Map%Primary)%Conc(I,J,L) + State%Species(Map%Secondary)%Conc(I,J,L)
    ELSE
       Value = State%Species(N)%Conc(I,J,L)
    ENDIF
  END FUNCTION

  PURE FUNCTION BRC_WTC_Missing_Restart_Error(Name, Found) RESULT(Message)
    CHARACTER(*), INTENT(IN) :: Name
    LOGICAL, INTENT(IN) :: Found
    CHARACTER(255) :: Message
    Message = ''
    IF ( .NOT. Found .AND. (Name == 'PWTC' .OR. Name == 'SWTC') ) &
         Message = 'Split-WTC restart requires explicit '//TRIM(Name)// &
         '; legacy WTC cannot be partitioned. Prepare an explicit fresh seed or a component restart.'
  END FUNCTION

  PURE LOGICAL FUNCTION BRC_WTC_Components_Valid(Primary, Secondary) RESULT(Valid)
    REAL(fp), INTENT(IN) :: Primary(:,:,:), Secondary(:,:,:)
    Valid = .FALSE.
    IF ( ANY(SHAPE(Primary) /= SHAPE(Secondary)) ) RETURN
    IF ( .NOT. ALL(IEEE_IS_FINITE(Primary)) .OR. .NOT. ALL(IEEE_IS_FINITE(Secondary)) ) RETURN
    IF ( ANY(Primary < 0.0_fp) .OR. ANY(Secondary < 0.0_fp) ) RETURN
    ! Check representability before evaluating the sum, including under traps.
    IF ( ANY(Primary > HUGE(1.0_fp) - Secondary) ) RETURN
    Valid = ALL(IEEE_IS_FINITE(Primary + Secondary))
  END FUNCTION

  PURE LOGICAL FUNCTION BRC_WTC_Total_Matches(Total, Primary, Secondary) RESULT(Valid)
    REAL(fp), INTENT(IN) :: Total(:,:,:), Primary(:,:,:), Secondary(:,:,:)
    REAL(fp), PARAMETER :: Relative = 8.0_fp * REAL(EPSILON(1.0_f4),fp)
    REAL(fp), PARAMETER :: Absolute = Relative * REAL(TINY(1.0_f4),fp)
    Valid = .FALSE.
    IF ( ANY(SHAPE(Total) /= SHAPE(Primary)) .OR. ANY(SHAPE(Total) /= SHAPE(Secondary)) ) RETURN
    IF ( .NOT. ALL(IEEE_IS_FINITE(Total)) ) RETURN
    IF ( ANY(Total < 0.0_fp) ) RETURN
    IF ( .NOT. BRC_WTC_Components_Valid(Primary, Secondary) ) RETURN
    Valid = ALL(ABS(Total - (Primary + Secondary)) <= &
         Absolute + Relative * MAX(ABS(Total), ABS(Primary + Secondary)))
  END FUNCTION
END MODULE
