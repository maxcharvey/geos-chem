! Minimal interfaces for exercising the actual History_Write routine without
! NetCDF or a model run.  These stubs observe writes; they do not model History_Update.
MODULE History_Fixture_Types
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: real64
  IMPLICIT NONE
  INTEGER, PARAMETER :: fp=real64, f8=real64
  REAL(fp), PARAMETER :: EPS=1.0e-6_fp
  TYPE OptInput
    LOGICAL :: unused=.FALSE.
  END TYPE
  TYPE GrdState
    INTEGER :: unused=0
  END TYPE
  TYPE SpeciesState
    REAL(fp), POINTER :: Conc(:,:,:)=>NULL()
    INTEGER :: Units=1
  END TYPE
  TYPE ChmState
    TYPE(SpeciesState), ALLOCATABLE :: Species(:)
  END TYPE
  TYPE DgnMap
    INTEGER :: nSlots=0
    INTEGER, ALLOCATABLE :: slot2id(:)
  END TYPE
  TYPE DgnState
    TYPE(DgnMap), POINTER :: Map_SpeciesBC=>NULL()
    REAL(f8), POINTER :: SpeciesBC(:,:,:,:)=>NULL()
    LOGICAL :: Archive_SatDiagnCount=.FALSE., Archive_SatDiagnEdgeCount=.FALSE.
    REAL(f8), POINTER :: SatDiagnCount=>NULL(), SatDiagnEdgeCount=>NULL()
  END TYPE
  TYPE HistContainer
    CHARACTER(LEN=32) :: Name=''
    CHARACTER(LEN=20) :: Spc_Units=''
    REAL(fp) :: ElapsedSec=0, FileWriteAlarm=3600, FileCloseAlarm=0
    REAL(fp) :: FileWriteIvalSec=3600, FileCloseIvalSec=3600
    INTEGER :: FileCloseYmd=0, FileWriteYmd=0
    INTEGER :: writes=0, defines=0, closes=0
    REAL(fp) :: captured(2)=-1
  END TYPE
  TYPE MetaHistContainer
    TYPE(HistContainer), POINTER :: Container=>NULL()
    TYPE(MetaHistContainer), POINTER :: Next=>NULL()
  END TYPE
  TYPE(MetaHistContainer), POINTER :: CollectionList=>NULL()
  INTEGER :: update_calls=0
END MODULE

MODULE Charpak_Mod
CONTAINS
  FUNCTION To_UpperCase(s) RESULT(out)
    CHARACTER(*), INTENT(IN) :: s
    CHARACTER(LEN(s)) :: out
    INTEGER :: i,n
    out=s
    DO i=1,LEN(s)
      n=IACHAR(s(i:i))
      IF(n>=IACHAR('a').AND.n<=IACHAR('z')) out(i:i)=ACHAR(n-32)
    ENDDO
  END FUNCTION
END MODULE
MODULE ErrCode_Mod
  INTEGER, PARAMETER :: GC_SUCCESS=0
CONTAINS
  SUBROUTINE GC_Error(message,rc,location)
    CHARACTER(*), INTENT(IN) :: message,location
    INTEGER, INTENT(OUT) :: rc
    rc=1
  END SUBROUTINE
END MODULE
MODULE HistContainer_Mod
  USE History_Fixture_Types, ONLY: HistContainer
END MODULE
MODULE HistItem_Mod
  TYPE HistItem
    INTEGER :: unused
  END TYPE
END MODULE
MODULE Input_Opt_Mod
  USE History_Fixture_Types, ONLY: OptInput
END MODULE
MODULE State_Diag_Mod
  USE History_Fixture_Types, ONLY: DgnMap,DgnState
END MODULE
MODULE State_Chm_Mod
  USE History_Fixture_Types, ONLY: ChmState
END MODULE
MODULE State_Grid_Mod
  USE History_Fixture_Types, ONLY: GrdState
END MODULE
MODULE MetaHistContainer_Mod
  USE History_Fixture_Types, ONLY: MetaHistContainer
END MODULE
MODULE MetaHistItem_Mod
  TYPE MetaHistItem
    INTEGER :: unused
  END TYPE
END MODULE
MODULE PhysConstants
  REAL, PARAMETER :: AIRMW=28.9644
END MODULE
MODULE Registry_Params_Mod
END MODULE
MODULE UnitConv_Mod
  USE History_Fixture_Types, ONLY: ChmState
  INTEGER, PARAMETER :: MOLES_SPECIES_PER_MOLES_DRY_AIR=1
  CHARACTER(LEN=20), PARAMETER :: UNIT_STR(1)=['mol mol-1 dry       ']
CONTAINS
  LOGICAL FUNCTION Check_Units(state,units)
    TYPE(ChmState), INTENT(IN) :: state
    INTEGER, INTENT(IN) :: units
    Check_Units=ALL(state%Species%Units==units)
  END FUNCTION
END MODULE
MODULE History_Util_Mod
  USE History_Fixture_Types, ONLY: OptInput,HistContainer
CONTAINS
  SUBROUTINE HistContainer_FileCloseIvalSet(input,container,rc)
    TYPE(OptInput), INTENT(IN) :: input
    TYPE(HistContainer), INTENT(INOUT) :: container
    INTEGER, INTENT(OUT) :: rc
    rc=0
  END SUBROUTINE
  SUBROUTINE HistContainer_FileWriteIvalSet(input,container,rc)
    TYPE(OptInput), INTENT(IN) :: input
    TYPE(HistContainer), INTENT(INOUT) :: container
    INTEGER, INTENT(OUT) :: rc
    rc=0
  END SUBROUTINE
  SUBROUTINE SatDiagn_or_SatDiagnEdge(name,sat,edge)
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL, INTENT(OUT) :: sat,edge
    sat=(TRIM(name)=='SATDIAGN')
    edge=(TRIM(name)=='SATDIAGNEDGE')
  END SUBROUTINE
END MODULE
MODULE History_Netcdf_Mod
  USE History_Fixture_Types, ONLY: OptInput,DgnState,GrdState,HistContainer
CONTAINS
  SUBROUTINE History_Netcdf_Close(container,rc)
    TYPE(HistContainer), INTENT(INOUT) :: container
    INTEGER, INTENT(OUT) :: rc
    container%closes=container%closes+1
    rc=0
  END SUBROUTINE
  SUBROUTINE History_Netcdf_Define(input_opt,state_grid,container,rc)
    TYPE(OptInput), INTENT(IN) :: input_opt
    TYPE(GrdState), INTENT(IN) :: state_grid
    TYPE(HistContainer), INTENT(INOUT) :: container
    INTEGER, INTENT(OUT) :: rc
    container%defines=container%defines+1
    rc=0
  END SUBROUTINE
  SUBROUTINE History_Netcdf_Write(input_opt,state_diag,container,rc)
    TYPE(OptInput), INTENT(IN) :: input_opt
    TYPE(DgnState), INTENT(IN) :: state_diag
    TYPE(HistContainer), INTENT(INOUT) :: container
    INTEGER, INTENT(OUT) :: rc
    container%writes=container%writes+1
    container%captured=state_diag%SpeciesBC(1,1,1,:)
    rc=0
  END SUBROUTINE
END MODULE
