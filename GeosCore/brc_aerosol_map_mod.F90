!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: brc_aerosol_map_mod.F90
!
! !DESCRIPTION: Validates the mapping from hygroscopic species to canonical
!  brown-carbon aerosol optical bins, and preserves the dry DBRC slot when
!  Fast-JX selects an optical-humidity record.
!
! !REVISION HISTORY:
!  23 Jul 2026 - M. Harvey - Initial version
!
! !INTERFACE:
!
MODULE BRC_AEROSOL_MAP_MOD

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: BRC_DRY_AEROSOL_BIN = 11
  INTEGER, PARAMETER, PUBLIC :: BRC_DRY_RH_BIN      = 1

  PUBLIC :: VALIDATE_BRC_AEROSOL_MAP
  PUBLIC :: BRC_FASTJX_RH_BIN

CONTAINS

  SUBROUTINE VALIDATE_BRC_AEROSOL_MAP( Bins, NumBins, DuplicateEntry, MissingBin )

    INTEGER, INTENT(IN)  :: Bins(:)
    INTEGER, INTENT(IN)  :: NumBins
    INTEGER, INTENT(OUT) :: DuplicateEntry
    INTEGER, INTENT(OUT) :: MissingBin

    INTEGER :: N
    LOGICAL :: Seen(NumBins)

    DuplicateEntry = 0
    MissingBin     = 0
    Seen(:)        = .FALSE.

    DO N = 1, SIZE(Bins)
       IF ( Bins(N) < 1 .OR. Bins(N) > NumBins ) THEN
          DuplicateEntry = N
          RETURN
       ENDIF
       IF ( Seen(Bins(N)) ) THEN
          DuplicateEntry = N
          RETURN
       ENDIF
       Seen(Bins(N)) = .TRUE.
    ENDDO

    DO N = 1, NumBins
       IF ( .NOT. Seen(N) ) THEN
          MissingBin = N
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE VALIDATE_BRC_AEROSOL_MAP


  PURE INTEGER FUNCTION BRC_FASTJX_RH_BIN( AerosolBin, AmbientRHBin ) &
       RESULT( OpticalRHBin )

    ! DBRC uses the dry optical carrier in aerosol_mod.  Fast-JX stores each
    ! aerosol/RH pair in a separate online-LUT slot, so route DBRC only to
    ! IR=1 rather than applying ambient-RH optical properties to dry AOD.
    INTEGER, INTENT(IN) :: AerosolBin
    INTEGER, INTENT(IN) :: AmbientRHBin

    OpticalRHBin = AmbientRHBin
    IF ( AerosolBin == BRC_DRY_AEROSOL_BIN ) THEN
       OpticalRHBin = BRC_DRY_RH_BIN
    ENDIF

  END FUNCTION BRC_FASTJX_RH_BIN

END MODULE BRC_AEROSOL_MAP_MOD
!EOP
!------------------------------------------------------------------------------
