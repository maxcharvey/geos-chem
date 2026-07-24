!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: brc_aerosol_map_mod.F90
!
! !DESCRIPTION: Validates the mapping from hygroscopic species to canonical
!  brown-carbon aerosol optical bins.
!
! !REVISION HISTORY:
!  23 Jul 2026 - M. Harvey - Initial version
!
! !INTERFACE:
!
MODULE BRC_AEROSOL_MAP_MOD

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: VALIDATE_BRC_AEROSOL_MAP

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

END MODULE BRC_AEROSOL_MAP_MOD
!EOP
!------------------------------------------------------------------------------
