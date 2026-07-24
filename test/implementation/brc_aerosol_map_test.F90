PROGRAM BRC_AEROSOL_MAP_TEST

  USE BRC_AEROSOL_MAP_MOD, ONLY : VALIDATE_BRC_AEROSOL_MAP

  IMPLICIT NONE

  INTEGER :: DuplicateEntry, MissingBin

  CALL Check( (/ 6, 1, 11, 3, 8, 4, 2, 10, 5, 9, 7 /), 0, 0 )
  CALL Check( (/ 1, 2, 3, 3, 5, 6, 7, 8, 9, 10, 11 /), 4, 0 )
  CALL Check( (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 11 /), 0, 10 )

  WRITE(*,'(a)') 'PASS: executable BrC aerosol-map regression'

CONTAINS

  SUBROUTINE Check( Bins, ExpectedDuplicate, ExpectedMissing )

    INTEGER, INTENT(IN) :: Bins(:)
    INTEGER, INTENT(IN) :: ExpectedDuplicate
    INTEGER, INTENT(IN) :: ExpectedMissing

    CALL VALIDATE_BRC_AEROSOL_MAP( Bins, 11, DuplicateEntry, MissingBin )
    IF ( DuplicateEntry /= ExpectedDuplicate ) ERROR STOP 1
    IF ( MissingBin /= ExpectedMissing ) ERROR STOP 2

  END SUBROUTINE Check

END PROGRAM BRC_AEROSOL_MAP_TEST
