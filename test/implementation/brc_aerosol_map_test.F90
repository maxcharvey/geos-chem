PROGRAM BRC_AEROSOL_MAP_TEST

  USE BRC_AEROSOL_MAP_MOD, ONLY : VALIDATE_BRC_AEROSOL_MAP, &
                                  BRC_DRY_AEROSOL_BIN, BRC_DRY_RH_BIN, &
                                  BRC_FASTJX_RH_BIN

  IMPLICIT NONE

  INTEGER :: DuplicateEntry, MissingBin

  CALL Check( (/ 6, 1, 11, 3, 8, 4, 2, 10, 5, 9, 7 /), 0, 0 )
  CALL Check( (/ 1, 2, 3, 3, 5, 6, 7, 8, 9, 10, 11 /), 4, 0 )
  CALL Check( (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 11 /), 0, 10 )
  CALL CheckFastJXDrySlot()

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


  SUBROUTINE CheckFastJXDrySlot()

    INTEGER :: AerosolBin, AmbientRHBin

    DO AmbientRHBin = 1, 5
       IF ( BRC_FASTJX_RH_BIN( BRC_DRY_AEROSOL_BIN, AmbientRHBin ) /= &
            BRC_DRY_RH_BIN ) ERROR STOP 3
    ENDDO

    DO AerosolBin = 1, 11
       IF ( AerosolBin == BRC_DRY_AEROSOL_BIN ) CYCLE
       DO AmbientRHBin = 1, 5
          IF ( BRC_FASTJX_RH_BIN( AerosolBin, AmbientRHBin ) /= &
               AmbientRHBin ) ERROR STOP 4
       ENDDO
    ENDDO

  END SUBROUTINE CheckFastJXDrySlot

END PROGRAM BRC_AEROSOL_MAP_TEST
