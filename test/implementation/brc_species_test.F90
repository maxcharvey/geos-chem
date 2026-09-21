PROGRAM BRC_SPECIES_TEST
  USE BRC_Species_Mod, ONLY : BRC_REQUIRED_SPECIES_ERROR
  IMPLICIT NONE
  INTEGER :: Mask, N, Ids(3)
  CHARACTER(LEN=255) :: Message, Expected
  CHARACTER(LEN=6), PARAMETER :: Names(3) = &
       (/ 'FSOAS ', 'BRCSOA', 'WTC   ' /)

  ! All required-species presence combinations; optional species are not
  ! inputs to this guard and therefore cannot suppress an independent step.
  DO Mask = 0, 7
     Ids = 1
     Expected = 'brown_carbon: true requires'
     DO N = 1, 3
        IF ( .NOT. BTEST(Mask,N-1) ) THEN
           Ids(N) = 0
           Expected = TRIM(Expected)//' '//TRIM(Names(N))
        ENDIF
     ENDDO
     IF ( Mask == 7 ) Expected = ''
     Message = BRC_REQUIRED_SPECIES_ERROR( Ids )
     IF ( Message /= Expected ) ERROR STOP 'Incorrect required-species error'
  ENDDO
  Message = BRC_REQUIRED_SPECIES_ERROR( (/ -1, 42, -99 /) )
  IF ( Message /= 'brown_carbon: true requires FSOAS WTC' ) &
       ERROR STOP 'Negative species indices must also fail'
  WRITE(*,'(a)') 'PASS: BrC required-species regression (9 cases)'
END PROGRAM BRC_SPECIES_TEST
