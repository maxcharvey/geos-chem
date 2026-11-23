PROGRAM BRC_SPECIES_TEST
  USE BRC_Species_Mod, ONLY : BRC_REQUIRED_SPECIES_ERROR
  IMPLICIT NONE
  INTEGER :: Mask, N, Ids(5)
  CHARACTER(LEN=255) :: Message, Expected
  CHARACTER(LEN=6), PARAMETER :: Names(5) = &
       (/ 'FSOAS ', 'BRCSOA', 'WTC   ', 'PWTC  ', 'SWTC  ' /)

  ! All required-species presence combinations; optional species are not
  ! inputs to this guard and therefore cannot suppress an independent step.
  DO Mask = 0, 31
     Ids = 1
     Expected = 'brown_carbon: true requires'
     DO N = 1, 5
        IF ( .NOT. BTEST(Mask,N-1) ) THEN
           Ids(N) = 0
           Expected = TRIM(Expected)//' '//TRIM(Names(N))
        ENDIF
     ENDDO
     IF ( Mask == 31 ) Expected = ''
     Message = BRC_REQUIRED_SPECIES_ERROR( Ids )
     IF ( Message /= Expected ) ERROR STOP 'Incorrect required-species error'
  ENDDO
  Message = BRC_REQUIRED_SPECIES_ERROR( (/ -1, 42, -99, 4, 5 /) )
  IF ( Message /= 'brown_carbon: true requires FSOAS WTC' ) &
       ERROR STOP 'Negative species indices must also fail'
  WRITE(*,'(a)') 'PASS: BrC required-species regression (33 cases)'
END PROGRAM BRC_SPECIES_TEST
