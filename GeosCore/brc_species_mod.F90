! Required chemistry species for the optional brown-carbon pathway.
! Kept independent of model state so reduced-species errors are testable.
MODULE BRC_Species_Mod
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_REQUIRED_SPECIES_ERROR
CONTAINS
  PURE FUNCTION BRC_REQUIRED_SPECIES_ERROR( Ids ) RESULT( Message )
    INTEGER, INTENT(IN) :: Ids(3)  ! FSOAS, BRCSOA, WTC, in that order
    CHARACTER(LEN=255) :: Message
    CHARACTER(LEN=6), PARAMETER :: Names(3) = &
         (/ 'FSOAS ', 'BRCSOA', 'WTC   ' /)
    INTEGER :: N

    Message = ''
    IF ( ALL( Ids > 0 ) ) RETURN
    Message = 'brown_carbon: true requires'
    DO N = 1, SIZE( Ids )
       IF ( Ids(N) <= 0 ) Message = TRIM(Message)//' '//TRIM(Names(N))
    ENDDO
  END FUNCTION BRC_REQUIRED_SPECIES_ERROR
END MODULE BRC_Species_Mod
