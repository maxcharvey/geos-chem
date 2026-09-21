PROGRAM BRC_CLOUDJ_MAP_TEST

  USE BRC_CLOUDJ_MAP_MOD, ONLY : BUILD_BRC_CLOUDJ_MAP

  IMPLICIT NONE

  CHARACTER(LEN=80)  :: Titles(74)
  CHARACTER(LEN=255) :: ErrMsg
  INTEGER            :: AerMap(11,5), J, RC

  Titles = ''

  CALL BUILD_BRC_CLOUDJ_MAP( .TRUE., 'ORGANIC', 63, Titles, &
                             AerMap, RC, ErrMsg )
  IF ( RC /= 0 ) ERROR STOP 1
  DO J = 1, 5
     IF ( ANY( AerMap(6:10,J) /= 35 + J ) ) ERROR STOP 2
     IF ( AerMap(11,J) /= 36 ) ERROR STOP 3
  ENDDO

  Titles(64:68) = (/ 'WB00', 'WB50', 'WB70', 'WB80', 'WB90' /)
  Titles(69:73) = (/ 'PB00', 'PB50', 'PB70', 'PB80', 'PB90' /)
  Titles(74)    = 'DB00'

  CALL BUILD_BRC_CLOUDJ_MAP( .TRUE., 'DEDICATED', 74, Titles, &
                             AerMap, RC, ErrMsg )
  IF ( RC /= 0 ) ERROR STOP 4
  DO J = 1, 5
     IF ( AerMap(6,J)  /= 63 + J ) ERROR STOP 5
     IF ( AerMap(7,J)  /= 63 + J ) ERROR STOP 6
     IF ( AerMap(8,J)  /= 35 + J ) ERROR STOP 7
     IF ( AerMap(9,J)  /= 63 + J ) ERROR STOP 8
     IF ( AerMap(10,J) /= 68 + J ) ERROR STOP 9
     IF ( AerMap(11,J) /= 74 ) ERROR STOP 10
  ENDDO

  CALL BUILD_BRC_CLOUDJ_MAP( .TRUE., 'DEDICATED', 63, Titles, &
                             AerMap, RC, ErrMsg )
  IF ( RC == 0 ) ERROR STOP 11

  Titles(69) = 'BAD0'
  CALL BUILD_BRC_CLOUDJ_MAP( .TRUE., 'DEDICATED', 74, Titles, &
                             AerMap, RC, ErrMsg )
  IF ( RC == 0 ) ERROR STOP 12

  CALL BUILD_BRC_CLOUDJ_MAP( .FALSE., 'DEDICATED', 63, Titles, &
                             AerMap, RC, ErrMsg )
  IF ( RC /= 0 ) ERROR STOP 13

  WRITE(*,'(a)') 'PASS: executable BrC Cloud-J map regression'

END PROGRAM BRC_CLOUDJ_MAP_TEST
