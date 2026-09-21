!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: brc_cloudj_map_mod.F90
!
! !DESCRIPTION: Defines explicit Cloud-J optical mappings for BrC aerosol
!  bins.  Organic mode supports optical-equivalence tests.  Dedicated mode
!  requires generated wet-BrC, persistent-BrC, and dry-BrC records.
!\\
!\\
! !INTERFACE:
!
MODULE BRC_CLOUDJ_MAP_MOD
!
! !USES:
!
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: BUILD_BRC_CLOUDJ_MAP
!
! !REVISION HISTORY:
!  23 Jul 2026 - M. Harvey - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

  INTEGER, PARAMETER :: N_BRC_AER = 11
  INTEGER, PARAMETER :: N_BRC_RH  = 5

  INTEGER, PARAMETER :: WET_BRC_FIRST  = 64
  INTEGER, PARAMETER :: PBRC_FIRST     = 69
  INTEGER, PARAMETER :: DRY_BRC_RECORD = 74

  CHARACTER(LEN=4), PARAMETER :: WET_BRC_TITLE(N_BRC_RH) =       &
     (/ 'WB00', 'WB50', 'WB70', 'WB80', 'WB90' /)
  CHARACTER(LEN=4), PARAMETER :: PBRC_TITLE(N_BRC_RH) =          &
     (/ 'PB00', 'PB50', 'PB70', 'PB80', 'PB90' /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: build_brc_cloudj_map
!
! !DESCRIPTION: Returns the Cloud-J record assigned to each GEOS-Chem
!  hygroscopic aerosol and RH slot.  DBRCPOA always uses one dry record.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BUILD_BRC_CLOUDJ_MAP( BrC_Active, Mode, NAA, Titles, &
                                   AerMap, RC, ErrMsg )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)  :: BrC_Active
    CHARACTER(LEN=*), INTENT(IN)  :: Mode
    INTEGER,          INTENT(IN)  :: NAA
    CHARACTER(LEN=*), INTENT(IN)  :: Titles(:)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: AerMap(:,:)
    INTEGER,          INTENT(OUT) :: RC
    CHARACTER(LEN=*), INTENT(OUT) :: ErrMsg
!
! !REVISION HISTORY:
!  23 Jul 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=80) :: MieTitle
    INTEGER           :: J, Record

    RC     = 0
    ErrMsg = ''

    IF ( SIZE(AerMap,1) /= N_BRC_AER .OR. &
         SIZE(AerMap,2) /= N_BRC_RH ) THEN
       WRITE( ErrMsg, '(a,i0,a,i0)' ) 'Cloud-J BrC map requires ', &
          N_BRC_AER, ' aerosol bins and ', N_BRC_RH
       RC = -1
       RETURN
    ENDIF

    DO J = 1, N_BRC_RH
       AerMap(:,J) = (/ 21+J, 28+J, 35+J, 42+J, 49+J, &
                        35+J, 35+J, 35+J, 35+J, 35+J, 36 /)
    ENDDO

    IF ( .NOT. BrC_Active ) RETURN

    SELECT CASE ( TRIM( Mode ) )
       CASE ( 'ORGANIC' )
          ! BrC wet bins use OC records 36-40.  Dry DBRC always uses OC00.
          RETURN

       CASE ( 'DEDICATED' )
          IF ( NAA < DRY_BRC_RECORD .OR. SIZE(Titles) < DRY_BRC_RECORD ) THEN
             WRITE( ErrMsg, '(a,i0,a,i0)' )                              &
                'Dedicated Cloud-J BrC optics require records 64-',      &
                DRY_BRC_RECORD, '; table has ', NAA
             RC = -1
             RETURN
          ENDIF

          DO J = 1, N_BRC_RH
             Record   = WET_BRC_FIRST + J - 1
             MieTitle = ADJUSTL( Titles(Record) )
             IF ( MieTitle(1:4) /= WET_BRC_TITLE(J) ) THEN
                WRITE( ErrMsg, '(a,i0,5a)' ) 'Cloud-J record ', Record, &
                   ' has title "', TRIM(Titles(Record)), '"; expected "', &
                   WET_BRC_TITLE(J), '"'
                RC = -1
                RETURN
             ENDIF

             Record   = PBRC_FIRST + J - 1
             MieTitle = ADJUSTL( Titles(Record) )
             IF ( MieTitle(1:4) /= PBRC_TITLE(J) ) THEN
                WRITE( ErrMsg, '(a,i0,5a)' ) 'Cloud-J record ', Record, &
                   ' has title "', TRIM(Titles(Record)), '"; expected "', &
                   PBRC_TITLE(J), '"'
                RC = -1
                RETURN
             ENDIF
          ENDDO

          MieTitle = ADJUSTL( Titles(DRY_BRC_RECORD) )
          IF ( MieTitle(1:4) /= 'DB00' ) THEN
             WRITE( ErrMsg, '(a,i0,5a)' ) 'Cloud-J record ',            &
                DRY_BRC_RECORD, ' has title "',                         &
                TRIM(Titles(DRY_BRC_RECORD)), '"; expected "',          &
                'DB00', '"'
             RC = -1
             RETURN
          ENDIF

          DO J = 1, N_BRC_RH
             AerMap(6,J)  = WET_BRC_FIRST + J - 1
             AerMap(7,J)  = WET_BRC_FIRST + J - 1
             AerMap(8,J)  = 35 + J
             AerMap(9,J)  = WET_BRC_FIRST + J - 1
             AerMap(10,J) = PBRC_FIRST + J - 1
             AerMap(11,J) = DRY_BRC_RECORD
          ENDDO

       CASE DEFAULT
          ErrMsg = 'Cloud-J BrC optics mode must be organic or dedicated'
          RC = -1
    END SELECT

  END SUBROUTINE BUILD_BRC_CLOUDJ_MAP
!EOC

END MODULE BRC_CLOUDJ_MAP_MOD
