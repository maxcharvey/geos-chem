!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: brc_mod.F90
!
! !DESCRIPTION: Module BRC\_MOD contains routines for brown carbon (BrC)
!  aerosol chemistry, including the darkening chain:
!
!    FSOAS  ---(rapid darkening)---> BRCSOA ---(photobleaching)---> WTC
!
!  This module is designed to be called from CHEMCARBON in carbon\_mod.F90
!  with minimal modifications to existing code.
!\\
!\\
! !INTERFACE:
!
MODULE BRC_MOD
!
! !USES:
!
  USE Precision_Mod    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: ChemBrC
  PUBLIC :: Init_BrC
  PUBLIC :: Cleanup_BrC
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHEM_FSOAS
  PRIVATE :: CHEM_BRCSOA
  PRIVATE :: CHEM_WTC
!
! !REVISION HISTORY:
!  24 Feb 2026 - M. Gooding - Initial version: moved BrC chemistry
!                              out of carbon_mod.F90 into standalone module
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Conversion arrays to pass mass between chemistry steps within
  ! a single timestep.  These are module-level so that CHEM_FSOAS
  ! can write FSOAS_CONV and CHEM_BRCSOA can read it, etc.
  REAL(fp), ALLOCATABLE :: FSOAS_CONV(:,:,:)
  REAL(fp), ALLOCATABLE :: BRCSOA_CONV(:,:,:)

  ! First-call flag for lazy initialisation
  LOGICAL, SAVE :: FIRST = .TRUE.

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chembrc
!
! !DESCRIPTION: Subroutine ChemBrC is the top-level driver for brown carbon
!  chemistry.  It looks up species indices at runtime using Ind\_() so that
!  the code does nothing if BrC species are not defined in the simulation.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE ChemBrC( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, State_Met, RC          )
!
! !USES:
!
   USE ErrCode_Mod
   USE ERROR_MOD,      ONLY : DEBUG_MSG
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Chm_Mod,  ONLY : Ind_
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  24 Feb 2026 - M. Gooding - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER            :: id_FSOAS, id_BRCSOA, id_WTC
   CHARACTER(LEN=255) :: ErrMsg, ThisLoc

   !=================================================================
   ! ChemBrC begins here!
   !=================================================================
   RC      = GC_SUCCESS
   ErrMsg  = ''
   ThisLoc = ' -> at ChemBrC (in module GeosCore/brc_mod.F90)'

   !-----------------------------------------------------------------
   ! Look up species IDs — exit gracefully if not defined
   !-----------------------------------------------------------------
   id_FSOAS  = Ind_('FSOAS' )
   id_BRCSOA = Ind_('BRCSOA')
   id_WTC    = Ind_('WTC'   )

   IF ( id_FSOAS  <= 0 ) RETURN
   IF ( id_BRCSOA <= 0 ) RETURN
   IF ( id_WTC    <= 0 ) RETURN

   !-----------------------------------------------------------------
   ! Lazy initialisation of conversion arrays on first call
   !-----------------------------------------------------------------
   IF ( FIRST ) THEN
      CALL Init_BrC( State_Grid, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Init_BrC"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF
      FIRST = .FALSE.
   ENDIF

   !-----------------------------------------------------------------
   ! Step 1: FSOAS -> BRCSOA  (rapid darkening, tau ~ 0.25 day)
   !-----------------------------------------------------------------
   CALL CHEM_FSOAS( Input_Opt,  State_Chm, State_Diag, &
                    State_Grid, id_FSOAS,  RC          )

   IF ( RC /= GC_SUCCESS ) THEN
      ErrMsg = 'Error encountered in "CHEM_FSOAS"!'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
   ENDIF

   IF ( Input_Opt%Verbose ) THEN
      CALL DEBUG_MSG( '### CHEMBRC: after CHEM_FSOAS' )
   ENDIF

   !-----------------------------------------------------------------
   ! Step 2: Add darkened mass to BRCSOA, then bleach BRCSOA -> WTC
   !-----------------------------------------------------------------
   CALL CHEM_BRCSOA( Input_Opt,  State_Chm, State_Diag, &
                     State_Grid, id_BRCSOA, RC          )

   IF ( RC /= GC_SUCCESS ) THEN
      ErrMsg = 'Error encountered in "CHEM_BRCSOA"!'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
   ENDIF

   IF ( Input_Opt%Verbose ) THEN
      CALL DEBUG_MSG( '### CHEMBRC: after CHEM_BRCSOA' )
   ENDIF

   !-----------------------------------------------------------------
   ! Step 3: Receive bleached mass into WTC
   !-----------------------------------------------------------------
   CALL CHEM_WTC( Input_Opt,  State_Chm, State_Diag, &
                  State_Grid, id_WTC,    RC          )

   IF ( RC /= GC_SUCCESS ) THEN
      ErrMsg = 'Error encountered in "CHEM_WTC"!'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
   ENDIF

   IF ( Input_Opt%Verbose ) THEN
      CALL DEBUG_MSG( '### CHEMBRC: after CHEM_WTC' )
   ENDIF

 END SUBROUTINE ChemBrC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_brc
!
! !DESCRIPTION: Subroutine Init\_BrC allocates and zeroes the module-level
!  conversion arrays used to pass mass between the BrC chemistry steps.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Init_BrC( State_Grid, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
   TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  24 Feb 2026 - M. Gooding - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! Init_BrC begins here!
   !=================================================================
   RC = GC_SUCCESS

   ! Allocate FSOAS_CONV
   ALLOCATE( FSOAS_CONV( State_Grid%NX, &
                          State_Grid%NY, &
                          State_Grid%NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:FSOAS_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   FSOAS_CONV = 0e+0_fp

   ! Allocate BRCSOA_CONV
   ALLOCATE( BRCSOA_CONV( State_Grid%NX, &
                           State_Grid%NY, &
                           State_Grid%NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:BRCSOA_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   BRCSOA_CONV = 0e+0_fp

 END SUBROUTINE Init_BrC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_brc
!
! !DESCRIPTION: Subroutine Cleanup\_BrC deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Cleanup_BrC( RC )
!
! !USES:
!
   USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
   INTEGER, INTENT(OUT) :: RC   ! Success or failure?
!
! !REVISION HISTORY:
!  24 Feb 2026 - M. Gooding - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! Cleanup_BrC begins here!
   !=================================================================
   RC = GC_SUCCESS

   IF ( ALLOCATED( FSOAS_CONV ) ) THEN
      DEALLOCATE( FSOAS_CONV, STAT=RC )
      CALL GC_CheckVar( 'brc_mod.F90:FSOAS_CONV', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
   ENDIF

   IF ( ALLOCATED( BRCSOA_CONV ) ) THEN
      DEALLOCATE( BRCSOA_CONV, STAT=RC )
      CALL GC_CheckVar( 'brc_mod.F90:BRCSOA_CONV', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
   ENDIF

 END SUBROUTINE Cleanup_BrC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_fsoas
!
! !DESCRIPTION: Subroutine CHEM\_FSOAS converts FSOAS to BRCSOA via a
!  first-order darkening process with e-folding time FSOAS\_LIFE days.
!  The converted mass is stored in the module array FSOAS\_CONV for
!  uptake by CHEM\_BRCSOA.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_FSOAS( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, spcId,     RC          )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   INTEGER,        INTENT(IN)    :: spcId        ! FSOAS species Id
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  Drydep is applied in mixing_mod.F90.
!
! !REVISION HISTORY:
!  19 Feb 2026 - M. Gooding - Initial version
!  24 Feb 2026 - M. Gooding - Moved to brc_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I,      J,   L
   REAL(fp)            :: DTCHEM, KFSOAS, FREQ, TC0, CNEW, RKT

   ! Pointers
   REAL(fp), POINTER   :: TC(:,:,:)
!
! !DEFINED PARAMETERS:
!
   ! E-folding lifetime for FSOAS darkening [days]
   REAL(fp), PARAMETER :: FSOAS_LIFE = 0.25e+0_fp

   !=================================================================
   ! CHEM_FSOAS begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize
   KFSOAS    = 1.e+0_fp / ( 86400e+0_fp * FSOAS_LIFE )
   DTCHEM    = GET_TS_CHEM()
   FSOAS_CONV = 0e+0_fp
   TC        => State_Chm%Species(spcId)%Conc

   !=================================================================
   ! Conversion from FSOAS to BRCSOA:
   !   First-order loss with e-folding time FSOAS_LIFE days
   !
   ! Both aerosols are dry-deposited via mixing_mod.F90
   !=================================================================
   !$OMP PARALLEL DO                                                         &
   !$OMP DEFAULT( SHARED                                                    )&
   !$OMP PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW                            )&
   !$OMP COLLAPSE( 3                                                        )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial FSOAS mass [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Amount of FSOAS left after chemistry [kg]
      RKT  = ( KFSOAS + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of FSOAS converted to BRCSOA [kg/timestep]
      FSOAS_CONV(I,J,L) = ( TC0 - CNEW ) &
                         * KFSOAS / ( KFSOAS + FREQ )

      ! Store new concentration back into species array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_FSOAS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_brcsoa
!
! !DESCRIPTION: Subroutine CHEM\_BRCSOA performs two steps:
!  (1) Adds newly-darkened BrC mass from FSOAS (stored in FSOAS\_CONV)
!      into the BRCSOA tracer.
!  (2) Photo-bleaches BRCSOA to WTC using a first-order loss with
!      e-folding time BRCSOA\_LIFE days.
!  The bleached mass is stored in BRCSOA\_CONV for uptake by CHEM\_WTC.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BRCSOA( Input_Opt,  State_Chm, State_Diag, &
                         State_Grid, spcId,     RC          )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   INTEGER,        INTENT(IN)    :: spcId        ! BRCSOA species Id
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  Drydep is applied in mixing_mod.F90.
!
! !REVISION HISTORY:
!  19 Feb 2026 - M. Gooding - Initial version
!  24 Feb 2026 - M. Gooding - Moved to brc_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I,      J,   L
   REAL(fp)            :: DTCHEM, KBRCSOA, FREQ, TC0, CNEW, RKT, CCV

   ! Pointers
   REAL(fp), POINTER   :: TC(:,:,:)
!
! !DEFINED PARAMETERS:
!
   ! E-folding lifetime for BRCSOA photobleaching [days]
   REAL(fp), PARAMETER :: BRCSOA_LIFE = 0.25e+0_fp

   !=================================================================
   ! CHEM_BRCSOA begins here!
   !=================================================================

   ! Assume success
   RC          = GC_SUCCESS

   ! Initialize
   KBRCSOA     = 1.e+0_fp / ( 86400e+0_fp * BRCSOA_LIFE )
   DTCHEM      = GET_TS_CHEM()

   ! IMPORTANT: Do NOT zero FSOAS_CONV here — it was set in CHEM_FSOAS
   !            and we consume it below, then zero it at the end.
   BRCSOA_CONV = 0e+0_fp

   TC          => State_Chm%Species(spcId)%Conc

   !=================================================================
   ! Photo-bleaching from BRCSOA to WTC (first-order):
   !   e-folding time BRCSOA_LIFE days
   !=================================================================
   !$OMP PARALLEL DO                                                         &
   !$OMP DEFAULT( SHARED                                                    )&
   !$OMP PRIVATE( I, J, L, CCV, TC0, FREQ, RKT, CNEW                       )&
   !$OMP COLLAPSE( 3                                                        )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      !==============================================================
      ! 1) Add newly formed BRCSOA from FSOAS (darkening step)
      !==============================================================
      CCV = FSOAS_CONV(I,J,L)

      ! BRCSOA mass available to bleach this timestep [kg]
      TC0 = TC(I,J,L) + CCV

      !==============================================================
      ! 2) Bleach BRCSOA -> WTC as first-order loss over DTCHEM
      !==============================================================

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Remaining BRCSOA after bleaching [kg]
      RKT  = ( KBRCSOA + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount bleached from BRCSOA to WTC [kg/timestep]
      BRCSOA_CONV(I,J,L) = ( TC0 - CNEW ) &
                          * KBRCSOA / ( KBRCSOA + FREQ )

      ! Store updated BRCSOA back into species array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! We have now consumed FSOAS_CONV for this timestep — zero it
   !=================================================================
   FSOAS_CONV = 0e+0_fp

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_BRCSOA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_wtc
!
! !DESCRIPTION: Subroutine CHEM\_WTC receives the bleached mass from
!  BRCSOA (stored in BRCSOA\_CONV) and adds it to the WTC tracer.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_WTC( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, spcId,     RC          )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   INTEGER,        INTENT(IN)    :: spcId        ! WTC species Id
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
   TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  19 Feb 2026 - M. Gooding - Initial version
!  24 Feb 2026 - M. Gooding - Moved to brc_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER  :: I,   J,   L
   REAL(fp) :: TC0, CNEW, CCV

   ! Pointers
   REAL(fp), POINTER :: TC(:,:,:)

   !=================================================================
   ! CHEM_WTC begins here!
   !=================================================================

   ! Assume success
   RC = GC_SUCCESS
   TC => State_Chm%Species(spcId)%Conc

   !$OMP PARALLEL DO                                                         &
   !$OMP DEFAULT( SHARED                                                    )&
   !$OMP PRIVATE( I, J, L, TC0, CCV, CNEW                                  )&
   !$OMP COLLAPSE( 3                                                        )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Current WTC mass [kg]
      TC0 = TC(I,J,L)

      ! Bleached mass arriving from BRCSOA [kg]
      CCV = BRCSOA_CONV(I,J,L)

      ! Add converted mass to WTC
      CNEW = TC0 + CCV

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Store modified concentration back in species array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! Zero BRCSOA_CONV array for next timestep
   !=================================================================
   BRCSOA_CONV = 0e+0_fp

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_WTC
!EOC
END MODULE BRC_MOD