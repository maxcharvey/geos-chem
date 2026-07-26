!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: brc_optics_test_mod.F90
!
! !DESCRIPTION: Module BRC\_OPTICS\_TEST\_MOD defines isolated online-AOD and
!  RRTMG tests for BrC optical size, OM:OC, and hygroscopic routing.  These
!  modes do not alter transported mass, Cloud-J, or aerosol surface area.
!\\
!\\
! !INTERFACE:
!
MODULE BRC_OPTICS_TEST_MOD
!
! !USES:
!
  USE PRECISION_MOD, ONLY : fp, f8

  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS:
!
#ifndef BRC_OPTICS_TEST_MODE
#define BRC_OPTICS_TEST_MODE 0
#endif

  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_MODE = BRC_OPTICS_TEST_MODE

  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_CURRENT   = 0
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_PDER      = 1
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_OMOC      = 2
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_WETDRY    = 3
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_FULL      = 4
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_BASE_PDER = 5
  INTEGER, PARAMETER, PUBLIC :: BRC_OPTICS_FIXED_FULL = 6
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: BrC_Uses_PDER
  PUBLIC :: BrC_Uses_Base_PDER
  PUBLIC :: BrC_Optical_Mass_Parts
  PUBLIC :: BrC_Base_Compatible_PDER
  PUBLIC :: BrC_Get_PDER_Optics
  PUBLIC :: BrC_Optics_Mode_Name
!
! !REVISION HISTORY:
!  26 Jul 2026 - OpenAI Codex - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Uses_PDER
!
! !DESCRIPTION: Return TRUE for modes using size-interpolated BrC optics.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE LOGICAL FUNCTION BrC_Uses_PDER()

    BrC_Uses_PDER = BRC_OPTICS_MODE == BRC_OPTICS_PDER      .OR. &
                    BRC_OPTICS_MODE == BRC_OPTICS_FULL      .OR. &
                    BRC_OPTICS_MODE == BRC_OPTICS_BASE_PDER

  END FUNCTION BrC_Uses_PDER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Uses_Base_PDER
!
! !DESCRIPTION: Return TRUE for modes using reconstructed base-compatible PDER.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE LOGICAL FUNCTION BrC_Uses_Base_PDER()

    BrC_Uses_Base_PDER = BRC_OPTICS_MODE == BRC_OPTICS_FULL .OR. &
                         BRC_OPTICS_MODE == BRC_OPTICS_BASE_PDER

  END FUNCTION BrC_Uses_Base_PDER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Optics_Mode_Name
!
! !DESCRIPTION: Return a stable label for the compiled test mode.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE FUNCTION BrC_Optics_Mode_Name() RESULT( Name )

    CHARACTER(LEN=24) :: Name

    SELECT CASE ( BRC_OPTICS_MODE )
    CASE ( BRC_OPTICS_CURRENT )
       Name = 'current'
    CASE ( BRC_OPTICS_PDER )
       Name = 'pder-current-mass'
    CASE ( BRC_OPTICS_OMOC )
       Name = 'base-omoc-only'
    CASE ( BRC_OPTICS_WETDRY )
       Name = 'base-wetdry-only'
    CASE ( BRC_OPTICS_FULL )
       Name = 'base-compatible-full'
    CASE ( BRC_OPTICS_BASE_PDER )
       Name = 'pder-base-mass'
    CASE ( BRC_OPTICS_FIXED_FULL )
       Name = 'fixed-base-full'
    CASE DEFAULT
       Name = 'invalid'
    END SELECT

  END FUNCTION BrC_Optics_Mode_Name
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Optical_Mass_Parts
!
! !DESCRIPTION: Return wet and dry optical mass for one BrC aerosol bin.
!  Primary bins N=7,8,10 are altered by OM:OC and wet/dry factorial modes.
!  Secondary bins N=6,9 remain wet and additive.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE SUBROUTINE BrC_Optical_Mass_Parts( N,       WetOMOC, DryOMOC, &
                                          MassIn,  MassWet, MassDry  )

    INTEGER,  INTENT(IN)  :: N
    REAL(fp), INTENT(IN)  :: WetOMOC
    REAL(fp), INTENT(IN)  :: DryOMOC
    REAL(fp), INTENT(IN)  :: MassIn
    REAL(fp), INTENT(OUT) :: MassWet
    REAL(fp), INTENT(OUT) :: MassDry

    LOGICAL :: IsPrimary
    REAL(fp) :: DryToWet

    IsPrimary = N == 7 .OR. N == 8 .OR. N == 10
    DryToWet  = 1.0_fp
    IF ( WetOMOC > 0.0_fp ) DryToWet = DryOMOC / WetOMOC

    MassWet = MassIn
    MassDry = 0.0_fp

    IF ( .NOT. IsPrimary ) RETURN

    SELECT CASE ( BRC_OPTICS_MODE )
    CASE ( BRC_OPTICS_OMOC )
       MassWet = 0.5_fp * MassIn * ( 1.0_fp + DryToWet )
    CASE ( BRC_OPTICS_WETDRY )
       MassWet = 0.5_fp * MassIn
       MassDry = 0.5_fp * MassIn
    CASE ( BRC_OPTICS_FULL, BRC_OPTICS_FIXED_FULL )
       MassWet = 0.5_fp * MassIn
       MassDry = 0.5_fp * MassIn * DryToWet
    END SELECT

  END SUBROUTINE BrC_Optical_Mass_Parts
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Base_Compatible_PDER
!
! !DESCRIPTION: Reconstruct PDER using base-compatible primary BrC OM mass.
!  BRCSOA and FSOAS remain additive secondary aerosol.  WTC is provisionally
!  treated as primary because its primary/secondary provenance is not retained.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE REAL(fp) FUNCTION BrC_Base_Compatible_PDER( SNA,      OCPO,   OCPISOA, &
                                                   BRCSOA,   NPBRC,  WTC,     &
                                                   FSOAS,    PBRC,   WetOMOC, &
                                                   DryOMOC                    )

    REAL(fp), INTENT(IN) :: SNA
    REAL(fp), INTENT(IN) :: OCPO
    REAL(fp), INTENT(IN) :: OCPISOA
    REAL(fp), INTENT(IN) :: BRCSOA
    REAL(fp), INTENT(IN) :: NPBRC
    REAL(fp), INTENT(IN) :: WTC
    REAL(fp), INTENT(IN) :: FSOAS
    REAL(fp), INTENT(IN) :: PBRC
    REAL(fp), INTENT(IN) :: WetOMOC
    REAL(fp), INTENT(IN) :: DryOMOC

    REAL(fp) :: OM
    REAL(fp) :: PrimaryScale
    REAL(fp) :: Ratio
    REAL(fp) :: SNAOM

    BrC_Base_Compatible_PDER = 0.005_fp
    IF ( SNA <= 0.0_fp ) RETURN

    PrimaryScale = 1.0_fp
    IF ( WetOMOC > 0.0_fp ) THEN
       PrimaryScale = 0.5_fp * ( WetOMOC + DryOMOC ) / WetOMOC
    ENDIF

    OM = OCPO + OCPISOA + BRCSOA + FSOAS + &
         ( NPBRC + WTC + PBRC ) * PrimaryScale
    IF ( OM <= 0.0_fp ) RETURN

    SNAOM = ( SNA + OM ) * 1.0e+9_fp
    Ratio = OM / SNA
    IF ( SNAOM <= 0.0_fp .OR. Ratio <= 0.0_fp ) RETURN

    BrC_Base_Compatible_PDER = &
         EXP( 4.36_fp + 0.20_fp * LOG( SNAOM ) + &
              0.065_fp * LOG( Ratio ) ) * 0.001_fp / 0.9_fp

    IF ( BrC_Base_Compatible_PDER <= 0.0_fp ) THEN
       BrC_Base_Compatible_PDER = 0.005_fp
    ENDIF

  END FUNCTION BrC_Base_Compatible_PDER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BrC_Get_PDER_Optics
!
! !DESCRIPTION: Interpolate one aerosol lookup table over dry effective radius
!  and relative humidity.  Return the dry radius/Q, hygroscopic AOD scale, wet
!  single-scattering albedo, and wet asymmetry factor.
!EOP
!------------------------------------------------------------------------------
!BOC
  PURE SUBROUTINE BrC_Get_PDER_Optics( PDER,     RelHum, RH,      Radius, &
                                       Q,        SSA,    Asym,    Rdry,   &
                                       Qdry,     ScaleOD, SSADry, AsymDry, &
                                       SSAWet,   AsymWet                    )

    REAL(fp), INTENT(IN)  :: PDER
    REAL(fp), INTENT(IN)  :: RelHum
    REAL(fp), INTENT(IN)  :: RH(:)
    REAL(f8), INTENT(IN)  :: Radius(:,:)
    REAL(f8), INTENT(IN)  :: Q(:,:)
    REAL(f8), INTENT(IN)  :: SSA(:,:)
    REAL(f8), INTENT(IN)  :: Asym(:,:)
    REAL(f8), INTENT(OUT) :: Rdry
    REAL(f8), INTENT(OUT) :: Qdry
    REAL(f8), INTENT(OUT) :: ScaleOD
    REAL(f8), INTENT(OUT) :: SSADry
    REAL(f8), INTENT(OUT) :: AsymDry
    REAL(f8), INTENT(OUT) :: SSAWet
    REAL(f8), INTENT(OUT) :: AsymWet

    INTEGER  :: G
    INTEGER  :: IRH
    INTEGER  :: R
    INTEGER  :: NRHLocal
    REAL(f8) :: Fraction
    REAL(f8) :: QWet
    REAL(f8) :: RadiusWet
    REAL(f8) :: RW( SIZE( RH ) )
    REAL(f8) :: QW( SIZE( RH ) )
    REAL(f8) :: SSW( SIZE( RH ) )
    REAL(f8) :: AsymW( SIZE( RH ) )

    NRHLocal = SIZE( RH )
    G        = 1

    DO WHILE ( REAL( PDER, f8 ) > Radius(1,G) .AND. G < SIZE( Radius, 2 ) )
       G = G + 1
    ENDDO

    IF ( G == 1 ) THEN
       RW    = Radius(:,G)
       QW    = Q(:,G)
       SSW   = SSA(:,G)
       AsymW = Asym(:,G)
    ELSE
       Fraction = ( REAL( PDER, f8 ) - Radius(1,G-1) ) / &
                  ( Radius(1,G) - Radius(1,G-1) )
       Fraction = MIN( Fraction, 1.0_f8 )
       DO R = 1, NRHLocal
          RW(R)    = Fraction * Radius(R,G) + ( 1.0_f8 - Fraction ) * Radius(R,G-1)
          QW(R)    = Fraction * Q(R,G)      + ( 1.0_f8 - Fraction ) * Q(R,G-1)
          SSW(R)   = Fraction * SSA(R,G)    + ( 1.0_f8 - Fraction ) * SSA(R,G-1)
          AsymW(R) = Fraction * Asym(R,G)   + ( 1.0_f8 - Fraction ) * Asym(R,G-1)
       ENDDO
    ENDIF

    IF      ( RelHum <= RH(2) ) THEN
       IRH = 1
    ELSE IF ( RelHum <= RH(3) ) THEN
       IRH = 2
    ELSE IF ( RelHum <= RH(4) ) THEN
       IRH = 3
    ELSE IF ( RelHum <= RH(5) ) THEN
       IRH = 4
    ELSE
       IRH = NRHLocal
    ENDIF

    IF ( IRH == NRHLocal ) THEN
       RadiusWet = RW(NRHLocal)
       QWet      = QW(NRHLocal)
       SSAWet    = SSW(NRHLocal)
       AsymWet   = AsymW(NRHLocal)
    ELSE
       Fraction = ( REAL( RelHum, f8 ) - REAL( RH(IRH), f8 ) ) / &
                  REAL( RH(IRH+1) - RH(IRH), f8 )
       Fraction = MIN( Fraction, 1.0_f8 )
       RadiusWet = Fraction * RW(IRH+1)    + ( 1.0_f8 - Fraction ) * RW(IRH)
       QWet      = Fraction * QW(IRH+1)    + ( 1.0_f8 - Fraction ) * QW(IRH)
       SSAWet    = Fraction * SSW(IRH+1)   + ( 1.0_f8 - Fraction ) * SSW(IRH)
       AsymWet   = Fraction * AsymW(IRH+1) + ( 1.0_f8 - Fraction ) * AsymW(IRH)
    ENDIF

    Rdry = RW(1)
    Qdry = QW(1)
    SSADry  = SSW(1)
    AsymDry = AsymW(1)

    IF ( Rdry > 0.0_f8 .AND. Qdry /= 0.0_f8 ) THEN
       ScaleOD = ( QWet / Qdry ) * ( RadiusWet / Rdry )**2
    ELSE
       ScaleOD = 0.0_f8
    ENDIF

  END SUBROUTINE BrC_Get_PDER_Optics
!EOC

END MODULE BRC_OPTICS_TEST_MOD
