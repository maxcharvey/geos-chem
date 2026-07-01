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
!    FSOAP ---(gas-to-particle)---> FSOAS ---(rapid darkening)--->
!    BRCSOA ---(photobleaching)---> WTC
!
!  In parallel, non-persistent BrC-POA (directly emitted from fires)
!  also photo-bleaches to WTC via the same viscosity-dependent scheme,
!  while persistent primary BrC-POA is carried as a non-bleaching tracer:
!
!    NPBRCPOA ---(photobleaching)---> WTC
!    PBRCPOA  ---(persistent tracer)---> PBRCPOA
!
!  The photobleaching rate (BRCSOA -> WTC) is parameterised as a function
!  of local temperature and relative humidity following the viscosity-
!  dependent kinetic framework of:
!
!    Schnitzler, E.G. et al. (2022), "Rate of atmospheric brown carbon
!    whitening governed by environmental conditions", PNAS, 119(38),
!    e2205610119. https://doi.org/10.1073/pnas.2205610119
!
!  The chain of calculation is:
!    (T, RH) -> viscosity eta  [VFT + Arrhenius mixing rule, SI Eq. S1-S6]
!            -> D_O3           [fractional Stokes-Einstein, SI Eq. S12-S13]
!            -> tau_BrC        [resistor-model lifetime, Main Eq. 3]
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
  PRIVATE :: CHEM_FSOAP
  PRIVATE :: CHEM_FSOAS
  PRIVATE :: CHEM_BRCSOA
  PRIVATE :: CHEM_NPBRCPOA
  PRIVATE :: CHEM_WTC
  PRIVATE :: VISC_BBOA_FUNC
  PRIVATE :: VISC_WATER_FUNC
  PRIVATE :: CALC_TAU_BRC

! !REVISION HISTORY:
!  24 Feb 2026 - M. Harvey - Initial version: BrC chemistry module
!  24 Feb 2026 - M. Harvey - Added Schnitzler et al. (2022) viscosity-
!                             dependent photobleaching parameterisation
!  25 Feb 2026 - M. Harvey - Added FSOAP gas-phase precursor pathway
!  25 Feb 2026 - M. Harvey - Added development diagnostics arrays
!  25 Feb 2026 - M. Harvey - Added NPBRCPOA (non-persistent BrC-POA)
!  25 Feb 2026 - M. Harvey - Updated FSOAS darkening lifetime to 1 day
!  11 Mar 2026 - M. Harvey - Replaced fixed P_O3_ATM with local O3 from
!                             State_Chm; CALC_TAU_BRC now takes P_O3 arg
!  19 Mar 2026 - M. Harvey - TAU_MIN changed from 3600 s to 21600 s (6 h)
!  19 Mar 2026 - M. Harvey - Added 25% stop-loss floor (FRAC_PERM) to
!                             BRCSOA and NPBRCPOA bleaching
!  19 Mar 2026 - M. Harvey - Added BLEACH_SCHEME runtime switch (0-4)
!                             for selecting bleaching parameterisation
!  30 Jun 2026 - M. Harvey - Removed timestep stop-loss; persistent
!                             primary BrC is now represented by PBRCPOA
!                             at emission time
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Small number to prevent underflow
  REAL(fp), PARAMETER :: SMALLNUM = 1e-20_fp

  ! Conversion arrays to pass mass between chemistry steps within
  ! a single timestep.  These are module-level so that CHEM_FSOAP
  ! can write FSOAP_CONV and CHEM_FSOAS can read it, etc.
  REAL(fp), ALLOCATABLE :: FSOAP_CONV(:,:,:)
  REAL(fp), ALLOCATABLE :: FSOAS_CONV(:,:,:)
  REAL(fp), ALLOCATABLE :: BRCSOA_CONV(:,:,:)
  REAL(fp), ALLOCATABLE :: NPBRC_CONV(:,:,:)

  ! First-call flag for lazy initialisation
  LOGICAL, SAVE :: FIRST = .TRUE.

  !=========================================================================
  ! Physical constants
  !=========================================================================
  REAL(fp), PARAMETER :: PI      = 3.14159265358979_fp
  REAL(fp), PARAMETER :: KB      = 1.380649e-23_fp   ! Boltzmann [J/K]
  REAL(fp), PARAMETER :: HPA2ATM = 1.0_fp / 1013.25_fp ! hPa -> atm

  !=========================================================================
  ! Parameters for Schnitzler et al. (2022) viscosity parameterisation
  ! Source: SI Appendix, Section S1 (Eq. S1-S6)
  !=========================================================================

  !--- Molecular properties of BBOA ---
  ! MW_BBOA: average molecular weight [g/mol]
  !   Based on BrC chromophores 138-358 g/mol in pine BBOA
  !   (Fleming et al. 2020, Atmos. Chem. Phys. 20, 1105-1129)
  REAL(fp), PARAMETER :: MW_BBOA = 248.0_fp

  ! MW_H2O: molecular weight of water [g/mol]
  REAL(fp), PARAMETER :: MW_H2O  = 18.015_fp

  ! MW_O3: molecular weight of ozone [g/mol]
  REAL(fp), PARAMETER :: MW_O3   = 48.0_fp

  ! MW_AIR: mean molecular weight of dry air [g/mol]
  REAL(fp), PARAMETER :: MW_AIR  = 28.97_fp

  !--- Hygroscopicity ---
  ! KAPPA: mass-based hygroscopicity parameter (SI Eq. S4)
  !   Fitted to poke-flow viscosity data (Fig. 2A); back-calculated
  !   from the fit curve to reproduce eta ~ 10^3 Pa s at 25% RH.
  !   Consistent with published BBOA kappa range 0.05-0.15
  !   (Rothfuss & Petters 2017, PCCP 19, 6532-6545)
  REAL(fp), PARAMETER :: KAPPA   = 0.07_fp

  !--- Reference viscosities at T_REF = 294 K ---
  ! ETA_DRY_294: viscosity of dry BBOA at 294 K [Pa s]
  !   From poke-flow measurements at 0% RH (Fig. 2A, SI Fig. S6):
  !   ~1-2 x 10^5 Pa s
  REAL(fp), PARAMETER :: ETA_DRY_294 = 1.0e+5_fp

  ! LOG10 of dry viscosity at 294 K for Arrhenius mixing rule (SI Eq. S2)
  REAL(fp), PARAMETER :: LOG_ETA_DRY = 5.0_fp         ! log10(1e5)

  ! LOG10 of water viscosity [Pa s] for mixing rule
  !   eta_H2O = 10^-3 Pa s (Korson et al. 1969, J. Phys. Chem. 73, 34)
  REAL(fp), PARAMETER :: LOG_ETA_H2O = -3.0_fp         ! log10(1e-3)

  ! Reference temperature for viscosity data [K]
  REAL(fp), PARAMETER :: T_REF       = 294.0_fp

  !--- Vogel-Fulcher-Tamman (VFT) parameters (SI Eq. S5) ---
  ! ETA_INF: viscosity at infinite temperature [Pa s]
  !   (Angell 2002, Chem. Rev. 102, 2627; Angell 1991, J. Non-Cryst.
  !   Solids 131, 13)
  REAL(fp), PARAMETER :: ETA_INF = 1.0e-5_fp

  ! DF: fragility parameter (dimensionless)
  !   (Shiraiwa et al. 2017, Nat. Commun. 8, 15002;
  !    DeRieux et al. 2018, ACP 18, 6331)
  REAL(fp), PARAMETER :: DF      = 10.0_fp

  !--- Water viscosity as f(T) (SI Eq. S14) ---
  !   log10(eta_water) = AW_VFT + BW_VFT / (T - T0W_VFT)
  !   (Maclean et al. 2021, ACS Earth Space Chem. 5, 3458)
  REAL(fp), PARAMETER :: AW_VFT  = -4.28_fp
  REAL(fp), PARAMETER :: BW_VFT  = 152.87_fp
  REAL(fp), PARAMETER :: T0W_VFT = 173.06_fp

  !=========================================================================
  ! Parameters for fractional Stokes-Einstein equation
  ! Source: SI Appendix, Section S3 (Eq. S12-S14)
  !
  ! D_O3(RH,T) = D_O3_water(T) * (eta_water(T) / eta_BBOA(RH,T))^xi
  !
  ! where D_O3_water(T) = k_B*T / (6*pi*eta_water(T)*R_O3)
  !       is the standard Stokes-Einstein diffusion in pure water.
  !=========================================================================

  ! XI_O3: fractional exponent for O3 diffusion in BBOA
  !   From SI Eq. S13 (Evoy et al. 2020, J. Phys. Chem. A 124, 2301):
  !     xi = 1 - A*exp(-B * R_O3 / R_matrix)
  !     A = 0.73, B = 1.79
  !     R_O3    = 0.198 nm  (van der Waals radius)
  !     R_matrix = 0.423 nm (from MW=248 g/mol, density=1.3 g/cm3,
  !                           assuming spherical molecules)
  !     => xi = 1 - 0.73*exp(-1.79 * 0.198/0.423) = 0.684
  REAL(fp), PARAMETER :: XI_O3   = 0.684_fp

  ! R_O3: hydrodynamic radius of ozone [m]
  !   Based on van der Waals radii using atomic increments
  !   (Evoy et al. 2020)
  REAL(fp), PARAMETER :: R_O3    = 0.198e-9_fp

  !=========================================================================
  ! Parameters for BrC whitening lifetime
  ! Source: Main text Eq. 3, SI Appendix Fig. S7
  !
  ! tau_BrC = C_FACTOR * 2*a / (3 * Hk(T) * P_O3 * sqrt(D_O3))
  !
  ! where Hk(T) = H * sqrt(k2/[BrC]_0) is nearly T-independent
  !       (~9 atm^-1 s^-1/2), with a weak linear fit from SI Fig. S7.
  !=========================================================================

  ! Hk(T) = HK_SLOPE * T + HK_INTERCEPT   [atm^-1 s^-1/2]
  !   Linear fit from SI Fig. S7:  y = -3.50e-4 * x + 8.84
  REAL(fp), PARAMETER :: HK_SLOPE     = -3.50e-4_fp
  REAL(fp), PARAMETER :: HK_INTERCEPT = 8.84_fp

  ! A_PARTICLE: assumed BBOA particle radius [m]
  !   Based on volume size distributions of atmospheric BBOA
  !   (Shi et al. 2019; Yu et al. 2013; Ogunjobi et al. 2004)
  !   Main text p.4: a = 150 nm
  REAL(fp), PARAMETER :: A_PARTICLE   = 150.0e-9_fp

  ! P_O3_ATM: assumed tropospheric O3 partial pressure [atm]
  !   35 ppb globally averaged, multiplied by mean pressure
  !   (Main text p.4; Tarasick et al. 2019; Lefohn et al. 2018)
  REAL(fp), PARAMETER :: P_O3_ATM     = 3.5e-8_fp

  ! C_FACTOR: prefactor in lifetime equation (Main Eq. 3)
  !   = 1 - 1/sqrt(e) = 1 - e^(-0.5) ≈ 0.3935
  REAL(fp), PARAMETER :: C_FACTOR     = 0.39346934_fp

  ! Maximum allowed lifetime [s] (~3.2 years; effectively no whitening)
  ! Changed value to be 1e10 s (317 years) for testing
  REAL(fp), PARAMETER :: TAU_MAX      = 1.0e+10_fp

  ! Minimum allowed lifetime [s]
  ! Set to 6 hours based on observational constraints on fastest
  ! plausible bleaching in warm, humid boundary layer conditions
  REAL(fp), PARAMETER :: TAU_MIN      = 21600.0_fp

  ! Maximum viscosity [Pa s] (glass transition cutoff; SI Fig. S9)
  !   Viscosities above 10^12 Pa s correspond to a glass state and
  !   are not modelled well by VFT (Fig. 2B caption)
  REAL(fp), PARAMETER :: ETA_MAX      = 1.0e+12_fp
  
  !=======================================================================
  !Conversion parameters 
  !=======================================================================

  ! OM:OC ratio for fresh fire SOA represented by FSOAS [kg_OM kgC-1].
  ! FSOAS is carried as OM mass like SOAS; conversion to BRCSOA divides
  ! by this value to store BRCSOA as carbon mass.
  !   Typical range 1.6-2.1 for fresh BBOA
  !   (Aiken et al. 2008, Environ. Sci. Technol. 42, 4478;
  !    Turpin & Lim 2001, Aerosol Sci. Technol. 35, 602)
  REAL(fp), PARAMETER :: OMOC_BBOA = 1.8_fp

  !=========================================================================
  ! Bleaching scheme selector
  !=========================================================================

  ! BrC_Bleach_Scheme in geoschem_config.yml selects photobleaching:
  !   0 = No bleaching (BRCSOA/NPBRCPOA persist indefinitely)
  !   1 = Fixed 1-day lifetime everywhere
  !   2 = Fixed 1-day lifetime below 1 km AGL only
  !   3 = Viscosity-dependent (Schnitzler et al. 2022), fixed 35 ppb O3
  !   4 = Viscosity-dependent (Schnitzler et al. 2022), local O3 (default)

  ! Fixed bleaching lifetime [s] for schemes 1 and 2
  REAL(fp), PARAMETER :: TAU_1DAY = 86400.0_fp

  ! Altitude threshold [m] for scheme 2 (bleach below this height AGL)
  REAL(fp), PARAMETER :: ALT_THRESH = 1000.0_fp
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: visc_water_func
!
! !DESCRIPTION: Returns the dynamic viscosity of pure water [Pa s] at
!  temperature T [K], using the VFT-form parameterisation from
!  Maclean et al. (2021), as cited in Schnitzler et al. (2022) SI Eq. S14:
!
!    log10(eta) = A + B / (T - T0)
!
!  where A = -4.28, B = 152.87, T0 = 173.06 K.
!\\
!\\
! !INTERFACE:
!
 FUNCTION VISC_WATER_FUNC( T ) RESULT( ETA_W )
!
! !INPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: T        ! Temperature [K]
!
! !RETURN VALUE:
!
   REAL(fp)             :: ETA_W    ! Viscosity of water [Pa s]
!EOP
!------------------------------------------------------------------------------
!BOC

   ! Guard against T <= T0 (unphysical regime)
   IF ( T <= T0W_VFT + 1.0_fp ) THEN
      ETA_W = 1.0e+2_fp   ! Very high but below glass transition
   ELSE
      ETA_W = 10.0_fp**( AW_VFT + BW_VFT / ( T - T0W_VFT ) )
   ENDIF

 END FUNCTION VISC_WATER_FUNC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: visc_bboa_func
!
! !DESCRIPTION: Returns the dynamic viscosity of BBOA [Pa s] at
!  temperature T [K] and water activity a\_w (= RH/100).
!
!  The calculation follows two steps from Schnitzler et al. (2022) SI:
!
!  Step 1 - Viscosity at T\_REF = 294 K via Arrhenius mixing rule:
!    (a) Compute BBOA mass fraction from kappa-Koehler (SI Eq. S4):
!        w_s = (1 + kappa * a_w / (1 - a_w))^(-1)
!    (b) Convert to BBOA mole fraction (SI Eq. S3):
!        chi = (w/MW_BBOA) / (w/MW_BBOA + (1-w)/MW_H2O)
!    (c) Apply Arrhenius mixing rule (SI Eq. S2):
!        log10(eta_mix) = chi*log10(eta_dry) + (1-chi)*log10(eta_H2O)
!
!  Step 2 - Extend to arbitrary T via VFT equation (SI Eq. S5-S6):
!    (a) Compute Vogel temperature T0(RH) from Eq. S6:
!        T0 = ln(eta_294/eta_inf)*294 / (Df + ln(eta_294/eta_inf))
!    (b) Apply VFT (Eq. S5):
!        eta(T) = eta_inf * exp(T0 * Df / (T - T0))
!\\
!\\
! !INTERFACE:
!
 FUNCTION VISC_BBOA_FUNC( T, AW ) RESULT( ETA )
!
! !INPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: T     ! Temperature [K]
   REAL(fp), INTENT(IN) :: AW    ! Water activity [0-1] (= RH/100)
!
! !RETURN VALUE:
!
   REAL(fp)             :: ETA   ! BBOA viscosity [Pa s]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: AW_SAFE         ! Clamped water activity
   REAL(fp) :: W_BBOA          ! BBOA mass fraction [0-1]
   REAL(fp) :: CHI_BBOA        ! BBOA mole fraction [0-1]
   REAL(fp) :: LOG_ETA_294     ! log10(viscosity) at 294 K
   REAL(fp) :: ETA_294         ! Viscosity at 294 K [Pa s]
   REAL(fp) :: LN_RATIO        ! ln(eta_294 / eta_inf)
   REAL(fp) :: T0              ! Vogel temperature [K]

   !=================================================================
   ! Step 1: Viscosity at T_REF = 294 K as a function of RH
   !=================================================================

   ! Clamp water activity to avoid division by zero
   AW_SAFE = MAX( AW, 0.0_fp )
   ! Altering for testing to match the Schnitzler et al. (2022) range
   AW_SAFE = MIN( AW_SAFE, 0.90_fp )

   ! BBOA mass fraction from kappa-Koehler (SI Eq. S4)
   !   w_s = (1 + kappa * a_w / (1 - a_w))^(-1)
   W_BBOA = 1.0_fp / ( 1.0_fp                                     &
          + KAPPA * AW_SAFE / ( 1.0_fp - AW_SAFE ) )

   ! BBOA mole fraction (SI Eq. S3)
   !   chi = (w/MW_BBOA) / (w/MW_BBOA + (1-w)/MW_H2O)
   CHI_BBOA = ( W_BBOA / MW_BBOA )                                &
            / ( W_BBOA / MW_BBOA                                  &
              + ( 1.0_fp - W_BBOA ) / MW_H2O )

   ! Arrhenius mixing rule (SI Eq. S2)
   !   log10(eta_mix) = chi*log10(eta_dry) + (1-chi)*log10(eta_H2O)
   LOG_ETA_294 = CHI_BBOA         * LOG_ETA_DRY                   &
               + ( 1.0_fp - CHI_BBOA ) * LOG_ETA_H2O

   ETA_294 = 10.0_fp**LOG_ETA_294

   !=================================================================
   ! Step 2: Extend to arbitrary T using VFT (SI Eq. S5-S6)
   !=================================================================

   ! Compute Vogel temperature T0 from SI Eq. S6:
   !   T0 = ln(eta_294 / eta_inf) * T_REF / (Df + ln(eta_294 / eta_inf))
   LN_RATIO = LOG( ETA_294 / ETA_INF )  ! natural log

   IF ( LN_RATIO <= 0.0_fp ) THEN
      ! If eta_294 <= eta_inf, no VFT correction needed
      ETA = ETA_294
      RETURN
   ENDIF

   T0 = LN_RATIO * T_REF / ( DF + LN_RATIO )

   ! VFT equation (SI Eq. S5):
   !   eta(RH,T) = eta_inf * exp( T0(RH) * Df / (T - T0(RH)) )
   IF ( T <= T0 + 1.0_fp ) THEN
      ! Below or at Vogel temperature: glass state
      ETA = ETA_MAX
   ELSE
      ETA = ETA_INF * EXP( T0 * DF / ( T - T0 ) )
   ENDIF

   ! Cap at glass transition
   ETA = MIN( ETA, ETA_MAX )

 END FUNCTION VISC_BBOA_FUNC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_tau_brc
!
! !DESCRIPTION: Returns the e-folding lifetime [seconds] for BrC
!  photobleaching (BRCSOA -> WTC) at given temperature T and water
!  activity AW, using the resistor-model framework of Schnitzler et al.
!  (2022), Main text Eq. 3:
!
!    tau = (1 - 1/sqrt(e)) * 2*a / (3 * Hk(T) * P_O3 * sqrt(D_O3))
!
!  where:
!    Hk(T) = H * sqrt(k2/[BrC]_0) from SI Fig. S7 linear fit
!    D_O3  = D_O3_water * (eta_water / eta_BBOA)^xi   [SI Eq. S12]
!    D_O3_water = k_B*T / (6*pi*eta_water*R_O3)       [Stokes-Einstein]
!    P_O3  = local ozone partial pressure [atm] (passed by caller)
!
!  The viscosity dependence means:
!    - Near surface (warm, humid): tau ~ hours to 1 day
!    - Free troposphere (cold, dry): tau >> 1 week (glass state)
!  This reproduces Fig. 3A: 1-day contour at ~1 km altitude.
!\\
!\\
! !INTERFACE:
!
 FUNCTION CALC_TAU_BRC( T, AW, P_O3 ) RESULT( TAU )
!
! !INPUT PARAMETERS:
!
   REAL(fp), INTENT(IN) :: T     ! Temperature [K]
   REAL(fp), INTENT(IN) :: AW    ! Water activity [0-1] (= RH/100)
   REAL(fp), INTENT(IN) :: P_O3  ! Local O3 partial pressure [atm]
!
! !RETURN VALUE:
!
   REAL(fp)             :: TAU   ! BrC bleaching lifetime [s]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   REAL(fp) :: ETA_BBOA      ! BBOA viscosity [Pa s]
   REAL(fp) :: ETA_W         ! Water viscosity [Pa s]
   REAL(fp) :: D_O3_WATER    ! D_O3 in pure water [m2/s]
   REAL(fp) :: D_O3          ! D_O3 in BBOA [m2/s]
   REAL(fp) :: HK            ! H*sqrt(k2/[BrC]_0) [atm^-1 s^-1/2]
   REAL(fp) :: DENOM         ! Denominator of Eq. 3

   !=================================================================
   ! Step 1: Compute BBOA viscosity and water viscosity at local T
   !=================================================================
   ETA_BBOA = VISC_BBOA_FUNC( T, AW )
   ETA_W    = VISC_WATER_FUNC( T )

   !=================================================================
   ! Step 2: Diffusion coefficient of O3 in BBOA
   !
   ! Fractional Stokes-Einstein (SI Eq. S12):
   !   D_O3 = D_O3_water(T) * (eta_water(T) / eta_BBOA(RH,T))^xi
   !
   ! where D_O3_water(T) = k_B*T / (6*pi*eta_water*R_O3)
   !   is the standard Stokes-Einstein diffusion coefficient for
   !   ozone in pure water.
   !
   ! The fractional exponent xi = 0.684 accounts for the fact that
   ! small molecules (O3) diffuse faster than predicted by standard
   ! Stokes-Einstein in viscous media (Evoy et al. 2020).
   !=================================================================

   ! D_O3 in pure water via standard Stokes-Einstein [m2/s]
   D_O3_WATER = KB * T / ( 6.0_fp * PI * ETA_W * R_O3 )

   ! Fractional SE: scale by viscosity ratio raised to xi
   IF ( ETA_BBOA < 1.0e-10_fp ) THEN
      ! Safety: BBOA less viscous than water (shouldn't happen)
      D_O3 = D_O3_WATER
   ELSE
      D_O3 = D_O3_WATER * ( ETA_W / ETA_BBOA )**XI_O3
   ENDIF

   ! Floor D_O3 to prevent division by zero in sqrt below
   D_O3 = MAX( D_O3, 1.0e-30_fp )

   !=================================================================
   ! Step 3: BrC photobleaching lifetime (Main text Eq. 3)
   !
   ! tau = C_FACTOR * 2*a / (3 * Hk * P_O3 * sqrt(D_O3))
   !
   ! Hk(T) = H*sqrt(k2/[BrC]_0) from linear fit of SI Fig. S7:
   !   Hk = -3.50e-4 * T + 8.84   [atm^-1 s^-1/2]
   !   ~9 atm^-1 s^-1/2 across 253-293 K (weak T dependence)
   !=================================================================

   ! Hk as function of temperature
   HK = HK_SLOPE * T + HK_INTERCEPT
   HK = MAX( HK, 1.0e-10_fp )   ! Safety floor

   ! Denominator: 3 * Hk * P_O3 * sqrt(D_O3)
   DENOM = 3.0_fp * HK * P_O3 * SQRT( D_O3 )

   IF ( DENOM < 1.0e-30_fp ) THEN
      TAU = TAU_MAX
   ELSE
      TAU = C_FACTOR * 2.0_fp * A_PARTICLE / DENOM
   ENDIF

   ! Clamp to reasonable range
   TAU = MIN( TAU, TAU_MAX )
   TAU = MAX( TAU, TAU_MIN )

 END FUNCTION CALC_TAU_BRC
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
!
!  The full chain is:
!    Step 0: FSOAP -> FSOAS   (gas-to-particle, tau ~ 1 day)
!    Step 1: FSOAS -> BRCSOA  (rapid darkening, tau ~ 1 day)
!    Step 2: BRCSOA -> WTC    (viscosity-dependent photobleaching)
!    Step 2b: NPBRCPOA -> WTC (viscosity-dependent photobleaching)
!    Step 2c: PBRCPOA persists as emitted primary BrC-POA
!    Step 3: WTC receives bleached mass from BRCSOA and NPBRCPOA
!
!  FSOAP and NPBRCPOA are optional: if not defined in the simulation,
!  their steps are skipped and the remaining chain operates as before.
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
!  24 Feb 2026 - M. Harvey - Initial version
!  25 Feb 2026 - M. Harvey - Added FSOAP step (optional)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER            :: id_FSOAP, id_FSOAS, id_BRCSOA, id_WTC
   INTEGER            :: id_NPBRCPOA
   CHARACTER(LEN=255) :: ErrMsg, ThisLoc

   !=================================================================
   ! ChemBrC begins here!
   !=================================================================
   RC      = GC_SUCCESS
   ErrMsg  = ''
   ThisLoc = ' -> at ChemBrC (in module GeosCore/brc_mod.F90)'

   !-----------------------------------------------------------------
   ! Look up species IDs - exit gracefully if not defined
   ! FSOAP and NPBRCPOA are optional; FSOAS, BRCSOA, WTC are required
   !-----------------------------------------------------------------
   id_FSOAP    = Ind_('FSOAP'   )
   id_FSOAS    = Ind_('FSOAS'   )
   id_BRCSOA   = Ind_('BRCSOA'  )
   id_WTC      = Ind_('WTC'     )
   id_NPBRCPOA = Ind_('NPBRCPOA')

   ! Required species: if any missing, return silently
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
   ! Step 0 (optional): FSOAP -> FSOAS  (gas-to-particle, tau ~ 1 d)
   !   Only runs if FSOAP is defined in the simulation.
   !   If FSOAP is not defined, FSOAS receives mass from direct
   !   emissions only (backward-compatible).
   !-----------------------------------------------------------------
   IF ( id_FSOAP > 0 ) THEN

      CALL CHEM_FSOAP( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, id_FSOAP,  RC          )

      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "CHEM_FSOAP"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( Input_Opt%Verbose ) THEN
         CALL DEBUG_MSG( '### CHEMBRC: after CHEM_FSOAP' )
      ENDIF

   ENDIF

   !-----------------------------------------------------------------
   ! Step 1: FSOAS -> BRCSOA  (rapid darkening, tau ~ 0.25 day)
   !   Also receives condensed mass from FSOAP (if present)
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
   !         Bleaching rate depends on local T and RH
   !-----------------------------------------------------------------
   CALL CHEM_BRCSOA( Input_Opt,  State_Chm, State_Diag,  &
                     State_Grid, State_Met, id_BRCSOA, RC )

   IF ( RC /= GC_SUCCESS ) THEN
      ErrMsg = 'Error encountered in "CHEM_BRCSOA"!'
      CALL GC_Error( ErrMsg, RC, ThisLoc )
      RETURN
   ENDIF

   IF ( Input_Opt%Verbose ) THEN
      CALL DEBUG_MSG( '### CHEMBRC: after CHEM_BRCSOA' )
   ENDIF

   !-----------------------------------------------------------------
   ! Step 2b (optional): NPBRCPOA -> WTC  (viscosity-dependent)
   !   Non-persistent BrC-POA bleaches to WTC using the same
   !   Schnitzler et al. (2022) parameterisation as BRCSOA.
   !   Only runs if NPBRCPOA is defined in the simulation.
   !-----------------------------------------------------------------
   IF ( id_NPBRCPOA > 0 ) THEN

      CALL CHEM_NPBRCPOA( Input_Opt,  State_Chm, State_Diag,  &
                          State_Grid, State_Met, id_NPBRCPOA, RC )

      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "CHEM_NPBRCPOA"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( Input_Opt%Verbose ) THEN
         CALL DEBUG_MSG( '### CHEMBRC: after CHEM_NPBRCPOA' )
      ENDIF

   ENDIF

   ! PBRCPOA is a persistent primary BrC-POA tracer.  It receives the
   ! Forrister-style persistent emission fraction in HEMCO and is not
   ! bleached here.

   !-----------------------------------------------------------------
   ! Step 3: Receive bleached mass into WTC
   !         (from both BRCSOA and NPBRCPOA)
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
!  conversion arrays and diagnostic arrays used by the BrC chemistry.
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
!  24 Feb 2026 - M. Harvey - Initial version
!  25 Feb 2026 - M. Harvey - Added FSOAP_CONV and diagnostic arrays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER :: NX, NY, NZ

   !=================================================================
   ! Init_BrC begins here!
   !=================================================================
   RC = GC_SUCCESS
   NX = State_Grid%NX
   NY = State_Grid%NY
   NZ = State_Grid%NZ

   !-----------------------------------------------------------------
   ! Conversion arrays
   !-----------------------------------------------------------------

   ! Allocate FSOAP_CONV (gas-to-particle conversion)
   ALLOCATE( FSOAP_CONV( NX, NY, NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:FSOAP_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   FSOAP_CONV = 0e+0_fp

   ! Allocate FSOAS_CONV
   ALLOCATE( FSOAS_CONV( NX, NY, NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:FSOAS_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   FSOAS_CONV = 0e+0_fp

   ! Allocate BRCSOA_CONV
   ALLOCATE( BRCSOA_CONV( NX, NY, NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:BRCSOA_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   BRCSOA_CONV = 0e+0_fp

   ! Allocate NPBRC_CONV (NPBRCPOA -> WTC conversion)
   ALLOCATE( NPBRC_CONV( NX, NY, NZ ), STAT=RC )
   CALL GC_CheckVar( 'brc_mod.F90:NPBRC_CONV', 0, RC )
   IF ( RC /= GC_SUCCESS ) RETURN
   NPBRC_CONV = 0e+0_fp

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
!  24 Feb 2026 - M. Harvey - Initial version
!  25 Feb 2026 - M. Harvey - Added FSOAP_CONV and diagnostic arrays
!EOP
!------------------------------------------------------------------------------
!BOC

   !=================================================================
   ! Cleanup_BrC begins here!
   !=================================================================
   RC = GC_SUCCESS

   ! Conversion arrays
   IF ( ALLOCATED( FSOAP_CONV ) ) THEN
      DEALLOCATE( FSOAP_CONV, STAT=RC )
      CALL GC_CheckVar( 'brc_mod.F90:FSOAP_CONV', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
   ENDIF

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

   IF ( ALLOCATED( NPBRC_CONV ) ) THEN
      DEALLOCATE( NPBRC_CONV, STAT=RC )
      CALL GC_CheckVar( 'brc_mod.F90:NPBRC_CONV', 2, RC )
      IF ( RC /= GC_SUCCESS ) RETURN
   ENDIF

 END SUBROUTINE Cleanup_BrC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_fsoap
!
! !DESCRIPTION: Subroutine CHEM\_FSOAP converts the gas-phase fire SOA
!  precursor FSOAP to particle-phase FSOAS via a first-order process
!  with e-folding time FSOAP\_LIFE days.
!
!  This mimics the standard GEOS-Chem treatment of fire SOA precursor
!  gases (Pai et al. 2020) where gaseous SVOCs condense to form SOA
!  with a ~1 day timescale.
!
!  The converted mass is stored in FSOAP\_CONV for uptake by CHEM\_FSOAS.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_FSOAP( Input_Opt,  State_Chm, State_Diag, &
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
   INTEGER,        INTENT(IN)    :: spcId        ! FSOAP species Id
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
!  Drydep is applied in mixing_mod.F90 (gas-phase dry dep for FSOAP).
!  The 1-day lifetime is consistent with Pai et al. (2020) fire SOA
!  precursor aging timescale used in standard GEOS-Chem.
!
! !REVISION HISTORY:
!  25 Feb 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I,      J,   L
   REAL(fp)            :: DTCHEM, KFSOAP, FREQ, TC0, CNEW, RKT

   ! Pointers
   REAL(fp), POINTER   :: TC(:,:,:)
!
! !DEFINED PARAMETERS:
!
   ! E-folding lifetime for FSOAP gas-to-particle conversion [days]
   !   Consistent with Pai et al. (2020) fire SOA precursor aging
   REAL(fp), PARAMETER :: FSOAP_LIFE = 1.0e+0_fp

   !=================================================================
   ! CHEM_FSOAP begins here!
   !=================================================================

   ! Assume success
   RC         = GC_SUCCESS

   ! Initialize
   KFSOAP     = 1.e+0_fp / ( 86400e+0_fp * FSOAP_LIFE )
   DTCHEM     = GET_TS_CHEM()
   FSOAP_CONV = 0e+0_fp
   TC         => State_Chm%Species(spcId)%Conc

   !=================================================================
   ! Gas-to-particle conversion from FSOAP to FSOAS:
   !   First-order loss with e-folding time FSOAP_LIFE days
   !=================================================================
   !$OMP PARALLEL DO                                                &
   !$OMP DEFAULT( SHARED                                           )&
   !$OMP PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW                   )&
   !$OMP COLLAPSE( 3                                               )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial FSOAP mass [kg]
      TC0  = TC(I,J,L)

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Amount of FSOAP left after chemistry [kg]
      RKT  = ( KFSOAP + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of FSOAP converted to FSOAS [kg/timestep]
      FSOAP_CONV(I,J,L) = ( TC0 - CNEW )                          &
                         * KFSOAP / ( KFSOAP + FREQ )

      ! Store diagnostic: FSOAP->FSOAS flux [kg/timestep]
      IF ( State_Diag%Archive_BrCFluxFSOAP2FSOAS ) THEN
         State_Diag%BrCFluxFSOAP2FSOAS(I,J,L) = FSOAP_CONV(I,J,L)
      ENDIF 

      ! Store new concentration back into species array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_FSOAP
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
!  It also receives any condensed mass from the gas-phase precursor
!  FSOAP (stored in FSOAP\_CONV).
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
!  19 Feb 2026 - M. Harvey - Initial version
!  24 Feb 2026 - M. Harvey - Moved to brc_mod.F90
!  25 Feb 2026 - M. Harvey - Added FSOAP_CONV intake and diagnostics
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
   !   Updated from 0.25 d to 1.0 d based on Wong et al. (2019)
   !   and Hems et al. (2021) review of BrC darkening timescales
   REAL(fp), PARAMETER :: FSOAS_LIFE = 1.0e+0_fp

   !=================================================================
   ! CHEM_FSOAS begins here!
   !=================================================================

   ! Assume success
   RC         = GC_SUCCESS

   ! Initialize
   KFSOAS     = 1.e+0_fp / ( 86400e+0_fp * FSOAS_LIFE )
   DTCHEM     = GET_TS_CHEM()
   FSOAS_CONV = 0e+0_fp
   TC         => State_Chm%Species(spcId)%Conc

   !=================================================================
   ! Conversion from FSOAS to BRCSOA:
   !   First-order loss with e-folding time FSOAS_LIFE days
   !   Also receive condensed mass from FSOAP (if any)
   !=================================================================
   !$OMP PARALLEL DO                                                &
   !$OMP DEFAULT( SHARED                                           )&
   !$OMP PRIVATE( I, J, L, TC0, FREQ, RKT, CNEW                   )&
   !$OMP COLLAPSE( 3                                               )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Initial FSOAS mass [kg] + any condensed mass from FSOAP
      TC0  = TC(I,J,L) + FSOAP_CONV(I,J,L)

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Amount of FSOAS left after chemistry [kg]
      RKT  = ( KFSOAS + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount of FSOAS converted to BRCSOA [kg/timestep]
      FSOAS_CONV(I,J,L) = ( TC0 - CNEW )                          &
                         * KFSOAS / ( KFSOAS + FREQ )

      ! Store diagnostic: FSOAS->BRCSOA flux [kg/timestep]
      IF ( State_Diag%Archive_BrCFluxFSOAS2BRC ) THEN
         State_Diag%BrCFluxFSOAS2BRC(I,J,L) = FSOAS_CONV(I,J,L)
      ENDIF 

      ! Store new concentration back into species array
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! Zero FSOAP_CONV -- we have consumed it
   !=================================================================
   FSOAP_CONV = 0e+0_fp

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
!  (2) Photo-bleaches BRCSOA to WTC using a first-order loss whose
!      rate constant varies with local temperature and relative humidity,
!      following the viscosity-dependent parameterisation of
!      Schnitzler et al. (2022), PNAS 119(38), e2205610119.
!
!  The bleached mass is stored in BRCSOA\_CONV for uptake by CHEM\_WTC.
!
!  Physics summary:
!    T, RH -> viscosity (VFT + Arrhenius mixing, SI Eq. S1-S6)
!          -> D_O3     (fractional Stokes-Einstein, SI Eq. S12-S13)
!          -> tau_BrC  (resistor-model lifetime, Main Eq. 3)
!          -> k_bleach = 1/tau_BrC  [s^-1]
!
!  This gives tau ~ hours at surface, >> 1 week in upper troposphere,
!  reproducing Fig. 3A (1-day contour at ~1 km altitude).
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_BRCSOA( Input_Opt,  State_Chm, State_Diag,  &
                         State_Grid, State_Met, spcId, RC     )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Chm_Mod,  ONLY : Ind_
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Met_Mod,  ONLY : MetState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
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
!  The bleaching lifetime is computed locally at each grid cell from
!  temperature and RH via the Schnitzler et al. (2022) parameterisation.
!
! !REVISION HISTORY:
!  19 Feb 2026 - M. Harvey - Initial version (fixed-rate bleaching)
!  24 Feb 2026 - M. Harvey - Moved to brc_mod.F90; replaced fixed rate
!                             with T/RH-dependent Schnitzler et al. (2022)
!                             viscosity parameterisation
!  25 Feb 2026 - M. Harvey - Added development diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I, J, L
   INTEGER             :: id_O3          ! O3 species index
   REAL(fp)            :: DTCHEM
   REAL(fp)            :: KBRCSOA_LOCAL  ! Local bleaching rate [s^-1]
   REAL(fp)            :: FREQ           ! Drydep freq (zero here)
   REAL(fp)            :: TC0, CNEW, RKT, CCV
   REAL(fp)            :: T_LOCAL        ! Local temperature [K]
   REAL(fp)            :: AW_LOCAL       ! Local water activity [0-1]
   REAL(fp)            :: TAU_LOCAL      ! Local bleaching lifetime [s]
   REAL(fp)            :: ETA_LOCAL      ! Local BBOA viscosity [Pa s]
   REAL(fp)            :: P_O3_LOCAL     ! Local O3 partial pressure [atm]
   REAL(fp)            :: X_O3           ! Local O3 mixing ratio [mol/mol]
   REAL(fp)            :: ALT_TOP        ! Height of box top AGL [m]

   ! Pointers
   REAL(fp), POINTER   :: TC(:,:,:)

   !=================================================================
   ! CHEM_BRCSOA begins here!
   !=================================================================

   ! Assume success
   RC          = GC_SUCCESS

   ! Chemistry timestep [s]
   DTCHEM      = GET_TS_CHEM()

   ! IMPORTANT: Do NOT zero FSOAS_CONV here -- it was set in CHEM_FSOAS
   !            and we consume it below, then zero it at the end.
   BRCSOA_CONV = 0e+0_fp

   TC          => State_Chm%Species(spcId)%Conc

   ! Look up O3 species index for local partial pressure
   ! Falls back to global constant P_O3_ATM (35 ppb) if O3 not found
   id_O3 = Ind_('O3')

   !=================================================================
   ! Photo-bleaching from BRCSOA to WTC
   ! Rate constant varies with local T and RH following the
   ! viscosity-dependent parameterisation (Schnitzler et al. 2022)
   !=================================================================
   !$OMP PARALLEL DO                                                &
   !$OMP DEFAULT( SHARED                                           )&
   !$OMP PRIVATE( I, J, L, T_LOCAL, AW_LOCAL, TAU_LOCAL            )&
   !$OMP PRIVATE( ETA_LOCAL, P_O3_LOCAL, X_O3, ALT_TOP             )&
   !$OMP PRIVATE( KBRCSOA_LOCAL, CCV, TC0, FREQ, RKT, CNEW        )&
   !$OMP COLLAPSE( 3                                               )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      !==============================================================
      ! Compute local bleaching rate
      ! Scheme is selected by runtime BrC_Bleach_Scheme
      !==============================================================

      ! Local temperature [K] -- always needed for diagnostics
      T_LOCAL = State_Met%T(I,J,L)

      ! Local water activity [0-1] -- always needed for diagnostics
      AW_LOCAL = State_Met%RH(I,J,L) / 100.0_fp
      AW_LOCAL = MAX( AW_LOCAL, 0.0_fp  )
      AW_LOCAL = MIN( AW_LOCAL, 0.99_fp )

      ! Determine bleaching lifetime based on selected scheme
      SELECT CASE ( Input_Opt%BrC_Bleach_Scheme )

      CASE ( 0 )
         !--- No bleaching ---
         TAU_LOCAL = TAU_MAX

      CASE ( 1 )
         !--- Fixed 1-day lifetime everywhere ---
         TAU_LOCAL = TAU_1DAY

      CASE ( 2 )
         !--- Fixed 1-day below 1 km AGL, no bleaching above ---
         ALT_TOP = SUM( State_Met%BXHEIGHT(I,J,1:L) )
         IF ( ALT_TOP <= ALT_THRESH ) THEN
            TAU_LOCAL = TAU_1DAY
         ELSE
            TAU_LOCAL = TAU_MAX
         ENDIF

      CASE ( 3 )
         !--- Viscosity-dependent, fixed 35 ppb O3 ---
         TAU_LOCAL = CALC_TAU_BRC( T_LOCAL, AW_LOCAL, P_O3_ATM )

      CASE ( 4 )
         !--- Viscosity-dependent, local O3 (default) ---
         IF ( id_O3 > 0 ) THEN
            X_O3 = ( State_Chm%Species(id_O3)%Conc(I,J,L) * MW_AIR ) &
                 / ( State_Met%AD(I,J,L) * MW_O3 )
            X_O3 = MAX( X_O3, 0.0_fp )
            P_O3_LOCAL = X_O3 * State_Met%PMID(I,J,L) * HPA2ATM
            P_O3_LOCAL = MAX( P_O3_LOCAL, 1.0e-12_fp )
         ELSE
            P_O3_LOCAL = P_O3_ATM
         ENDIF
         TAU_LOCAL = CALC_TAU_BRC( T_LOCAL, AW_LOCAL, P_O3_LOCAL )

      CASE DEFAULT
         TAU_LOCAL = TAU_MAX

      END SELECT

      ! First-order rate constant [s^-1]
      KBRCSOA_LOCAL = 1.0_fp / TAU_LOCAL

      ! Local BBOA viscosity [Pa s] for diagnostics
      ETA_LOCAL = VISC_BBOA_FUNC( T_LOCAL, AW_LOCAL )

      !==============================================================
      ! Store development diagnostics
      !==============================================================
      IF ( State_Diag%Archive_BrCTauBleach ) THEN
         State_Diag%BrCTauBleach(I,J,L) = TAU_LOCAL
      ENDIF

      IF ( State_Diag%Archive_BrCKBleach ) THEN
         State_Diag%BrCKBleach(I,J,L) = KBRCSOA_LOCAL
      ENDIF

      IF ( State_Diag%Archive_BrCEtaBBOA ) THEN
         State_Diag%BrCEtaBBOA(I,J,L) = ETA_LOCAL
      ENDIF      

      !==============================================================
      ! 1) Add newly formed BRCSOA from FSOAS (darkening step)
      !==============================================================
      ! Convert FSOAS from OM mass to carbon mass.  BRCSOA aerosol mass
      ! later uses the model's oxidized OC OM:OC, representing additional
      ! oxygenated mass gained during darkening/aging.
      CCV = FSOAS_CONV(I,J,L) / OMOC_BBOA

      ! BRCSOA mass available to bleach this timestep [kg]
      TC0 = TC(I,J,L) + CCV

      !==============================================================
      ! 2) Bleach BRCSOA -> WTC as first-order loss over DTCHEM
      !==============================================================

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Remaining BRCSOA after bleaching [kgC].  No timestep floor is
      ! applied; persistent primary BrC is represented by PBRCPOA.
      RKT  = ( KBRCSOA_LOCAL + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount bleached from BRCSOA to WTC [kg/timestep]
      BRCSOA_CONV(I,J,L) = ( TC0 - CNEW )                         &
                          * KBRCSOA_LOCAL / ( KBRCSOA_LOCAL + FREQ )

      ! Store diagnostic: BRCSOA->WTC flux [kg/timestep]
      IF ( State_Diag%Archive_BrCFluxBRC2WTC ) THEN
         State_Diag%BrCFluxBRC2WTC(I,J,L) = BRCSOA_CONV(I,J,L)
      ENDIF 

      ! Store updated BRCSOA back into species array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   !=================================================================
   ! We have now consumed FSOAS_CONV for this timestep -- zero it
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
! !IROUTINE: chem_npbrcpoa
!
! !DESCRIPTION: Subroutine CHEM\_NPBRCPOA photo-bleaches non-persistent
!  BrC-POA to WTC using the same viscosity-dependent parameterisation of
!  Schnitzler et al. (2022) as used for BRCSOA.
!
!  NPBRCPOA represents the fraction of fire primary organic aerosol that
!  contains low-MW BrC chromophores (e.g. methoxyphenols, nitrophenols,
!  lignin fragments) that are susceptible to O3 bleaching.  Unlike
!  DBRCPOA (persistent/dark BrC), NPBRCPOA bleaches on timescales
!  governed by local environmental conditions.
!
!  The bleached mass is stored in NPBRC\_CONV for uptake by CHEM\_WTC.
!\\ 
!\\
! !INTERFACE:
!
 SUBROUTINE CHEM_NPBRCPOA( Input_Opt,  State_Chm, State_Diag,  &
                           State_Grid, State_Met, spcId, RC     )
!
! !USES:
!
   USE ErrCode_Mod
   USE Input_Opt_Mod,  ONLY : OptInput
   USE State_Chm_Mod,  ONLY : ChmState
   USE State_Chm_Mod,  ONLY : Ind_
   USE State_Diag_Mod, ONLY : DgnState
   USE State_Grid_Mod, ONLY : GrdState
   USE State_Met_Mod,  ONLY : MetState
   USE TIME_MOD,       ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
   TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
   INTEGER,        INTENT(IN)    :: spcId        ! NPBRCPOA species Id
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
!  The bleaching lifetime uses the same Schnitzler et al. (2022)
!  viscosity-dependent parameterisation as CHEM_BRCSOA.
!  NPBRCPOA is directly emitted from fires (no darkening precursor).
!
! !REVISION HISTORY:
!  25 Feb 2026 - M. Harvey - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER             :: I, J, L
   INTEGER             :: id_O3          ! O3 species index
   REAL(fp)            :: DTCHEM
   REAL(fp)            :: KNPBRC_LOCAL   ! Local bleaching rate [s^-1]
   REAL(fp)            :: FREQ           ! Drydep freq (zero here)
   REAL(fp)            :: TC0, CNEW, RKT
   REAL(fp)            :: T_LOCAL        ! Local temperature [K]
   REAL(fp)            :: AW_LOCAL       ! Local water activity [0-1]
   REAL(fp)            :: TAU_LOCAL      ! Local bleaching lifetime [s]
   REAL(fp)            :: P_O3_LOCAL     ! Local O3 partial pressure [atm]
   REAL(fp)            :: X_O3           ! Local O3 mixing ratio [mol/mol]
   REAL(fp)            :: ALT_TOP        ! Height of box top AGL [m]

   ! Pointers
   REAL(fp), POINTER   :: TC(:,:,:)

   !=================================================================
   ! CHEM_NPBRCPOA begins here!
   !=================================================================

   ! Assume success
   RC          = GC_SUCCESS

   ! Chemistry timestep [s]
   DTCHEM      = GET_TS_CHEM()

   ! Zero the conversion array for this timestep
   NPBRC_CONV  = 0e+0_fp

   TC          => State_Chm%Species(spcId)%Conc

   ! Look up O3 species index for local partial pressure
   ! Falls back to global constant P_O3_ATM (35 ppb) if O3 not found
   id_O3 = Ind_('O3')

   !=================================================================
   ! Photo-bleaching from NPBRCPOA to WTC
   ! Uses the same viscosity-dependent Schnitzler et al. (2022)
   ! parameterisation as BRCSOA -> WTC
   !=================================================================
   !$OMP PARALLEL DO                                                &
   !$OMP DEFAULT( SHARED                                           )&
   !$OMP PRIVATE( I, J, L, T_LOCAL, AW_LOCAL, TAU_LOCAL            )&
   !$OMP PRIVATE( P_O3_LOCAL, X_O3, ALT_TOP                        )&
   !$OMP PRIVATE( KNPBRC_LOCAL, TC0, FREQ, RKT, CNEW              )&
   !$OMP COLLAPSE( 3                                               )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      !==============================================================
      ! Compute local bleaching rate
      ! Scheme is selected by runtime BrC_Bleach_Scheme
      ! (identical scheme selection to CHEM_BRCSOA)
      !==============================================================

      ! Local temperature [K]
      T_LOCAL = State_Met%T(I,J,L)

      ! Local water activity [0-1]
      AW_LOCAL = State_Met%RH(I,J,L) / 100.0_fp
      AW_LOCAL = MAX( AW_LOCAL, 0.0_fp  )
      AW_LOCAL = MIN( AW_LOCAL, 0.99_fp )

      ! Determine bleaching lifetime based on selected scheme
      SELECT CASE ( Input_Opt%BrC_Bleach_Scheme )

      CASE ( 0 )
         !--- No bleaching ---
         TAU_LOCAL = TAU_MAX

      CASE ( 1 )
         !--- Fixed 1-day lifetime everywhere ---
         TAU_LOCAL = TAU_1DAY

      CASE ( 2 )
         !--- Fixed 1-day below 1 km AGL, no bleaching above ---
         ALT_TOP = SUM( State_Met%BXHEIGHT(I,J,1:L) )
         IF ( ALT_TOP <= ALT_THRESH ) THEN
            TAU_LOCAL = TAU_1DAY
         ELSE
            TAU_LOCAL = TAU_MAX
         ENDIF

      CASE ( 3 )
         !--- Viscosity-dependent, fixed 35 ppb O3 ---
         TAU_LOCAL = CALC_TAU_BRC( T_LOCAL, AW_LOCAL, P_O3_ATM )

      CASE ( 4 )
         !--- Viscosity-dependent, local O3 (default) ---
         IF ( id_O3 > 0 ) THEN
            X_O3 = ( State_Chm%Species(id_O3)%Conc(I,J,L) * MW_AIR ) &
                 / ( State_Met%AD(I,J,L) * MW_O3 )
            X_O3 = MAX( X_O3, 0.0_fp )
            P_O3_LOCAL = X_O3 * State_Met%PMID(I,J,L) * HPA2ATM
            P_O3_LOCAL = MAX( P_O3_LOCAL, 1.0e-12_fp )
         ELSE
            P_O3_LOCAL = P_O3_ATM
         ENDIF
         TAU_LOCAL = CALC_TAU_BRC( T_LOCAL, AW_LOCAL, P_O3_LOCAL )

      CASE DEFAULT
         TAU_LOCAL = TAU_MAX

      END SELECT

      ! First-order rate constant [s^-1]
      KNPBRC_LOCAL = 1.0_fp / TAU_LOCAL

      !==============================================================
      ! Bleach NPBRCPOA -> WTC as first-order loss over DTCHEM
      !==============================================================

      ! Current NPBRCPOA mass [kg]
      TC0 = TC(I,J,L)

      ! Zero drydep freq (drydep handled in mixing_mod.F90)
      FREQ = 0e+0_fp

      ! Remaining NPBRCPOA after bleaching [kgC].  The persistent share
      ! of primary BrC emissions is emitted separately as PBRCPOA.
      RKT  = ( KNPBRC_LOCAL + FREQ ) * DTCHEM
      CNEW = TC0 * EXP( -RKT )

      ! Prevent underflow condition
      IF ( CNEW < SMALLNUM ) CNEW = 0e+0_fp

      ! Amount bleached from NPBRCPOA to WTC [kg/timestep]
      NPBRC_CONV(I,J,L) = ( TC0 - CNEW )                          &
                         * KNPBRC_LOCAL / ( KNPBRC_LOCAL + FREQ )

      ! Store diagnostic: NPBRCPOA->WTC flux [kg/timestep]
      IF ( State_Diag%Archive_BrCFluxNPBRC2WTC ) THEN
         State_Diag%BrCFluxNPBRC2WTC(I,J,L) = NPBRC_CONV(I,J,L)
      ENDIF 

      ! Store updated NPBRCPOA back into species array [kg]
      TC(I,J,L) = CNEW

   ENDDO
   ENDDO
   ENDDO
   !$OMP END PARALLEL DO

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_NPBRCPOA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_wtc
!
! !DESCRIPTION: Subroutine CHEM\_WTC receives the bleached mass from
!  BRCSOA (stored in BRCSOA\_CONV) and from NPBRCPOA (stored in
!  NPBRC\_CONV, if present) and adds them to the WTC tracer.
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
!  19 Feb 2026 - M. Harvey - Initial version
!  24 Feb 2026 - M. Harvey - Moved to brc_mod.F90
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

   !$OMP PARALLEL DO                                                &
   !$OMP DEFAULT( SHARED                                           )&
   !$OMP PRIVATE( I, J, L, TC0, CCV, CNEW                         )&
   !$OMP COLLAPSE( 3                                               )
   DO L = 1, State_Grid%NZ
   DO J = 1, State_Grid%NY
   DO I = 1, State_Grid%NX

      ! Current WTC mass [kg]
      TC0 = TC(I,J,L)

      ! Bleached mass arriving from BRCSOA [kg]
      CCV = BRCSOA_CONV(I,J,L)

      ! Also add bleached mass arriving from NPBRCPOA [kg]
      CCV = CCV + NPBRC_CONV(I,J,L)

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
   ! Zero conversion arrays for next timestep
   !=================================================================
   BRCSOA_CONV = 0e+0_fp
   NPBRC_CONV  = 0e+0_fp

   ! Free pointer
   TC => NULL()

 END SUBROUTINE CHEM_WTC
!EOC
END MODULE BRC_MOD
