!------------------------------------------------------------------------------
! Pure arithmetic for the default-off exact-band dry absorption diagnostic.
! This module does not read or mutate GEOS-Chem state.
!------------------------------------------------------------------------------
MODULE Dry_Absorption_Diag_Mod

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY : IEEE_IS_FINITE
  USE Precision_Mod, ONLY : fp

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: DRYABS_OK              = 0
  INTEGER, PARAMETER, PUBLIC :: DRYABS_INVALID_INPUT   = 1
  INTEGER, PARAMETER, PUBLIC :: DRYABS_INVALID_BRACKET = 2
  INTEGER, PARAMETER, PUBLIC :: DRYABS_MIXED_ZERO      = 3

  PUBLIC :: Dry_Absorption_Optical_Depth
  PUBLIC :: Interpolate_Dry_Absorption
  PUBLIC :: Layer_Absorption_Mm1
  PUBLIC :: Ambient_To_Reference_Mm1
  PUBLIC :: Absorption_Angstrom_Exponent
  PUBLIC :: Check_Component_Closure
  PUBLIC :: Compute_Dry_Absorption_Diagnostics

CONTAINS

  PURE SUBROUTINE Dry_Absorption_Optical_Depth( Extinction_OD, SSA, &
                                                Absorption_OD, Status )

    REAL(fp), INTENT(IN)  :: Extinction_OD
    REAL(fp), INTENT(IN)  :: SSA
    REAL(fp), INTENT(OUT) :: Absorption_OD
    INTEGER,  INTENT(OUT) :: Status

    Absorption_OD = 0.0_fp
    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Extinction_OD ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( SSA ) ) RETURN
    IF ( Extinction_OD < 0.0_fp ) RETURN
    IF ( SSA < 0.0_fp .or. SSA > 1.0_fp ) RETURN

    Absorption_OD = Extinction_OD * ( 1.0_fp - SSA )
    Status = DRYABS_OK

  END SUBROUTINE Dry_Absorption_Optical_Depth

  PURE SUBROUTINE Interpolate_Dry_Absorption( Target_Wavelength_Nm, &
                                               Lower_Wavelength_Nm,  &
                                               Upper_Wavelength_Nm,  &
                                               Lower_Absorption_OD,  &
                                               Upper_Absorption_OD,  &
                                               Absorption_OD, Status )

    REAL(fp), INTENT(IN)  :: Target_Wavelength_Nm
    REAL(fp), INTENT(IN)  :: Lower_Wavelength_Nm
    REAL(fp), INTENT(IN)  :: Upper_Wavelength_Nm
    REAL(fp), INTENT(IN)  :: Lower_Absorption_OD
    REAL(fp), INTENT(IN)  :: Upper_Absorption_OD
    REAL(fp), INTENT(OUT) :: Absorption_OD
    INTEGER,  INTENT(OUT) :: Status
    REAL(fp)              :: Weight

    Absorption_OD = 0.0_fp
    Status = DRYABS_INVALID_BRACKET
    IF ( .not. IEEE_IS_FINITE( Target_Wavelength_Nm ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Lower_Wavelength_Nm ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Upper_Wavelength_Nm ) ) RETURN
    IF ( Lower_Wavelength_Nm <= 0.0_fp ) RETURN
    IF ( Upper_Wavelength_Nm <= Lower_Wavelength_Nm ) RETURN
    IF ( Target_Wavelength_Nm < Lower_Wavelength_Nm ) RETURN
    IF ( Target_Wavelength_Nm > Upper_Wavelength_Nm ) RETURN

    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Lower_Absorption_OD ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Upper_Absorption_OD ) ) RETURN
    IF ( Lower_Absorption_OD < 0.0_fp ) RETURN
    IF ( Upper_Absorption_OD < 0.0_fp ) RETURN

    IF ( Target_Wavelength_Nm == Lower_Wavelength_Nm ) THEN
       Absorption_OD = Lower_Absorption_OD
       Status = DRYABS_OK
       RETURN
    ENDIF
    IF ( Target_Wavelength_Nm == Upper_Wavelength_Nm ) THEN
       Absorption_OD = Upper_Absorption_OD
       Status = DRYABS_OK
       RETURN
    ENDIF

    IF ( Lower_Absorption_OD == 0.0_fp .and. &
         Upper_Absorption_OD == 0.0_fp ) THEN
       Status = DRYABS_OK
       RETURN
    ENDIF
    IF ( Lower_Absorption_OD == 0.0_fp .or. &
         Upper_Absorption_OD == 0.0_fp ) THEN
       Status = DRYABS_MIXED_ZERO
       RETURN
    ENDIF

    Weight = LOG( Target_Wavelength_Nm / Lower_Wavelength_Nm ) / &
             LOG( Upper_Wavelength_Nm / Lower_Wavelength_Nm )
    Absorption_OD = EXP( LOG( Lower_Absorption_OD ) + Weight * &
                    ( LOG( Upper_Absorption_OD ) - &
                      LOG( Lower_Absorption_OD ) ) )
    Status = DRYABS_OK

  END SUBROUTINE Interpolate_Dry_Absorption

  PURE SUBROUTINE Layer_Absorption_Mm1( Absorption_OD, Box_Height_M, &
                                        Absorption_Mm1, Status )

    REAL(fp), INTENT(IN)  :: Absorption_OD
    REAL(fp), INTENT(IN)  :: Box_Height_M
    REAL(fp), INTENT(OUT) :: Absorption_Mm1
    INTEGER,  INTENT(OUT) :: Status

    Absorption_Mm1 = 0.0_fp
    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Absorption_OD ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Box_Height_M ) ) RETURN
    IF ( Absorption_OD < 0.0_fp ) RETURN
    IF ( Box_Height_M <= 0.0_fp ) RETURN

    Absorption_Mm1 = Absorption_OD / Box_Height_M * 1.0e6_fp
    Status = DRYABS_OK

  END SUBROUTINE Layer_Absorption_Mm1

  PURE SUBROUTINE Ambient_To_Reference_Mm1( Ambient_Absorption_Mm1, &
                                             Ambient_Pressure_Hpa,   &
                                             Ambient_Temperature_K,  &
                                             Reference_Pressure_Hpa, &
                                             Reference_Temperature_K,&
                                             Reference_Absorption_Mm1,&
                                             Status )

    REAL(fp), INTENT(IN)  :: Ambient_Absorption_Mm1
    REAL(fp), INTENT(IN)  :: Ambient_Pressure_Hpa
    REAL(fp), INTENT(IN)  :: Ambient_Temperature_K
    REAL(fp), INTENT(IN)  :: Reference_Pressure_Hpa
    REAL(fp), INTENT(IN)  :: Reference_Temperature_K
    REAL(fp), INTENT(OUT) :: Reference_Absorption_Mm1
    INTEGER,  INTENT(OUT) :: Status

    Reference_Absorption_Mm1 = 0.0_fp
    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Ambient_Absorption_Mm1 ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Ambient_Pressure_Hpa ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Ambient_Temperature_K ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Reference_Pressure_Hpa ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Reference_Temperature_K ) ) RETURN
    IF ( Ambient_Absorption_Mm1 < 0.0_fp ) RETURN
    IF ( Ambient_Pressure_Hpa <= 0.0_fp ) RETURN
    IF ( Ambient_Temperature_K <= 0.0_fp ) RETURN
    IF ( Reference_Pressure_Hpa <= 0.0_fp ) RETURN
    IF ( Reference_Temperature_K <= 0.0_fp ) RETURN

    Reference_Absorption_Mm1 = Ambient_Absorption_Mm1 * &
         ( Reference_Pressure_Hpa * Ambient_Temperature_K ) / &
         ( Ambient_Pressure_Hpa * Reference_Temperature_K )
    Status = DRYABS_OK

  END SUBROUTINE Ambient_To_Reference_Mm1

  PURE SUBROUTINE Absorption_Angstrom_Exponent( Short_Absorption_Mm1, &
                                                 Long_Absorption_Mm1,  &
                                                 Short_Wavelength_Nm,  &
                                                 Long_Wavelength_Nm,   &
                                                 AAE, Status )

    REAL(fp), INTENT(IN)  :: Short_Absorption_Mm1
    REAL(fp), INTENT(IN)  :: Long_Absorption_Mm1
    REAL(fp), INTENT(IN)  :: Short_Wavelength_Nm
    REAL(fp), INTENT(IN)  :: Long_Wavelength_Nm
    REAL(fp), INTENT(OUT) :: AAE
    INTEGER,  INTENT(OUT) :: Status

    AAE = 0.0_fp
    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Short_Absorption_Mm1 ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Long_Absorption_Mm1 ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Short_Wavelength_Nm ) ) RETURN
    IF ( .not. IEEE_IS_FINITE( Long_Wavelength_Nm ) ) RETURN
    IF ( Short_Absorption_Mm1 <= 0.0_fp ) RETURN
    IF ( Long_Absorption_Mm1 <= 0.0_fp ) RETURN
    IF ( Short_Wavelength_Nm <= 0.0_fp ) RETURN
    IF ( Long_Wavelength_Nm <= Short_Wavelength_Nm ) RETURN

    AAE = -LOG( Short_Absorption_Mm1 / Long_Absorption_Mm1 ) / &
           LOG( Short_Wavelength_Nm / Long_Wavelength_Nm )
    Status = DRYABS_OK

  END SUBROUTINE Absorption_Angstrom_Exponent

  PURE SUBROUTINE Check_Component_Closure( Total_Absorption_Mm1, &
                                            Component_Absorption_Mm1, &
                                            Relative_Tolerance,       &
                                            Absolute_Tolerance_Mm1,   &
                                            Status )

    REAL(fp), INTENT(IN)  :: Total_Absorption_Mm1
    REAL(fp), INTENT(IN)  :: Component_Absorption_Mm1(:)
    REAL(fp), INTENT(IN)  :: Relative_Tolerance
    REAL(fp), INTENT(IN)  :: Absolute_Tolerance_Mm1
    INTEGER,  INTENT(OUT) :: Status
    REAL(fp)              :: Component_Total, Tolerance

    Status = DRYABS_INVALID_INPUT
    IF ( .not. IEEE_IS_FINITE( Total_Absorption_Mm1 ) ) RETURN
    IF ( .not. ALL( IEEE_IS_FINITE( Component_Absorption_Mm1 ) ) ) RETURN
    IF ( Total_Absorption_Mm1 < 0.0_fp ) RETURN
    IF ( ANY( Component_Absorption_Mm1 < 0.0_fp ) ) RETURN
    IF ( SIZE( Component_Absorption_Mm1 ) < 1 ) RETURN
    IF ( Relative_Tolerance < 0.0_fp ) RETURN
    IF ( Absolute_Tolerance_Mm1 < 0.0_fp ) RETURN

    Component_Total = SUM( Component_Absorption_Mm1 )
    Tolerance = Absolute_Tolerance_Mm1 + Relative_Tolerance * &
                ABS( Total_Absorption_Mm1 )
    IF ( ABS( Total_Absorption_Mm1 - Component_Total ) > Tolerance ) RETURN
    Status = DRYABS_OK

  END SUBROUTINE Check_Component_Closure

  SUBROUTINE Compute_Dry_Absorption_Diagnostics( Input_Opt, State_Chm, &
                                                  State_Diag, State_Grid, &
                                                  State_Met, RC )

    USE CMN_SIZE_Mod,    ONLY : NDUST, NRHAER, NSTRATAER
    USE ErrCode_Mod,     ONLY : GC_FAILURE, GC_SUCCESS
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState

    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(GrdState), INTENT(IN)    :: State_Grid
    TYPE(MetState), INTENT(IN)    :: State_Met
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
    INTEGER,        INTENT(OUT)   :: RC

    INTEGER  :: I, J, L, N, NA, W, IWV1, IWV2, Status, Category
    INTEGER  :: Idst, Istrat, SpcId
    INTEGER  :: Optical_Bin(NRHAER)
    REAL(fp) :: Density(NRHAER)
    REAL(fp) :: Components(5), Total
    REAL(fp) :: Ext1, Ext2, Abs1, Abs2, AbsTarget
    REAL(fp) :: Radius1, Radius2, Q1, Q2, Ssa1, Ssa2
    REAL(fp) :: Hydrophobic_Mass, Hydrophobic_Radius
    REAL(fp) :: Bc_Factor
    LOGICAL  :: Archive_Any
    CHARACTER(LEN=32) :: SpcName

    RC = GC_SUCCESS
    Archive_Any = State_Diag%Archive_AerosolDryAbsWL1 .or. &
                  State_Diag%Archive_AerosolDryAbsWL2 .or. &
                  State_Diag%Archive_AerosolDryAbsWL3
    IF ( .not. Archive_Any ) RETURN
    IF ( Input_Opt%NWVSELECT < 3 ) THEN
       RC = GC_FAILURE
       RETURN
    ENDIF

    CALL Clear_Dry_Absorption_Diagnostics( State_Diag )

    Optical_Bin = 0
    Density = 0.0_fp
    DO NA = 1, NRHAER
       SpcId = State_Chm%Map_HygGrth(NA)
       SpcName = TRIM( State_Chm%SpcData(SpcId)%Info%Name )
       SELECT CASE ( TRIM( SpcName ) )
          CASE ( 'SO4' )
             N = 1
          CASE ( 'BCPI' )
             N = 2
          CASE ( 'OCPI', 'POA1' )
             N = 3
          CASE ( 'SALA' )
             N = 4
          CASE ( 'SALC' )
             N = 5
          CASE ( 'BRCSOA' )
             N = 6
          CASE ( 'NPBRCPOA' )
             N = 7
          CASE ( 'WTC' )
             N = 8
          CASE ( 'FSOAS' )
             N = 9
          CASE ( 'PBRCPOA' )
             N = 10
          CASE ( 'DBRCPOA' )
             N = 11
          CASE DEFAULT
             RC = GC_FAILURE
             RETURN
       END SELECT
       Optical_Bin(NA) = N
       Density(N) = State_Chm%SpcData(SpcId)%Info%Density
    ENDDO
    IF ( ANY( Density <= 0.0_fp ) ) THEN
       RC = GC_FAILURE
       RETURN
    ENDIF

    Idst = NRHAER + NSTRATAER + 1
    DO W = 1, 3
       IWV1 = State_Chm%Phot%IWVSELECT(1,W)
       IWV2 = State_Chm%Phot%IWVSELECT(2,W)
       IF ( IWV1 < 1 .or. IWV2 < 1 ) THEN
          RC = GC_FAILURE
          RETURN
       ENDIF

       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE
          Components = 0.0_fp

          DO NA = 1, NRHAER
             N = Optical_Bin(NA)
             CALL Get_Dry_Lut_Optics( State_Chm, I, J, L, N, IWV1, &
                                      Radius1, Q1, Ssa1, Status )
             IF ( Status /= DRYABS_OK ) THEN
                RC = GC_FAILURE
                RETURN
             ENDIF
             CALL Get_Dry_Lut_Optics( State_Chm, I, J, L, N, IWV2, &
                                      Radius2, Q2, Ssa2, Status )
             IF ( Status /= DRYABS_OK ) THEN
                RC = GC_FAILURE
                RETURN
             ENDIF

             Ext1 = Dry_Extinction_From_Mass(                           &
                  State_Chm%AerMass%WAERSL(I,J,L,N), Density(N),        &
                  Radius1, Q1, State_Met%BXHEIGHT(I,J,L) )
             Ext2 = Dry_Extinction_From_Mass(                           &
                  State_Chm%AerMass%WAERSL(I,J,L,N), Density(N),        &
                  Radius2, Q2, State_Met%BXHEIGHT(I,J,L) )
             Bc_Factor = 1.0_fp
             IF ( N == 2 .and. Input_Opt%LBCAE ) Bc_Factor = Input_Opt%BCAE_1
             Abs1 = Ext1 * Bc_Factor * ( 1.0_fp - Ssa1 )
             Abs2 = Ext2 * Bc_Factor * ( 1.0_fp - Ssa2 )

             Hydrophobic_Mass = 0.0_fp
             Hydrophobic_Radius = Radius1
             IF ( N == 2 ) THEN
                Hydrophobic_Mass = State_Chm%AerMass%DAERSL(I,J,L,1)
             ELSE IF ( N == 3 ) THEN
                Hydrophobic_Mass = State_Chm%AerMass%DAERSL(I,J,L,2)
                Hydrophobic_Radius = State_Chm%AerMass%PDER(I,J,L)
             ELSE IF ( N == 11 ) THEN
                Hydrophobic_Mass = State_Chm%AerMass%DAERSL(I,J,L,3)
             ENDIF

             IF ( Hydrophobic_Mass > 0.0_fp ) THEN
                Bc_Factor = 1.0_fp
                IF ( N == 2 .and. Input_Opt%LBCAE ) Bc_Factor = Input_Opt%BCAE_2
                IF ( N == 3 ) THEN
                   IF ( Hydrophobic_Radius <= 0.0_fp ) THEN
                      RC = GC_FAILURE
                      RETURN
                   ENDIF
                   Ext1 = Dry_Extinction_From_Mass( Hydrophobic_Mass, &
                        Density(N), Hydrophobic_Radius, Q1,            &
                        State_Met%BXHEIGHT(I,J,L) )
                   Ext2 = Dry_Extinction_From_Mass( Hydrophobic_Mass, &
                        Density(N), Hydrophobic_Radius, Q2,            &
                        State_Met%BXHEIGHT(I,J,L) )
                ELSE
                   Ext1 = Dry_Extinction_From_Mass( Hydrophobic_Mass, &
                        Density(N), Radius1, Q1, State_Met%BXHEIGHT(I,J,L) )
                   Ext2 = Dry_Extinction_From_Mass( Hydrophobic_Mass, &
                        Density(N), Radius2, Q2, State_Met%BXHEIGHT(I,J,L) )
                ENDIF
                Abs1 = Abs1 + Ext1 * Bc_Factor * ( 1.0_fp - Ssa1 )
                Abs2 = Abs2 + Ext2 * Bc_Factor * ( 1.0_fp - Ssa2 )
             ENDIF

             CALL Interpolate_Requested_Absorption( Input_Opt%WVSELECT(W), &
                  IWV1, IWV2, State_Chm%Phot%WVAA(IWV1,1),                &
                  State_Chm%Phot%WVAA(IWV2,1), Abs1, Abs2, AbsTarget, Status )
             IF ( Status /= DRYABS_OK ) THEN
                RC = GC_FAILURE
                RETURN
             ENDIF
             Category = Dry_Absorption_Category( N )
             Components(Category) = Components(Category) + &
                  AbsTarget / State_Met%BXHEIGHT(I,J,L) * 1.0e6_fp
          ENDDO

          DO Istrat = 1, NSTRATAER
             N = NRHAER + Istrat
             Abs1 = State_Chm%Phot%ODAER(I,J,L,IWV1,N) * &
                    ( 1.0_fp - State_Chm%Phot%SSAA( &
                      IWV1,1,N,State_Chm%Phot%DRg) )
             Abs2 = State_Chm%Phot%ODAER(I,J,L,IWV2,N) * &
                    ( 1.0_fp - State_Chm%Phot%SSAA( &
                      IWV2,1,N,State_Chm%Phot%DRg) )
             CALL Interpolate_Requested_Absorption( Input_Opt%WVSELECT(W), &
                  IWV1, IWV2, State_Chm%Phot%WVAA(IWV1,1),                &
                  State_Chm%Phot%WVAA(IWV2,1), Abs1, Abs2, AbsTarget, Status )
             IF ( Status /= DRYABS_OK ) THEN
                RC = GC_FAILURE
                RETURN
             ENDIF
             Components(5) = Components(5) + &
                  AbsTarget / State_Met%BXHEIGHT(I,J,L) * 1.0e6_fp
          ENDDO

          IF ( Input_Opt%LDUST ) THEN
             DO N = 1, NDUST
                Abs1 = State_Chm%Phot%ODMDUST(I,J,L,IWV1,N) * &
                     ( 1.0_fp - State_Chm%Phot%SSAA( &
                       IWV1,N,Idst,State_Chm%Phot%DRg) )
                Abs2 = State_Chm%Phot%ODMDUST(I,J,L,IWV2,N) * &
                     ( 1.0_fp - State_Chm%Phot%SSAA( &
                       IWV2,N,Idst,State_Chm%Phot%DRg) )
                CALL Interpolate_Requested_Absorption( Input_Opt%WVSELECT(W), &
                     IWV1, IWV2, State_Chm%Phot%WVAA(IWV1,1),                &
                     State_Chm%Phot%WVAA(IWV2,1), Abs1, Abs2, AbsTarget, Status )
                IF ( Status /= DRYABS_OK ) THEN
                   RC = GC_FAILURE
                   RETURN
                ENDIF
                Components(4) = Components(4) + &
                     AbsTarget / State_Met%BXHEIGHT(I,J,L) * 1.0e6_fp
             ENDDO
          ENDIF

          Total = SUM( Components )
          CALL Store_Dry_Absorption( State_Diag, W, 1, I, J, L, Total )
          CALL Store_Dry_Absorption( State_Diag, W, 2, I, J, L, Components(1) )
          CALL Store_Dry_Absorption( State_Diag, W, 3, I, J, L, Components(2) )
          CALL Store_Dry_Absorption( State_Diag, W, 4, I, J, L, Components(3) )
          CALL Store_Dry_Absorption( State_Diag, W, 5, I, J, L, Components(4) )
          CALL Store_Dry_Absorption( State_Diag, W, 6, I, J, L, Components(5) )
       ENDDO
       ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE Compute_Dry_Absorption_Diagnostics

  PURE SUBROUTINE Interpolate_Requested_Absorption( Target_Wavelength_Nm, &
                                                     Lower_Index, Upper_Index, &
                                                     Lower_Wavelength_Nm,      &
                                                     Upper_Wavelength_Nm,      &
                                                     Lower_Absorption_OD,      &
                                                     Upper_Absorption_OD,      &
                                                     Absorption_OD, Status )

    REAL(fp), INTENT(IN)  :: Target_Wavelength_Nm
    INTEGER,  INTENT(IN)  :: Lower_Index, Upper_Index
    REAL(fp), INTENT(IN)  :: Lower_Wavelength_Nm, Upper_Wavelength_Nm
    REAL(fp), INTENT(IN)  :: Lower_Absorption_OD, Upper_Absorption_OD
    REAL(fp), INTENT(OUT) :: Absorption_OD
    INTEGER,  INTENT(OUT) :: Status

    IF ( Lower_Index == Upper_Index ) THEN
       Absorption_OD = Lower_Absorption_OD
       Status = DRYABS_OK
    ELSE
       CALL Interpolate_Dry_Absorption( Target_Wavelength_Nm, &
            Lower_Wavelength_Nm, Upper_Wavelength_Nm,         &
            Lower_Absorption_OD, Upper_Absorption_OD,         &
            Absorption_OD, Status )
    ENDIF

  END SUBROUTINE Interpolate_Requested_Absorption

  PURE FUNCTION Dry_Extinction_From_Mass( Mass_Kg_M3, Density_Kg_M3, &
                                           Radius_Um, Qext, Box_Height_M ) &
                                           RESULT( Extinction_OD )

    REAL(fp), INTENT(IN) :: Mass_Kg_M3, Density_Kg_M3, Radius_Um
    REAL(fp), INTENT(IN) :: Qext, Box_Height_M
    REAL(fp)             :: Extinction_OD

    Extinction_OD = 0.75_fp * Box_Height_M * Mass_Kg_M3 * Qext / &
                    ( Density_Kg_M3 * Radius_Um * 1.0e-6_fp )

  END FUNCTION Dry_Extinction_From_Mass

  PURE INTEGER FUNCTION Dry_Absorption_Category( Optical_Bin )

    INTEGER, INTENT(IN) :: Optical_Bin

    SELECT CASE ( Optical_Bin )
       CASE ( 2 )
          Dry_Absorption_Category = 1
       CASE ( 6, 7, 9, 10, 11 )
          Dry_Absorption_Category = 2
       CASE ( 3, 8 )
          Dry_Absorption_Category = 3
       CASE DEFAULT
          Dry_Absorption_Category = 5
    END SELECT

  END FUNCTION Dry_Absorption_Category

  SUBROUTINE Get_Dry_Lut_Optics( State_Chm, I, J, L, N, IWV, &
                                  Radius, Qext, SSA, Status )

    USE State_Chm_Mod, ONLY : ChmState

    TYPE(ChmState), INTENT(IN) :: State_Chm
    INTEGER, INTENT(IN)        :: I, J, L, N, IWV
    REAL(fp), INTENT(OUT)      :: Radius, Qext, SSA
    INTEGER, INTENT(OUT)       :: Status
    INTEGER                    :: G
    REAL(fp)                   :: Fraction

    Status = DRYABS_INVALID_INPUT
    IF ( N == 1 .or. N == 3 ) THEN
       G = 1
       DO WHILE ( State_Chm%AerMass%PDER(I,J,L) > &
                  State_Chm%Phot%REAA(1,N,G) .and. G < State_Chm%Phot%NDRg )
          G = G + 1
       ENDDO
       IF ( G == 1 ) THEN
          Radius = State_Chm%Phot%REAA(1,N,G)
          Qext   = State_Chm%Phot%QQAA(IWV,1,N,G)
          SSA    = State_Chm%Phot%SSAA(IWV,1,N,G)
       ELSE
          Fraction = ( State_Chm%AerMass%PDER(I,J,L) - &
                       State_Chm%Phot%REAA(1,N,G-1) ) / &
                     ( State_Chm%Phot%REAA(1,N,G) - &
                       State_Chm%Phot%REAA(1,N,G-1) )
          Fraction = MIN( Fraction, 1.0_fp )
          Radius = Fraction * State_Chm%Phot%REAA(1,N,G) + &
                   ( 1.0_fp - Fraction ) * State_Chm%Phot%REAA(1,N,G-1)
          Qext = Fraction * State_Chm%Phot%QQAA(IWV,1,N,G) + &
                 ( 1.0_fp - Fraction ) * State_Chm%Phot%QQAA(IWV,1,N,G-1)
          SSA = Fraction * State_Chm%Phot%SSAA(IWV,1,N,G) + &
                ( 1.0_fp - Fraction ) * State_Chm%Phot%SSAA(IWV,1,N,G-1)
       ENDIF
    ELSE
       Radius = State_Chm%Phot%REAA(1,N,State_Chm%Phot%DRg)
       Qext   = State_Chm%Phot%QQAA(IWV,1,N,State_Chm%Phot%DRg)
       SSA    = State_Chm%Phot%SSAA(IWV,1,N,State_Chm%Phot%DRg)
    ENDIF
    IF ( Radius <= 0.0_fp .or. Qext < 0.0_fp ) RETURN
    IF ( SSA < 0.0_fp .or. SSA > 1.0_fp ) RETURN
    Status = DRYABS_OK

  END SUBROUTINE Get_Dry_Lut_Optics

  SUBROUTINE Clear_Dry_Absorption_Diagnostics( State_Diag )

    USE State_Diag_Mod, ONLY : DgnState
    TYPE(DgnState), INTENT(INOUT) :: State_Diag

    IF ( State_Diag%Archive_AerosolDryAbsWL1 ) &
         State_Diag%AerosolDryAbsWL1 = 0.0_fp
    IF ( State_Diag%Archive_AerosolDryAbsWL2 ) &
         State_Diag%AerosolDryAbsWL2 = 0.0_fp
    IF ( State_Diag%Archive_AerosolDryAbsWL3 ) &
         State_Diag%AerosolDryAbsWL3 = 0.0_fp

  END SUBROUTINE Clear_Dry_Absorption_Diagnostics

  SUBROUTINE Store_Dry_Absorption( State_Diag, W, Component, I, J, L, Value )

    USE State_Diag_Mod, ONLY : DgnState
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
    INTEGER, INTENT(IN)           :: W, Component, I, J, L
    REAL(fp), INTENT(IN)          :: Value
    INTEGER                       :: Slot

    SELECT CASE ( W )
       CASE ( 1 )
          IF ( .not. State_Diag%Archive_AerosolDryAbsWL1 ) RETURN
          Slot = State_Diag%Map_AerosolDryAbsWL1%id2slot(Component)
          IF ( Slot > 0 ) State_Diag%AerosolDryAbsWL1(I,J,L,Slot) = Value
       CASE ( 2 )
          IF ( .not. State_Diag%Archive_AerosolDryAbsWL2 ) RETURN
          Slot = State_Diag%Map_AerosolDryAbsWL2%id2slot(Component)
          IF ( Slot > 0 ) State_Diag%AerosolDryAbsWL2(I,J,L,Slot) = Value
       CASE ( 3 )
          IF ( .not. State_Diag%Archive_AerosolDryAbsWL3 ) RETURN
          Slot = State_Diag%Map_AerosolDryAbsWL3%id2slot(Component)
          IF ( Slot > 0 ) State_Diag%AerosolDryAbsWL3(I,J,L,Slot) = Value
    END SELECT

  END SUBROUTINE Store_Dry_Absorption

END MODULE Dry_Absorption_Diag_Mod
