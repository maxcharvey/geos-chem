! Experimental FINNv25 origin overlays, outside physical aerosol/optical bins.
! Requires explicit complete 7 x 4 registration and a qualified origin restart.
MODULE BRC_ORIGIN_MOD
  USE Precision_Mod, ONLY: fp, f8
  USE State_Chm_Mod, ONLY: ChmState, Ind_
  USE State_Grid_Mod, ONLY: GrdState
  USE State_Met_Mod, ONLY: MetState
  USE Input_Opt_Mod, ONLY: OptInput
  USE ErrCode_Mod
  USE UnitConv_Mod, ONLY: KG_SPECIES, KG_SPECIES_PER_KG_DRY_AIR, &
       KG_SPECIES_PER_M2, MOLES_SPECIES_PER_MOLES_DRY_AIR
  USE PhysConstants, ONLY: AIRMW
  USE BRC_ORIGIN_KERNEL_MOD, ONLY: BRC_ORIGIN_SPLIT_LOSS
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: BRC_ORIGIN_BIND, BRC_ORIGIN_TRANSFER, BRC_ORIGIN_AUDIT
  CHARACTER(LEN=8), PARAMETER :: Names(7) = [ CHARACTER(LEN=8) :: &
       'FSOAP','FSOAS','BRCSOA','NPBRCPOA','WTC','PBRCPOA','DBRCPOA' ]
  CHARACTER(LEN=3), PARAMETER :: Origins(4) = ['USA','CAN','ROW','UNT']
  INTEGER, SAVE :: Parents(7)=0, Tags(7,4)=0
  LOGICAL, SAVE :: Bound=.FALSE., Enabled=.FALSE., BadTransition=.FALSE.
CONTAINS
  SUBROUTINE BRC_ORIGIN_BIND(RC)
    INTEGER, INTENT(OUT) :: RC
    INTEGER :: S, O, Count
    RC=GC_SUCCESS
    IF (Bound) RETURN
    Count=0
    DO S=1,7
      Parents(S)=Ind_(TRIM(Names(S)))
      DO O=1,4
        Tags(S,O)=Ind_(TRIM(Names(S))//'_'//Origins(O))
        IF (Tags(S,O)>0) Count=Count+1
      ENDDO
    ENDDO
    IF (Count /= 0 .AND. (Count /= 28 .OR. ANY(Parents<=0))) THEN
      WRITE(6,*) 'BRC_ORIGIN_ERROR: require complete seven-parent/four-origin set'
      RC=GC_FAILURE
      RETURN
    ENDIF
    Enabled=Count==28
    Bound=.TRUE.
  END SUBROUTINE BRC_ORIGIN_BIND

  ! Called within each parent's existing cell loop after the actual CNEW is
  ! known. Condensed/darkened origin material reaches the next source before
  ! that source's own loop, just as the parent conversion arrays do.
  SUBROUTINE BRC_ORIGIN_TRANSFER(State_Chm, SourceId, DestId, I,J,L,Before,After,Factor)
    TYPE(ChmState), INTENT(INOUT) :: State_Chm
    INTEGER, INTENT(IN) :: SourceId, DestId, I,J,L
    REAL(fp), INTENT(IN) :: Before, After, Factor
    REAL(fp) :: P(4), Lost(4)
    INTEGER :: S,D,O,Status
    IF (.NOT. Enabled) RETURN
    IF (.NOT. IEEE_IS_FINITE(Factor) .OR. Factor<0.0_fp) THEN
      !$OMP ATOMIC WRITE
      BadTransition=.TRUE.
      RETURN
    ENDIF
    S=0; D=0
    DO O=1,7
      IF (Parents(O)==SourceId) S=O
      IF (Parents(O)==DestId) D=O
    ENDDO
    IF (S==0 .OR. D==0) THEN
      !$OMP ATOMIC WRITE
      BadTransition=.TRUE.
      RETURN
    ENDIF
    DO O=1,4
      P(O)=State_Chm%Species(Tags(S,O))%Conc(I,J,L)
    ENDDO
    CALL BRC_ORIGIN_SPLIT_LOSS(Before,After,P,Lost,Status)
    IF (Status/=0) THEN
      !$OMP ATOMIC WRITE
      BadTransition=.TRUE.
      RETURN
    ENDIF
    DO O=1,4
      State_Chm%Species(Tags(S,O))%Conc(I,J,L)=P(O)
      State_Chm%Species(Tags(D,O))%Conc(I,J,L)= &
        State_Chm%Species(Tags(D,O))%Conc(I,J,L)+Lost(O)*Factor
    ENDDO
  END SUBROUTINE BRC_ORIGIN_TRANSFER

  ! Logs mass-weighted aggregate closure without modifying any state. Native
  ! concentration units are verified and converted to kg only in scalar copies.
  SUBROUTINE BRC_ORIGIN_AUDIT(Label,Input_Opt,State_Chm,State_Grid,State_Met,RC)
    CHARACTER(LEN=*), INTENT(IN) :: Label
    TYPE(OptInput), INTENT(IN) :: Input_Opt
    TYPE(ChmState), INTENT(IN) :: State_Chm
    TYPE(GrdState), INTENT(IN) :: State_Grid
    TYPE(MetState), INTENT(IN) :: State_Met
    INTEGER, INTENT(OUT) :: RC
    INTEGER :: S,O,I,J,L,U,N
    REAL(f8) :: Weight, ParentMass, PartMass(4), Residual
    REAL(f8) :: Tot, L1, Signed, MaxAbs, ArcticL1, OriginMass(4)
    CALL BRC_ORIGIN_BIND(RC)
    IF (RC/=GC_SUCCESS .OR. .NOT. Enabled) RETURN
    IF (BadTransition) THEN
      WRITE(6,*) 'BRC_ORIGIN_ERROR: invalid chemistry transition'
      RC=GC_FAILURE
      RETURN
    ENDIF
    DO S=1,7
      N=Parents(S); U=State_Chm%Species(N)%Units
      DO O=1,4
        IF (State_Chm%Species(Tags(S,O))%Units/=U) THEN
          WRITE(6,*) 'BRC_ORIGIN_ERROR: origin and parent unit mismatch'
          RC=GC_FAILURE
          RETURN
        ENDIF
      ENDDO
      Tot=0; L1=0; Signed=0; MaxAbs=0; ArcticL1=0; OriginMass=0
      DO L=1,State_Grid%NZ
      DO J=1,State_Grid%NY
      DO I=1,State_Grid%NX
        SELECT CASE(U)
        CASE(KG_SPECIES)
          Weight=1.0_f8
        CASE(KG_SPECIES_PER_KG_DRY_AIR)
          Weight=State_Met%AD(I,J,L)
        CASE(KG_SPECIES_PER_M2)
          Weight=State_Met%AREA_M2(I,J)
        CASE(MOLES_SPECIES_PER_MOLES_DRY_AIR)
          Weight=State_Met%AD(I,J,L)*State_Chm%SpcData(N)%Info%MW_g/AIRMW
        CASE DEFAULT
          WRITE(6,*) 'BRC_ORIGIN_ERROR: unsupported audit units',U
          RC=GC_FAILURE
          RETURN
        END SELECT
        ParentMass=REAL(State_Chm%Species(N)%Conc(I,J,L),f8)*Weight
        DO O=1,4
          PartMass(O)=REAL(State_Chm%Species(Tags(S,O))%Conc(I,J,L),f8)*Weight
        ENDDO
        IF (.NOT. IEEE_IS_FINITE(ParentMass) .OR. ParentMass<0.0_f8 .OR. &
            ANY(.NOT. IEEE_IS_FINITE(PartMass)) .OR. ANY(PartMass<0.0_f8)) THEN
          WRITE(6,*) 'BRC_ORIGIN_ERROR: nonfinite or negative origin state',TRIM(Names(S)),I,J,L
          RC=GC_FAILURE
          RETURN
        ENDIF
        Residual=SUM(PartMass)-ParentMass
        Tot=Tot+ABS(ParentMass); Signed=Signed+Residual; L1=L1+ABS(Residual)
        MaxAbs=MAX(MaxAbs,ABS(Residual)); OriginMass=OriginMass+PartMass
        IF (State_Grid%Lat(J)>=66.5_f8) ArcticL1=ArcticL1+ABS(Residual)
      ENDDO
      ENDDO
      ENDDO
      IF (Input_Opt%amIRoot) WRITE(6,'(a,1x,a,1x,a,1x,i1,9(1x,es22.14))') &
        'BRC_ORIGIN_AUDIT',TRIM(Label),TRIM(Names(S)),U,Tot,Signed,L1,MaxAbs,ArcticL1,OriginMass
    ENDDO
  END SUBROUTINE BRC_ORIGIN_AUDIT
END MODULE BRC_ORIGIN_MOD
