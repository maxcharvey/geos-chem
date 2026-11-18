! Test actual extracted production statements. Model interfaces are stubbed;
! this does not qualify native transport, deposition, optical tables or IO.
! Mock diagnostics are REAL8 to isolate arithmetic. Production archives are
! REAL4 and require precision-aware bounds in later full-model budgets.
PROGRAM WTC_SPLIT_TEST
 USE Actual_WTC_Code
 USE Input_Opt_Mod
 IMPLICIT NONE
 TYPE(ChmState) :: State
 TYPE(DgnState) :: Diag
 TYPE(GrdState) :: Grid
 TYPE(OptInput) :: Opt
 INTEGER :: N,RC,Mode
 REAL(fp) :: ExpectedP(2,2,2), ExpectedS(2,2,2), OldPri(2,2,2), OldSec(2,2,2)
 REAL(fp) :: OldTotal(2,2,2), Legacy(2,2,2), OA,ExpectedRadius,RadiusReference
 DO N=1,9
    ALLOCATE(State%Species(N)%Conc(2,2,2))
    State%Species(N)%Conc=REAL(N,fp)
 ENDDO
 ! Deliberately nonconsecutive IDs and unequal sources expose swapped routing.
 Legacy=State%Species(1)%Conc
 NPBRC_CONV=2.0_fp
 BRCSOA_CONV=3.0_fp
 ExpectedP=State%Species(4)%Conc+NPBRC_CONV
 ExpectedS=State%Species(2)%Conc+BRCSOA_CONV
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 IF(RC/=0) ERROR STOP 'Receiver failed'
 IF(ANY(State%Species(4)%Conc/=ExpectedP)) ERROR STOP 'Primary route'
 IF(ANY(State%Species(2)%Conc/=ExpectedS)) ERROR STOP 'Secondary route'
 IF(ANY(State%Species(1)%Conc/=Legacy)) ERROR STOP 'Legacy changed'
 IF(ANY(State%Species(9)%Conc/=9.0_fp)) ERROR STOP 'Unrelated species changed'
 IF(ANY(BRCSOA_CONV/=0.0_fp).OR.ANY(NPBRC_CONV/=0.0_fp)) ERROR STOP 'Conversion not cleared'
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 IF(ANY(State%Species(4)%Conc/=ExpectedP)) ERROR STOP 'Primary double consumption'
 IF(ANY(State%Species(2)%Conc/=ExpectedS)) ERROR STOP 'Secondary double consumption'
 ! Each formation pathway works independently, including the zero-primary case.
 BRCSOA_CONV=7.0_fp
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 ExpectedS=ExpectedS+7.0_fp
 IF(ANY(State%Species(4)%Conc/=ExpectedP)) ERROR STOP 'Secondary leaked to primary'
 IF(ANY(State%Species(2)%Conc/=ExpectedS)) ERROR STOP 'Secondary-only receiver'
 NPBRC_CONV=11.0_fp
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 ExpectedP=ExpectedP+11.0_fp
 IF(ANY(State%Species(4)%Conc/=ExpectedP)) ERROR STOP 'Primary-only receiver'
 IF(ANY(State%Species(2)%Conc/=ExpectedS)) ERROR STOP 'Primary leaked to secondary'
 ! Production cutoff is applied independently, with equality retained.
 State%Species(4)%Conc=0.0_fp
 State%Species(2)%Conc=0.0_fp
 NPBRC_CONV=SMALLNUM/2.0_fp
 BRCSOA_CONV=SMALLNUM
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 IF(ANY(State%Species(4)%Conc/=0.0_fp)) ERROR STOP 'Primary cutoff'
 IF(ANY(State%Species(2)%Conc/=SMALLNUM)) ERROR STOP 'Secondary cutoff equality'
 ! Carbon family accounting uses actual code, including OM-to-carbon conversion.
 DO N=2,8
    State%Species(N)%Conc=REAL(N,fp)
 ENDDO
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 CALL CLOSE(Diag%BrCPriMass,6.0_fp+7.0_fp+8.0_fp+4.0_fp,'Primary family')
 CALL CLOSE(Diag%BrCSecMass,3.0_fp/OMOC_BBOA+5.0_fp+2.0_fp,'Secondary family')
 CALL CLOSE(Diag%BrCWhiteMass,4.0_fp+2.0_fp,'White total')
 IF(ANY(ABS(Diag%BrCTotMass-Diag%BrCPriMass-Diag%BrCSecMass)>1e-13_fp)) &
    ERROR STOP 'Total attribution identity'
 IF(ANY(ABS(Diag%BrCBleachedFrac-Diag%BrCWhiteMass/Diag%BrCTotMass)>1e-14_fp)) &
    ERROR STOP 'Bleached fraction'
 OldPri=Diag%BrCPriMass; OldSec=Diag%BrCSecMass; OldTotal=Diag%BrCTotMass
 ! Explicitly remove carbon from both parents before receiving the same amounts.
 State%Species(6)%Conc=State%Species(6)%Conc-0.75_fp
 State%Species(5)%Conc=State%Species(5)%Conc-1.25_fp
 NPBRC_CONV=0.75_fp; BRCSOA_CONV=1.25_fp
 CALL CHEM_WTC(Opt,State,Diag,Grid,4,2,RC)
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(ABS(Diag%BrCPriMass-OldPri)>1e-13_fp)) ERROR STOP 'Primary carbon conservation'
 IF(ANY(ABS(Diag%BrCSecMass-OldSec)>1e-13_fp)) ERROR STOP 'Secondary carbon conservation'
 IF(ANY(ABS(Diag%BrCTotMass-OldTotal)>1e-13_fp)) ERROR STOP 'Total carbon conservation'
 ! A poisoned dummy WTC field must not affect any physical family or total.
 State%Species(1)%Conc=State%Species(1)%Conc+100.0_fp
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(ABS(Diag%BrCPriMass-OldPri)>1e-13_fp)) ERROR STOP 'Legacy primary ancestry'
 IF(ANY(ABS(Diag%BrCSecMass-OldSec)>1e-13_fp)) ERROR STOP 'Legacy secondary ancestry'
 IF(ANY(ABS(Diag%BrCTotMass-OldTotal)>1e-13_fp)) ERROR STOP 'Diagnostic WTC double counting'
 ! Optional parent IDs must never be indexed when absent.
 CALL FAMILY_MASS(State,Diag,Grid,.false.)
 IF(ANY(Diag%BrCPriMass/=State%Species(4)%Conc)) ERROR STOP 'Optional primary parents'
 ! Disabled archives leave data unchanged; each new archive activates accounting alone.
 Diag%Archive_BrCAbsMass=.false.; Diag%Archive_BrCTotMass=.false.
 Diag%Archive_BrCBleachedFrac=.false.; Diag%Archive_BrCWhiteMass=.false.
 Diag%Archive_BrCPriMass=.false.; Diag%Archive_BrCSecMass=.false.
 Diag%BrCPriMass=-1.0_fp; Diag%BrCSecMass=-2.0_fp; Diag%BrCWhiteMass=-3.0_fp
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(Diag%BrCPriMass/=-1.0_fp).OR.ANY(Diag%BrCSecMass/=-2.0_fp) &
    .OR.ANY(Diag%BrCWhiteMass/=-3.0_fp)) ERROR STOP 'Disabled archive modified'
 Diag%Archive_BrCPriMass=.true.
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(ABS(Diag%BrCPriMass-OldPri)>1e-13_fp)) ERROR STOP 'Primary-only archive'
 Diag%Archive_BrCPriMass=.false.; Diag%Archive_BrCSecMass=.true.
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(ABS(Diag%BrCSecMass-OldSec)>1e-13_fp)) ERROR STOP 'Secondary-only archive'
 Diag%Archive_BrCSecMass=.false.; Diag%Archive_BrCWhiteMass=.true.
 CALL FAMILY_MASS(State,Diag,Grid,.true.)
 IF(ANY(ABS(Diag%BrCWhiteMass-8.0_fp)>1e-13_fp)) ERROR STOP 'White-only archive'
 ! Exercise actual aggregate and PDER expressions with physical components only.
 State%AerMass%OCFOPOA=2.1_fp
 State%AerMass%OCFPOA=1.4_fp
 State%AerMass%SO4_NH4_NIT=10.0_fp
 State%AerMass%OCPO=2.0_fp
 State%AerMass%OCPISOA=3.0_fp
 State%AerMass%FSOAS=3.0_fp
 State%AerMass%PBRC=14.7_fp
 State%Species(3)%Conc=3.0_fp
 State%Species(7)%Conc=7.0_fp
 State%Species(5)%Conc=5.0_fp
 State%Species(6)%Conc=6.0_fp
 OA=2.0_fp+3.0_fp+5.0_fp*2.1_fp+6.0_fp*2.1_fp+2.0_fp*7.0_fp*2.1_fp+3.0_fp
 ExpectedRadius=EXP(4.36_fp+0.20_fp*LOG((10.0_fp+OA)*1e9_fp) &
                   +0.065_fp*LOG(OA/10.0_fp))*0.001_fp/0.9_fp
 DO Mode=0,2
    State%Species(1)%Conc=0.0_fp
    State%Species(4)%Conc=0.0_fp
    State%Species(2)%Conc=0.0_fp
    SELECT CASE(Mode)
    CASE(0)
       State%Species(1)%Conc=1.0_fp
       State%Species(4)%Conc=4.0_fp
       State%Species(2)%Conc=3.0_fp
    CASE(1)
       State%Species(4)%Conc=7.0_fp
    CASE(2)
       State%Species(2)%Conc=7.0_fp
    END SELECT
    CALL MASS_AND_PDER(State,Opt)
    CALL CLOSE(State%AerMass%WTCPI,14.7_fp,'Aggregate white OM')
    CALL CLOSE(State%AerMass%R_OMSNA,OA/10.0_fp,'PDER OM to SNA')
    IF(ANY(ABS(State%AerMass%SNAOM-(10.0_fp+OA)*1e9_fp)>1e-12_fp*(10.0_fp+OA)*1e9_fp)) &
       ERROR STOP 'PDER total SNA plus OM'
    CALL CLOSE(State%AerMass%PDER,ExpectedRadius,'PDER radius with all white origins')
    IF(Mode==0) RadiusReference=State%AerMass%PDER(1,1,1)
    CALL CLOSE(State%AerMass%PDER,RadiusReference,'PDER origin repartition neutrality')
 ENDDO
 State%Species(1)%Conc=1.0e9_fp
 CALL MASS_AND_PDER(State,Opt)
 CALL CLOSE(State%AerMass%PDER,RadiusReference,'Dummy WTC excluded from PDER')
 State%Species(4)%Conc=0.0_fp
 State%Species(2)%Conc=0.0_fp
 CALL MASS_AND_PDER(State,Opt)
 IF(ANY(State%AerMass%PDER>=RadiusReference)) ERROR STOP 'White OM omitted from PDER sensitivity'
 State%AerMass%SO4_NH4_NIT=0.0_fp
 CALL MASS_AND_PDER(State,Opt)
 CALL CLOSE(State%AerMass%PDER,0.005_fp,'Zero SNA PDER fallback')
 DO N=1,9
    DEALLOCATE(State%Species(N)%Conc)
 ENDDO
 WRITE(*,'(a)') 'PASS: actual WTC routing, carbon families, aggregate OA and PDER'
CONTAINS
 SUBROUTINE CLOSE(Value,Expected,Label)
 REAL(fp),INTENT(IN) :: Value(:,:,:),Expected
 CHARACTER(*),INTENT(IN) :: Label
 IF(ANY(ABS(Value-Expected)>1e-13_fp)) THEN
    WRITE(*,*) Label,Value,Expected
    ERROR STOP 'Mass expectation failed'
 ENDIF
 END SUBROUTINE
END PROGRAM
