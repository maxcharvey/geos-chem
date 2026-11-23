! Mock only external interfaces; include the complete production species union.
MODULE ErrCode_Mod
 INTEGER, PARAMETER :: GC_SUCCESS=0
CONTAINS
 SUBROUTINE GC_Error(Message,RC,Location)
 CHARACTER(*),INTENT(IN) :: Message,Location
 INTEGER,INTENT(OUT) :: RC
 RC=1
 END SUBROUTINE
END MODULE
MODULE Input_Opt_Mod
 TYPE OptInput
 LOGICAL :: LBRC=.true.,ITS_A_FULLCHEM_SIM=.true.
 LOGICAL :: ITS_A_MERCURY_SIM=.false.,ITS_A_CARBON_SIM=.false.
 INTEGER :: N_Advect=3
 CHARACTER(31),ALLOCATABLE :: AdvectSpc_Name(:)
 END TYPE
END MODULE
MODULE GcKpp_Parameters
 INTEGER,PARAMETER :: NFIX=1,NSPEC=3,NVAR=2
END MODULE
MODULE GcKpp_Monitor
 CHARACTER(31) :: Spc_Names(3)=[CHARACTER(31)::'CO','OH','FIX']
END MODULE
MODULE Species_Union_Fixture
 IMPLICIT NONE
 INTEGER,PARAMETER :: MISSING_INT=-999
 CHARACTER(31),ALLOCATABLE :: Species_Names(:)
 INTEGER,ALLOCATABLE :: KppSpcId(:),KppVarId(:),KppFixId(:)
CONTAINS
 INCLUDE 'wtc_species_union_under_test.inc'
 SUBROUTINE GC_CheckVar(Name,Mode,RC)
 CHARACTER(*),INTENT(IN) :: Name
 INTEGER,INTENT(IN) :: Mode
 INTEGER,INTENT(INOUT) :: RC
 IF(RC/=0) ERROR STOP 'Fixture allocation failed'
 END SUBROUTINE
 SUBROUTINE Clear()
 IF(ALLOCATED(Species_Names)) DEALLOCATE(Species_Names)
 IF(ALLOCATED(KppSpcId)) DEALLOCATE(KppSpcId)
 IF(ALLOCATED(KppVarId)) DEALLOCATE(KppVarId)
 IF(ALLOCATED(KppFixId)) DEALLOCATE(KppFixId)
 END SUBROUTINE
END MODULE
PROGRAM WTC_SPECIES_UNION_TEST
 USE Species_Union_Fixture
 USE Input_Opt_Mod
 USE GcKpp_Monitor
 IMPLICIT NONE
 TYPE(OptInput) :: Opt
 INTEGER :: N,RC
 Opt%AdvectSpc_Name=[CHARACTER(31)::'CO','PWTC','SWTC']
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC/=0.OR.N/=6) ERROR STOP 'Fullchem diagnostic union'
 IF(ANY(Species_Names(:3)/=Opt%AdvectSpc_Name)) ERROR STOP 'Advected prefix changed'
 IF(Species_Names(N)/='WTC') ERROR STOP 'Diagnostic total missing'
 IF(KppSpcId(N)/=MISSING_INT.OR.KppVarId(N)/=MISSING_INT.OR.KppFixId(N)/=MISSING_INT) &
      ERROR STOP 'Diagnostic total registered as chemistry'
 IF(KppSpcId(1)/=1.OR.KppVarId(1)/=1.OR.KppVarId(4)/=2.OR.KppFixId(5)/=1) &
      ERROR STOP 'Existing KPP mapping changed'
 CALL Clear()
 Opt%ITS_A_FULLCHEM_SIM=.false.
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC/=0.OR.N/=4) ERROR STOP 'Non-KPP diagnostic union'
 IF(Species_Names(N)/='WTC'.OR.ANY(KppSpcId/=MISSING_INT)) ERROR STOP 'Non-KPP total indices'
 CALL Clear()
 Opt%LBRC=.false.;Opt%ITS_A_FULLCHEM_SIM=.true.;Opt%N_Advect=1
 Opt%AdvectSpc_Name=[CHARACTER(31)::'CO']
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC/=0.OR.N/=3.OR.ANY(Species_Names=='WTC')) ERROR STOP 'Ordinary model changed'
 CALL Clear()
 Opt%LBRC=.true.;Opt%N_Advect=4
 Opt%AdvectSpc_Name=[CHARACTER(31)::'CO','PWTC','SWTC','WTC']
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC==0.OR.ALLOCATED(Species_Names)) ERROR STOP 'Advected WTC accepted'
 Opt%N_Advect=2;Opt%AdvectSpc_Name=[CHARACTER(31)::'CO','SWTC']
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC==0.OR.ALLOCATED(Species_Names)) ERROR STOP 'Missing PWTC accepted'
 Opt%AdvectSpc_Name=[CHARACTER(31)::'CO','PWTC']
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC==0.OR.ALLOCATED(Species_Names)) ERROR STOP 'Missing SWTC accepted'
 Opt%N_Advect=3;Opt%AdvectSpc_Name=[CHARACTER(31)::'CO','PWTC','SWTC']
 Spc_Names(3)='WTC'
 CALL Unique_Species_Names(Opt,N,RC)
 IF(RC==0) ERROR STOP 'Chemical WTC accepted'
 CALL Clear()
 WRITE(*,'(a)') 'PASS: actual fullchem/non-KPP species union, ordinary state and invalid WTC configurations'
END PROGRAM
