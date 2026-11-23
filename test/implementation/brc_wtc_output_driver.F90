PROGRAM BRC_WTC_OUTPUT_TEST
  USE Precision_Mod, ONLY: fp, f4
  USE State_Chm_Mod, ONLY: ChmState
  USE BRC_WTC_Output_Mod
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_VALUE, IEEE_QUIET_NAN
  IMPLICIT NONE
  TYPE(ChmState) :: State
  TYPE(BRC_WTC_Map) :: Map, Ordinary
  REAL(fp) :: Cube(2,2,2), Slab(2,2), P(2,2,2), S(2,2,2), Total(2,2,2), Before(2,2,2)
  REAL(fp) :: Average(2,2,2)
  INTEGER :: N
  DO N=1,9
     ALLOCATE(State%Species(N)%Conc(2,2,2))
     State%Species(N)%Conc=REAL(N,fp)
  ENDDO
  Map=BRC_WTC_Output_Ids()
  IF(Map%Total/=1.OR.Map%Primary/=4.OR.Map%Secondary/=2) ERROR STOP 'Output map'
  ! The dummy is deliberately unrelated to both physical stocks.
  State%Species(1)%Conc=1e9_fp
  Before=State%Species(1)%Conc
  CALL BRC_Copy_Species_Output(State,1,Cube,Map)
  IF(ANY(Cube/=6.0_fp)) ERROR STOP '3D derived output'
  CALL BRC_Copy_Species_Output(State,3,Cube,Map)
  IF(ANY(Cube/=3.0_fp)) ERROR STOP 'Ordinary species copy'
  IF(BRC_Species_Output_Cell(State,1,2,1,2,Map)/=6.0_fp) ERROR STOP 'PF/surface derived cell'
  CALL BRC_Copy_Species_Slab(State,1,1,Slab,Map)
  IF(ANY(Slab/=6.0_fp)) ERROR STOP 'Satellite derived slab'
  IF(ANY(State%Species(1)%Conc/=Before)) ERROR STOP 'Output modified dummy state'
  IF(ANY(State%Species(4)%Conc/=4.0_fp).OR.ANY(State%Species(2)%Conc/=2.0_fp)) &
       ERROR STOP 'Output modified physical stock'
  ! Mimic a transport/deposition change between sampling instants.
  State%Species(4)%Conc(1,1,1)=0.5_fp
  State%Species(2)%Conc(1,1,1)=0.25_fp
  CALL BRC_Copy_Species_Output(State,1,Cube,Map)
  IF(Cube(1,1,1)/=0.75_fp.OR.Cube(2,2,2)/=6.0_fp) ERROR STOP 'Stale process output'
  IF(BRC_Species_Output_Cell(State,1,1,1,1,Map)/=0.75_fp) ERROR STOP 'Stale PF output'
  ! Each consumer sees its components' actual current units, not dummy units.
  State%Species(4)%Conc=State%Species(4)%Conc*1e-9_fp
  State%Species(2)%Conc=State%Species(2)%Conc*1e-9_fp
  CALL BRC_Copy_Species_Output(State,1,Cube,Map)
  IF(ANY(ABS(Cube-(State%Species(4)%Conc+State%Species(2)%Conc))>1e-25_fp)) &
       ERROR STOP 'Unit lifecycle output'
  CALL BRC_Copy_Species_Slab(State,1,1,Slab,Map)
  IF(ANY(ABS(Slab-Cube(1,:,:))>1e-25_fp)) ERROR STOP 'Slab units'
  ! Time averaging commutes with the sum when the samples have equal cadence.
  Average=Cube*0.25_fp
  State%Species(4)%Conc=State%Species(4)%Conc*2.0_fp
  State%Species(2)%Conc=State%Species(2)%Conc*3.0_fp
  CALL BRC_Copy_Species_Output(State,1,Cube,Map)
  Average=Average+Cube*0.75_fp
  P=State%Species(4)%Conc/2.0_fp*0.25_fp+State%Species(4)%Conc*0.75_fp
  S=State%Species(2)%Conc/3.0_fp*0.25_fp+State%Species(2)%Conc*0.75_fp
  IF(ANY(ABS(Average-(P+S))>1e-22_fp)) ERROR STOP 'Sampled mean total'
  CALL BRC_Copy_Species_Output(State,1,Cube,Ordinary)
  IF(ANY(Cube/=Before)) ERROR STOP 'No-WTC map changed ordinary output'

  IF(LEN_TRIM(BRC_WTC_Missing_Restart_Error('PWTC',.false.))==0) ERROR STOP 'Missing primary accepted'
  IF(LEN_TRIM(BRC_WTC_Missing_Restart_Error('SWTC',.false.))==0) ERROR STOP 'Missing secondary accepted'
  IF(LEN_TRIM(BRC_WTC_Missing_Restart_Error('SWTC',.true.))/=0) ERROR STOP 'Present secondary refused'
  IF(LEN_TRIM(BRC_WTC_Missing_Restart_Error('NO2',.false.))/=0) ERROR STOP 'Ordinary fallback refused'
  P=4.0_fp; S=2.0_fp; Total=6.0_fp
  IF(.NOT.BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Matching restart refused'
  Total=7.0_fp
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Mismatched total accepted'
  P=REAL(REAL(1.1_fp,f4),fp); S=REAL(REAL(2.2_fp,f4),fp)
  Total=REAL(REAL(3.3_fp,f4),fp)
  IF(.NOT.BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'REAL4 serialization refused'
  P=0.0_fp; S=0.0_fp; Total=0.0_fp
  IF(.NOT.BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Explicit fresh zero refused'
  Total=1e-18_fp
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Unpartitioned legacy mass accepted'
  Total=0.0_fp; P=-1e-30_fp
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Negative stock accepted'
  P=IEEE_VALUE(0.0_fp,IEEE_QUIET_NAN)
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'NaN accepted'
  P=HUGE(1.0_fp); S=HUGE(1.0_fp); Total=HUGE(1.0_fp)
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'Overflowing component sum accepted'
  IF(BRC_WTC_Components_Valid(P,S)) ERROR STOP 'Missing-total overflow accepted'
  P=0.0_fp; S=0.0_fp; Total=IEEE_VALUE(0.0_fp,IEEE_QUIET_NAN)
  IF(BRC_WTC_Total_Matches(Total,P,S)) ERROR STOP 'NaN total accepted'
  Total=0.0_fp
  IF(BRC_WTC_Total_Matches(Total(:1,:1,:1),P,S)) ERROR STOP 'Shape mismatch accepted'
  DO N=1,9
     DEALLOCATE(State%Species(N)%Conc)
  ENDDO
  WRITE(*,'(a)') 'PASS: derived WTC sampling, live units, no double stock, restart provenance/serialization'
END PROGRAM
