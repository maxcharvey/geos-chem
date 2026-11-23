! Exercise the production routines with an isolated three-layer column.
program brc_bleaching_test
   use actual_bleaching
   implicit none
   type(OptInput) :: opt
   type(ChmState) :: chm
   type(GrdState) :: grid
   type(MetState) :: met
   type(DgnState) :: diag
   integer :: species
   real(f4) :: expected_rate

   do species=1,size(chm%Species)
      allocate(chm%Species(species)%Conc(nx,ny,nz))
   enddo

   met%T = 286.0_fp
   met%RH = 75.0_fp
   met%AD = 100.0_fp
   met%PMID(1,:,:) = 1013.25_fp
   met%PMID(2,:,:) = 506.625_fp
   met%BXHEIGHT(:,:,1) = 999.0_fp
   met%BXHEIGHT(:,:,2:3) = 1.0_fp
   chm%Species(3)%Conc = met%AD * MW_O3 / MW_AIR * 35.0e-9_fp

   call run_scheme(2)
   if (any(diag%BrCTauBleach(:,:,1:2) /= 86400.0_f4)) error stop 'Scheme 2 lifetime at or below 1 km'
   if (any(diag%BrCKBleach(:,:,1:2) /= real(1.0_fp/86400.0_fp,f4))) error stop 'Scheme 2 rate at or below 1 km'
   if (any(diag%BrCKBleach(:,:,3) /= 0.0_f4)) error stop 'Scheme 2 rate above 1 km'
   if (any(BRCSOA_CONV(:,:,1:2) <= 0.0_fp)) error stop 'Scheme 2 SWTC transfer at or below 1 km'
   if (any(NPBRC_CONV(:,:,1:2) <= 0.0_fp)) error stop 'Scheme 2 PWTC transfer at or below 1 km'
   if (any(BRCSOA_CONV(:,:,3) /= 0.0_fp)) error stop 'Scheme 2 SWTC transfer above 1 km'
   if (any(NPBRC_CONV(:,:,3) /= 0.0_fp)) error stop 'Scheme 2 PWTC transfer above 1 km'
   if (any(chm%Species(1)%Conc(:,:,3) /= 15.0_fp)) error stop 'Scheme 2 BRCSOA stock above 1 km'
   if (any(chm%Species(2)%Conc(:,:,3) /= 10.0_fp)) error stop 'Scheme 2 NPBRCPOA stock above 1 km'

   call run_scheme(0)
   if (any(BRCSOA_CONV /= 0.0_fp)) error stop 'Scheme 0 SWTC transfer'
   if (any(NPBRC_CONV /= 0.0_fp)) error stop 'Scheme 0 PWTC transfer'

   call run_scheme(3)
   expected_rate = real(1.0_fp / CALC_TAU_BRC(286.0_fp,0.75_fp,3.5e-8_fp),f4)
   if (any(diag%BrCKBleach /= expected_rate)) error stop 'Scheme 3 fixed-partial-pressure rate'
   if (any(diag%BrCO3ppbv(1,:,:) /= 35.0_f4)) error stop 'Scheme 3 O3 at 1 atm'
   if (any(diag%BrCO3ppbv(2,:,:) /= 70.0_f4)) error stop 'Scheme 3 O3 at 0.5 atm'
   if (any(NPBRC_CONV(1,:,:) /= NPBRC_CONV(2,:,:))) error stop 'Scheme 3 pressure-independent PWTC transfer'

   call run_scheme(4)
   if (any(diag%BrCO3ppbv /= 35.0_f4)) error stop 'Scheme 4 local O3 mixing ratio'
   if (any(diag%BrCKBleach(1,:,:) <= diag%BrCKBleach(2,:,:))) error stop 'Scheme 4 pressure-dependent rate'

   ozone_index = 0
   call run_scheme(4)
   if (any(diag%BrCKBleach /= expected_rate)) error stop 'Scheme 4 fallback rate'
   if (any(diag%BrCO3ppbv(2,:,:) /= 70.0_f4)) error stop 'Scheme 4 fallback O3 at 0.5 atm'

   print '(a)', 'PASS: BrC bleaching boundary, transfers, and ozone pressure diagnostics'

contains

   subroutine run_scheme(scheme)
      integer, intent(in) :: scheme
      integer :: rc

      opt%BrC_Bleach_Scheme = scheme
      chm%Species(1)%Conc = 10.0_fp
      chm%Species(2)%Conc = 10.0_fp
      FSOAS_CONV = 5.0_fp * OMOC_BBOA

      call CHEM_BRCSOA(opt,chm,diag,grid,met,1,rc)
      if (rc /= GC_SUCCESS) error stop 'BRCSOA chemistry'
      call CHEM_NPBRCPOA(opt,chm,diag,grid,met,2,rc)
      if (rc /= GC_SUCCESS) error stop 'NPBRCPOA chemistry'
   end subroutine run_scheme
end program brc_bleaching_test
