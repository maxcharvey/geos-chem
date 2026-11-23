#!/usr/bin/env python3
"""Extract production fire OC aging, aerosol aggregation and GFED routing."""
from pathlib import Path
import re
import sys
ROOT = Path(__file__).resolve().parents[2]
carbon = (ROOT / 'GeosCore/carbon_mod.F90').read_text()
aging = re.search(r'^ SUBROUTINE CHEM_FOC\(.*?^ END SUBROUTINE CHEM_FOC', carbon, re.M | re.S).group()
aerosol = (ROOT / 'GeosCore/aerosol_mod.F90').read_text()
aggregate = aerosol.split('          ! Source-tagged residual fire OC shares the existing OC carriers.\n')[1].split('       ENDIF ! LCARB')[0]
gfed = (ROOT.parent / 'HEMCO/src/Extensions/hcox_gfed_mod.F90').read_text()
gates = gfed.split('       ! Skip branch-specific BrC emissions when disabled.')[1].split('       ! SpcArr are the total biomass')[0]
gates = gates[gates.index('       IF ( .NOT. Inst%UseBrC )'):]
speciation = gfed.split('       SELECT CASE ( Inst%SpcNames(N) )')[2].split('       END SELECT')[0]
speciation = '       SELECT CASE ( Inst%SpcNames(N) )' + speciation + '       END SELECT\n'
inherit = gfed.split('    ! Residual fire tags replace OC, so inherit its inventory scaling and mask.')[1].split('    ! SOAP and FSOAP are specified')[0]
stubs = '''
module Precision_Mod
 integer, parameter :: fp=kind(1d0), sp=kind(1d0), hp=kind(1d0)
end module
module ErrCode_Mod
 integer, parameter :: GC_SUCCESS=0
end module
module State_Grid_Mod
 type GrdState
 integer :: NX=2,NY=1,NZ=2
 end type
end module
module State_Chm_Mod
 use Precision_Mod
 type species
 real(fp), allocatable :: Conc(:,:,:)
 end type
 type aerosol
 real(fp) :: OCPI(2,1,2),OCPO(2,1,2),OCFOPOA(2,1),OCFPOA(2,1)
 end type
 type ChmState
 type(species) :: Species(4)
 type(aerosol) :: AerMass
 end type
end module
module Time_Mod
 use Precision_Mod
 real(fp) :: dt=3600d0
 contains
 function GET_TS_CHEM() result(value)
 real(fp) :: value
 value=dt
 end function
end module
module actual
 use Precision_Mod
 real(fp), parameter :: SMALLNUM=1d-30
 contains
'''
wrappers = '''
 subroutine aggregate(State_Chm, enabled)
 use State_Chm_Mod
 type(ChmState), intent(inout) :: State_Chm
 logical, intent(in) :: enabled
 type options
 logical :: LBRC
 end type
 type(options) :: Input_Opt
 type(species) :: Spc(4)
 integer :: I,J,L, id_FOCPI=1,id_FOCPO=2
 real(fp) :: AIRVOL(2,1,2)
 Input_Opt%LBRC=enabled
 Spc=State_Chm%Species
 AIRVOL=2d0
 DO L=1,2
 DO J=1,1
 DO I=1,2
''' + aggregate + '''
 ENDDO
 ENDDO
 ENDDO
 end subroutine
 subroutine routing(enabled, replacement, values)
 logical,intent(in) :: enabled
 real(fp),intent(in) :: replacement
 real(fp),intent(out) :: values(8)
 type config
 logical :: UseBrC
 character(31) :: SpcNames(8)
 real(fp) :: OCPIfrac=.5d0,NPBRCPOAfrac=.5d0,PBRCPOAfrac=.25d0
 real(fp) :: BCPIfrac=.2d0,POG1frac=.5d0,SOAPfrac=.013d0,FSOAPfrac,DBRCPOAfrac=4d0
 end type
 type(config) :: Inst
 integer :: N
 real(fp) :: SpcArr
 Inst%UseBrC=enabled
 Inst%FSOAPfrac=replacement
 Inst%SpcNames=[character(31)::'OCPI','OCPO','FOCPI','FOCPO','NPBRCPOA','PBRCPOA','SOAP','FSOAP']
 values=0d0
 DO N=1,8
''' + gates + '''
 SpcArr=1d0
''' + speciation + '''
 values(N)=SpcArr
 ENDDO
 end subroutine
 subroutine HCO_ERROR(message,rc)
 character(*),intent(in) :: message
 integer,intent(out) :: rc
 rc=1
 end subroutine
 subroutine check_inheritance(conflict,rc)
 logical,intent(in) :: conflict
 integer,intent(out) :: rc
 character(31) :: tmpName
 character(61),parameter :: HCOX_NOSCALE='none'
 logical :: Matched
 integer :: N,M
 type cfg
 logical :: UseBrC=.true.
 integer :: nSpc=4
 character(31) :: SpcNames(4)
 real(fp) :: SpcScal(4)
 character(61) :: SpcScalFldNme(4)
 end type
 type(cfg) :: Inst
 Inst%SpcNames=[character(31)::'OCPI','OCPO','FOCPI','FOCPO']
 Inst%SpcScal=[2d0,3d0,1d0,1d0]
 Inst%SpcScalFldNme=[character(61)::'PI-mask','PO-mask','none','none']
 if(conflict) Inst%SpcScal(3)=4d0
 rc=0
''' + inherit + '''
 if(any(Inst%SpcScal(3:4)/=Inst%SpcScal(1:2))) stop 14
 if(any(Inst%SpcScalFldNme(3:4)/=Inst%SpcScalFldNme(1:2))) stop 15
 end subroutine
end module
'''
driver = '''
program check
 use actual
 use State_Chm_Mod
 use State_Grid_Mod
 use Time_Mod
 type(ChmState) :: chm
 type(GrdState) :: grid
 real(fp) :: vals(8), before(2,1,2), donor(2,1,2), expected(2,1,2)
 real(fp) :: rate, initial(2,1,2), tagPI(2,1,2), tagPO(2,1,2)
 integer :: n,rc,k
 do n=1,4
 allocate(chm%Species(n)%Conc(2,1,2))
 chm%Species(n)%Conc=real(n,fp)
 enddo
 do k=0,3
 dt=real(k,fp)*1800d0
 chm%Species(1)%Conc=2d0
 chm%Species(2)%Conc=3d0
 before=chm%Species(1)%Conc+chm%Species(2)%Conc
 donor=chm%Species(2)%Conc
 expected=donor*exp(-dt/(86400d0*1.15d0))
 call CHEM_FOC(chm,grid,2,1,rc)
 if(rc/=0) stop 1
 if(maxval(abs(chm%Species(2)%Conc-expected))>2d-15) stop 2
 if(maxval(abs(chm%Species(1)%Conc+chm%Species(2)%Conc-before))>2d-15) stop 3
 enddo
 chm%Species(1)%Conc=0d0
 chm%Species(2)%Conc=SMALLNUM/2d0
 call CHEM_FOC(chm,grid,2,1,rc)
 if(any(chm%Species(1)%Conc/=0d0).or.any(chm%Species(2)%Conc/=0d0)) stop 4
 chm%Species(1)%Conc=4d0
 chm%Species(2)%Conc=6d0
 tagPI=chm%Species(1)%Conc
 tagPO=chm%Species(2)%Conc
 chm%AerMass%OCPI=10d0
 chm%AerMass%OCPO=20d0
 chm%AerMass%OCFOPOA=2.1d0
 chm%AerMass%OCFPOA=1.4d0
 call aggregate(chm,.false.)
 if(any(chm%AerMass%OCPI/=10d0).or.any(chm%AerMass%OCPO/=20d0)) stop 5
 call aggregate(chm,.true.)
 if(maxval(abs(chm%AerMass%OCPI-(10d0+tagPI*2.1d0/2d0)))>2d-15) stop 6
 if(maxval(abs(chm%AerMass%OCPO-(20d0+tagPO*1.4d0/2d0)))>2d-15) stop 7
 if(any(chm%Species(1)%Conc/=tagPI).or.any(chm%Species(2)%Conc/=tagPO)) stop 8
 call routing(.false.,.013d0,vals)
 if(any(vals(3:6)/=0d0).or.vals(8)/=0d0) stop 9
 if(abs(sum(vals(1:2))-1d0)>1d-15.or.vals(7)/=.013d0) stop 10
 call routing(.true.,.013d0,vals)
 if(any(vals(1:2)/=0d0).or.vals(7)/=0d0) stop 11
 if(abs(sum(vals(3:6))-1d0)>1d-15.or.vals(8)/=.013d0) stop 12
 call routing(.true.,0d0,vals)
 if(vals(7)/=.013d0.or.vals(8)/=0d0) stop 13
 call check_inheritance(.false.,rc)
 if(rc/=0) stop 16
 call check_inheritance(.true.,rc)
 if(rc/=1) stop 17
 print *, 'PASS: production fire OC aging, conservation, OC aggregation and GFED exclusive routing'
end program
'''
smallnum = re.search(r'PARAMETER\s*:: SMALLNUM\s*=\s*(\S+)', carbon).group(1)
Path(sys.argv[1]).write_text(stubs.replace('SMALLNUM=1d-30', 'SMALLNUM='+smallnum)+aging+'\n'+wrappers+driver)
