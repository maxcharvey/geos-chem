#!/usr/bin/env python3
"""Static lifecycle checks and extraction of actual WTC code for native tests."""
from pathlib import Path
import re
import sys
import unittest

ROOT = Path(__file__).resolve().parents[2]


def blocks():
    text = (ROOT / 'run/shared/species_database.yml').read_text()
    result = {}
    for name, body in re.findall(r'^([A-Za-z0-9_]+):\n(.*?)(?=^[A-Za-z0-9_]+:|\Z)', text, re.M | re.S):
        result[name] = dict(re.findall(r'^  (\w+):\s*(.*?)\s*$', body, re.M))
    return result


class SplitContractTest(unittest.TestCase):
    def test_species_physics_and_zero_provenance_backgrounds(self):
        db = blocks()
        # Is_HygroGrowth is deliberately a single aggregate-bin registration.
        metadata = {'FullName', 'Is_HygroGrowth', 'Background_VV'}
        reference = {k: v for k, v in db['PWTC'].items() if k not in metadata}
        for name in ('PWTC', 'SWTC'):
            self.assertEqual({k: v for k, v in db[name].items() if k not in metadata}, reference)
            self.assertEqual(db[name]['Is_HygroGrowth'], 'false')
            self.assertEqual(float(db[name]['Background_VV']), 0.0)
        self.assertEqual(db['WTC']['Is_HygroGrowth'], 'true')
        for key in ('Density','MW_g'):
            self.assertEqual(db['WTC'][key], db['PWTC'][key])
            self.assertEqual(db['WTC'][key], db['SWTC'][key])
        for key in ('Is_Gas', 'Is_Aerosol', 'Is_Advected', 'Is_DryDep', 'Is_WetDep'):
            self.assertEqual(db['WTC'][key], 'false')
        self.assertEqual(float(db['WTC']['Background_VV']), 0.0)

    def test_templates_register_both_tracers(self):
        for simulation in ('fullchem', 'aerosol'):
            text = (ROOT / f'run/GCClassic/geoschem_config.yml.templates/geoschem_config.yml.{simulation}').read_text()
            for name in ('PWTC', 'SWTC'):
                self.assertEqual(len(re.findall(rf'^      - {name}(?:\s|$)', text, re.M)), 1)
        for simulation in ('fullchem', 'carbon'):
            text = (ROOT / f'run/GCClassic/HISTORY.rc.templates/HISTORY.rc.{simulation}').read_text()
            for name in ('BrCWhiteMass', 'BrCPriMass', 'BrCSecMass'):
                self.assertIn("'" + name, text)

    def test_total_is_registered_without_transport(self):
        for simulation in ('fullchem','aerosol'):
            text = (ROOT / f'run/GCClassic/geoschem_config.yml.templates/geoschem_config.yml.{simulation}').read_text()
            self.assertNotRegex(text, re.compile(r'^      - WTC(?:\s|$)', re.M))
        reader = (ROOT / 'Headers/species_database_mod.F90').read_text()
        self.assertIn("Species_Names = [CHARACTER(LEN=31) :: Species_Names, 'WTC']", reader)
        for key in ('KppSpcId','KppVarId','KppFixId'):
            self.assertIn(f'{key} = [{key}, MISSING_INT]', reader)
        self.assertIn("ANY(Input_Opt%AdvectSpc_Name == 'WTC')", reader)

    def test_diagnostic_lifecycle(self):
        text = (ROOT / 'Headers/state_diag_mod.F90').read_text()
        for name in ('BrCWhiteMass', 'BrCPriMass', 'BrCSecMass'):
            self.assertRegex(text, rf'POINTER :: {name}\(:,:,:\)')
            self.assertRegex(text, rf'State_Diag%{name}\s*=> NULL\(\)')
            self.assertRegex(text, rf'State_Diag%Archive_{name}\s*= .FALSE.')
            self.assertIn(f"diagID = '{name}'", text)
            self.assertIn(f"CALL Finalize( diagId   = '{name}'", text)
            self.assertIn(f"Name_AllCaps ) == '{name.upper()}'", text)

    def test_optics_aggregate_and_pm_origin(self):
        text = (ROOT / 'GeosCore/aerosol_mod.F90').read_text()
        self.assertIn('State_Chm%AerMass%PWTCPI(I,J,L) + State_Chm%AerMass%SWTCPI(I,J,L)', text)
        self.assertIn('State_Chm%AerMass%WAERSL(I,J,L,8) = State_Chm%AerMass%WTCPI(I,J,L)', text)
        primary = text.split('       ! PM2.5 OC\n')[1].split('       ! PM2.5 dust')[0]
        secondary = text.split('       ! PM2.5 SOA\n')[1].split('       ! PM2.5 nitrate')[0]
        self.assertIn('%PWTCPI(', primary)
        self.assertNotIn('%SWTCPI(', primary)
        self.assertNotIn('%WTCPI(', primary)
        self.assertIn('%SWTCPI(', secondary)
        self.assertNotIn('%PWTCPI(', secondary)
        self.assertNotIn('%WTCPI(', secondary)
        container = (ROOT / 'Headers/aermass_container_mod.F90').read_text()
        for name in ('PWTCPI', 'SWTCPI'):
            self.assertIn(f'ALLOCATE( Aer%{name}(', container)
            self.assertIn(f'DEALLOCATE( Aer%{name},', container)
            self.assertIn(f'Aer%{name} => NULL()', container)

    def test_pder_and_radiation_consumers(self):
        aerosol = (ROOT / 'GeosCore/aerosol_mod.F90').read_text()
        pder = aerosol.split('      ! PDER [um]')[1].split('#ifdef MODEL_GEOS')[0]
        self.assertEqual(pder.count('%WTCPI(I,J,L)'), 2)
        self.assertNotIn('%PWTCPI(', pder)  # Already included in WTCPI; never add twice.
        self.assertNotIn('%SWTCPI(', pder)
        cloudj = (ROOT / 'GeosCore/cldj_interface_mod.F90').read_text()
        self.assertIn('AERSP(L,S_rhx) = ( State_Chm%AerMass%WTCPI(I,J,L)', cloudj)
        self.assertIn('RTODAER(I,J,L,IWV,NRT)     = ODAER(I,J,L,IWV,N)', aerosol)
        self.assertIn('NRT=N+2', aerosol)
        rrtmg = (ROOT / 'GeosCore/rrtmg_rad_transfer_mod.F90').read_text()
        white = rrtmg.split('       CASE( 20 )')[1].split('       CASE( 21 )')[0]
        self.assertIn('SPECMASK(10) = 20', white)
        total_brc = rrtmg.split('       CASE( 24 )')[1].split('       END SELECT')[0]
        self.assertIn('DO II = 8, 13', total_brc)
        total_pm = rrtmg.split('       CASE( 15 )')[1].split('       CASE( 16 )')[0]
        self.assertIn('DO II = 1, State_Chm%Phot%NASPECRAD', total_pm)

    def test_sampling_paths_and_nonphysical_total(self):
        for relative, patterns in {
            'GeosCore/diagnostics_mod.F90': ('State_Diag%SpeciesConcVV', 'State_Diag%SpeciesBC', 'State_Diag%SpeciesRst', 'State_Diag%SpeciesConcMND'),
            'History/history_mod.F90': ('State_Diag%SpeciesBC',),
            'GeosCore/fullchem_mod.F90': ('State_Diag%ConcBeforeChem', 'State_Diag%ConcAfterChem'),
        }.items():
            text = (ROOT / relative).read_text()
            for pattern in patterns:
                self.assertRegex(text, 'CALL BRC_Copy_Species_Output\\(State_Chm, N, ' + re.escape(pattern))
        text = (ROOT / 'GeosCore/planeflight_mod.F90').read_text()
        self.assertEqual(text.count('VARI(V) = BRC_Species_Output_Cell(State_Chm, N, I, J, L, WtcOutput)'), 2)
        aerosol = (ROOT / 'GeosCore/aerosol_mod.F90').read_text()
        self.assertNotIn('Spc(id_BRC_WTC)%Conc', aerosol)
        brc = (ROOT / 'GeosCore/brc_mod.F90').read_text()
        self.assertNotIn('State_Chm%Species(id_WTC)%Conc', brc)
        self.assertNotIn('State_Chm%Species(id_WTC )%Conc', brc)

    def test_restart_and_conditional_boundary_contract(self):
        text = (ROOT / 'GeosCore/hco_utilities_gc_mod.F90').read_text()
        a = text.index('Msg = BRC_WTC_Missing_Restart_Error')
        self.assertLess(a, text.index('! Check if species data is in file', a))
        self.assertIn('BRC_WTC_Total_Matches', text)
        self.assertIn('Spc(WtcOutput%Total)%Conc = 0.0_fp', text)
        meta = (ROOT / 'Headers/state_diag_mod.F90').read_text()
        self.assertIn("CASE( 'ADVWTC' )", meta)
        self.assertIn("IF ( Ind_('WTC') > State_Chm%nAdvect ) numTags = numTags + 1", meta)
        self.assertIn("IF ( N > State_Chm%nAdvect ) D = Ind_('WTC')", meta)
        bc = meta.split("diagID  = 'SpeciesBC'", 1)[1].split('IF ( RC /= GC_SUCCESS )',1)[0]
        self.assertIn("diagFlag       = 'S'", bc)
        for sim in ('fullchem','aerosol','carbon'):
            history = (ROOT / f'run/GCClassic/HISTORY.rc.templates/HISTORY.rc.{sim}').read_text()
            self.assertIn('SpeciesBC_?ADVWTC?', history)
            self.assertNotIn("'SpeciesBC_WTC", history)

    def test_no_legacy_receiver(self):
        text = (ROOT / 'GeosCore/brc_mod.F90').read_text()
        receiver = text.split(' SUBROUTINE CHEM_WTC(')[1].split(' END SUBROUTINE CHEM_WTC')[0]
        self.assertNotIn('Species(id_WTC)', receiver)
        self.assertEqual(receiver.count('BRCSOA_CONV = 0e+0_fp'), 1)
        self.assertEqual(receiver.count('NPBRC_CONV  = 0e+0_fp'), 1)
        self.assertIn('id_PWTC, id_SWTC', receiver)


def extract_native(target):
    """Extract unchanged production statements; only model interfaces are mocked."""
    text = (ROOT / 'GeosCore/brc_mod.F90').read_text()
    receiver = ' SUBROUTINE CHEM_WTC(' + text.split(' SUBROUTINE CHEM_WTC(', 1)[1].split(' END SUBROUTINE CHEM_WTC', 1)[0] + ' END SUBROUTINE CHEM_WTC\n'
    body = text[text.index('   IF ( State_Diag%Archive_BrCAbsMass'):]
    body = body.split('   IF ( Input_Opt%amIRoot ) THEN', 1)[0]
    small = re.search(r'^  REAL\(fp\), PARAMETER :: SMALLNUM = .*$', text, re.M).group(0)
    omoc = re.search(r'^  REAL\(fp\).*:: OMOC_BBOA\s*=\s*[^!\n]+', text, re.M).group(0).replace('SAVE', 'PARAMETER')
    prefix = '''MODULE Precision_Mod
 USE ISO_FORTRAN_ENV, ONLY: fp=>real64, f4=>real32
END MODULE
MODULE ErrCode_Mod
 INTEGER, PARAMETER :: GC_SUCCESS=0
END MODULE
MODULE Input_Opt_Mod
 TYPE OptInput
 LOGICAL :: LBRC=.true.
 END TYPE
END MODULE
MODULE State_Grid_Mod
 TYPE GrdState
 INTEGER :: NX=2, NY=2, NZ=2
 END TYPE
END MODULE
MODULE State_Chm_Mod
 USE ISO_FORTRAN_ENV, ONLY: real64
 TYPE SpcConc
 REAL(real64), POINTER :: Conc(:,:,:) => NULL()
 END TYPE
 TYPE AerMassType
 REAL(real64) :: BRCPI(2,2,2),BRCPO(2,2,2),NPBRC(2,2,2),WTCPI(2,2,2)
 REAL(real64) :: PWTCPI(2,2,2),SWTCPI(2,2,2),FSOAS(2,2,2),PBRC(2,2,2)
 REAL(real64) :: OCFOPOA(2,2),OCFPOA(2,2),OCPO(2,2,2),OCPISOA(2,2,2)
 REAL(real64) :: SO4_NH4_NIT(2,2,2),SNAOM(2,2,2),R_OMSNA(2,2,2),PDER(2,2,2)
 END TYPE
 TYPE ChmState
 TYPE(AerMassType) :: AerMass
 TYPE(SpcConc) :: Species(9)
 END TYPE
CONTAINS
 INTEGER FUNCTION Ind_(Name)
 CHARACTER(*), INTENT(IN) :: Name
 SELECT CASE(Name)
 CASE('WTC'); Ind_=1
 CASE('PWTC'); Ind_=4
 CASE('SWTC'); Ind_=2
 CASE DEFAULT; Ind_=0
 END SELECT
 END FUNCTION
END MODULE
MODULE State_Diag_Mod
 USE ISO_FORTRAN_ENV, ONLY: real64
 TYPE DgnState
 LOGICAL :: Archive_BrCAbsMass=.true., Archive_BrCTotMass=.true.
 LOGICAL :: Archive_BrCBleachedFrac=.true., Archive_BrCWhiteMass=.true.
 LOGICAL :: Archive_BrCPriMass=.true., Archive_BrCSecMass=.true.
 REAL(real64) :: BrCAbsMass(2,2,2), BrCTotMass(2,2,2), BrCBleachedFrac(2,2,2)
 REAL(real64) :: BrCWhiteMass(2,2,2), BrCPriMass(2,2,2), BrCSecMass(2,2,2)
 END TYPE
END MODULE
MODULE TIME_MOD
 CONTAINS
 INTEGER FUNCTION GET_TS_CHEM()
 GET_TS_CHEM=600
 END FUNCTION
END MODULE
MODULE Actual_WTC_Code
 USE ISO_FORTRAN_ENV, ONLY: fp=>real64
 USE State_Chm_Mod
 USE State_Diag_Mod
 USE State_Grid_Mod
 IMPLICIT NONE
 REAL(fp) :: BRCSOA_CONV(2,2,2), NPBRC_CONV(2,2,2)
'''
    helper = ''' SUBROUTINE FAMILY_MASS(State_Chm, State_Diag, State_Grid, optional_parents)
 TYPE(ChmState), INTENT(IN) :: State_Chm
 TYPE(DgnState), INTENT(INOUT) :: State_Diag
 TYPE(GrdState), INTENT(IN) :: State_Grid
 LOGICAL, INTENT(IN) :: optional_parents
 INTEGER :: I,J,L
 INTEGER :: id_FSOAS=3,id_BRCSOA=5,id_WTC=1,id_PWTC=4,id_SWTC=2
 INTEGER :: id_NPBRCPOA, id_PBRCPOA, id_DBRCPOA
 REAL(fp) :: AbsMass,TotMass,WhiteMass,PriMass,SecMass
 id_NPBRCPOA=6; id_PBRCPOA=7; id_DBRCPOA=8
 IF (.not.optional_parents) THEN
 id_NPBRCPOA=0; id_PBRCPOA=0; id_DBRCPOA=0
 ENDIF
'''
    aerosol = (ROOT / 'GeosCore/aerosol_mod.F90').read_text()
    aggregate = aerosol.split('          ! Zero all BrC mass fields before accumulation', 1)[1]
    aggregate = aggregate.split('          ENDIF ! LBRC', 1)[0]
    # Include actual FSOAS/PBRC/DBRC population after resetting the fields.
    aggregate += '          ENDIF\n'
    pder = aerosol.split('      ! PDER [um]', 1)[1]
    pder = pder[pder.index('      IF ( State_Chm%AerMass%SO4_NH4_NIT'):]
    pder = pder.split('#ifdef MODEL_GEOS', 1)[0]
    sizing = ''' SUBROUTINE MASS_AND_PDER(State_Chm,Input_Opt)
 USE Input_Opt_Mod, ONLY: OptInput
 TYPE(ChmState), TARGET, INTENT(INOUT) :: State_Chm
 TYPE(OptInput), INTENT(IN) :: Input_Opt
 TYPE(SpcConc), POINTER :: Spc(:)
 INTEGER :: I,J,L
 INTEGER,PARAMETER :: id_BRC_WTC=1,id_BRC_PWTC=4,id_BRC_SWTC=2,id_BRC_SOA=5,id_BRC_NPBRC=6
 INTEGER,PARAMETER :: id_FSOAS=3,id_BRC_PBRC=7,id_BRC_DBRC=8
 REAL(fp) :: AIRVOL(2,2,2)
 AIRVOL=1.0_fp
 Spc=>State_Chm%Species
 DO L=1,2
 DO J=1,2
 DO I=1,2
'''
    target.write_text(prefix + small + '\n' + omoc + '\nCONTAINS\n' + receiver + helper + body
                      + '\n END SUBROUTINE\n' + sizing + aggregate + pder
                      + '\n ENDDO\n ENDDO\n ENDDO\n END SUBROUTINE\nEND MODULE\n')


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--extract-native':
        extract_native(Path(sys.argv[2]))
    else:
        unittest.main()
