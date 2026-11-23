# Brown carbon mass basis

This reference describes the full-chemistry BrC tracers. `MW_g` in
`species_database.yml` is the molecular weight used in tracer unit
conversions. In the carbon-mass tracers, 12.01 represents the model's carbon
mass convention.

| Species | Source and chemistry mass basis | Aerosol mass basis |
| --- | --- | --- |
| FSOAP | kg of modeled gas precursor, MW 150. GFED scales the CO mass source by a configured emitted mass ratio. Chemistry transfers this mass one-for-one to FSOAS. | Gas precursor. |
| FSOAS | kg organic matter (OM), MW 150. Produced by FSOAP conversion. Darkening removes OM mass and delivers that mass divided by 1.8 to BRCSOA as kg carbon. | Uses its OM mass directly in bin 9. |
| BRCSOA | kg carbon, MW 12.01. Produced by FSOAS darkening. Bleaching transfers carbon one-for-one to SWTC. | Uses the model's hydrophilic OC-to-OM factor in bin 6. |
| NPBRCPOA | kg carbon, MW 12.01. Receives the nonpersistent share of the routed fire OC. Bleaching transfers carbon one-for-one to PWTC. | Uses the hydrophilic OC-to-OM factor in bin 7. |
| PBRCPOA | kg carbon, MW 12.01. Receives the persistent share of the routed fire OC. | Uses the hydrophilic OC-to-OM factor in bin 10. |
| DBRCPOA | kg carbon, MW 12.01. The GFED method scales the inventory BC source by a configured mass ratio. | Uses the hydrophobic OC-to-OM factor and the dry optics carrier in bin 11. |
| PWTC, SWTC | kg carbon, MW 12.01. Transported primary and secondary products of bleaching. | Their summed mass uses the hydrophilic OC-to-OM factor in bin 8. |
| WTC | Derived sum of PWTC and SWTC on the same carbon basis and MW. | Registers the shared bin 8 output. |
| FOCPI, FOCPO | kg carbon, MW 12.01. Hold residual primary fire OC when that routing is enabled. | Use the hydrophilic and hydrophobic OC-to-OM factors, respectively, in the ordinary OC aerosol bins. |

The table gives the physical basis of the native mass fields. The standard
`SpeciesConcVV` and `SpeciesRst` fields serialize the tracers in mol/mol dry
air. Recover grid-box mass at an instantaneous sample as
`mixing_ratio * Met_AD * MW_g / 28.9644`. Use dry air mass and tracer values
from the same time and processing stage, such as native restart or
`BudgetState` values paired with the corresponding `Met_AD`. Products of
independently time-averaged mixing ratio and dry air mass can differ from
the time average of mass.

`BrCFluxFSOAP2FSOAS` records precursor mass per chemistry step, equal to OM
added to FSOAS. `BrCFluxFSOAS2BRC` records the OM removed from FSOAS per step.
Divide it by 1.8 to reconstruct carbon added to BRCSOA.
`BrCFluxBRC2WTC` and `BrCFluxNPBRC2WTC` record carbon transferred per step.
A time-averaged collection reports the mean of the sampled step amounts.
For the chemistry-only ledger, express FSOAP and FSOAS in carbon-equivalent
units by dividing both by 1.8 and sum with the four receiving and bleaching
stocks BRCSOA, NPBRCPOA, PWTC and SWTC.

The reported family mass diagnostics use kg carbon. `BrCAbsMass` includes
FSOAS / 1.8, BRCSOA, NPBRCPOA, PBRCPOA and DBRCPOA. This is a mass classification.
`BrCWhiteMass` is PWTC + SWTC and `BrCTotMass` is the sum of these two
diagnostics. The primary family is NPBRCPOA + PBRCPOA + DBRCPOA + PWTC.
The secondary family is FSOAS / 1.8 + BRCSOA + SWTC. Residual FOC and the
gas precursor are outside these aerosol-family sums. The 1.8 darkening
factor and the aerosol hydrophilic OC-to-OM factor are distinct model
assumptions. Use the latter for aerosol mass instead of converting the
family diagnostic back to OM using one uniform factor.

For the GFED primary OC ledger, sum the residual hydrophilic and hydrophobic
OC with NPBRCPOA and PBRCPOA on their emitted carbon basis. The persistent
split partitions the routed OC fraction. DBRCPOA has a separately scaled
BC-based source. Recheck the inventory-specific source factors and emission
units when comparing FINN, GFAS or QFED cases.
