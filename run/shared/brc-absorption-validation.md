# BrC absorption diagnostics

The full-chemistry Aerosols collection can report species-resolved
extinction optical depth through `AODHygWL1_?HYG?` and the corresponding
configured wavelength slots WL2 and WL3. The resulting NetCDF name
includes the wavelength, for example `AODHyg550nm_BRCSOA`.
`BrCDryAOD550nm` reports dry DBRCPOA extinction. These quantities contain
scattering and absorption. `BrCAbsMass` reports carbon mass.

## Online diagnostic source

Use the online aerosol optics calculation in `aerosol_mod.F90` at the
same grid cell, layer, wavelength, humidity and radius as its extinction.
For each species and layer, absorption optical depth is
`AAOD = extinction_optical_depth * (1 - single_scattering_albedo)`.
Species-resolved AAOD should be calculated before HISTORY time averaging.
A mean extinction multiplied by one minus a separately averaged single
scattering albedo generally produces a different result.

For DBRCPOA use the dry-carrier extinction and dry single scattering
albedo. NPBRCPOA, PBRCPOA and BRCSOA use hydrophilic aerosol bins in the
current implementation. Preserve the species independently. WTC and
FSOAS are separately resolved. If a family sum includes either, state
that choice alongside the optical mode. When validating dedicated
absorption, record both the optical tables and the radius and humidity
prescriptions used. The Cloud-J lookup records and RRTMG flux differences
describe different outputs from species-resolved aerosol AAOD.

Sum layer AAOD to compare a column. For absorption Ångström exponent
at a pair of positive, distinct wavelengths, use
`-log(AAOD_short / AAOD_long) / log(lambda_short / lambda_long)`.
Sum the species and layers of interest before taking the ratio.
Mask nonpositive and below-detection-limit AAOD and record the exact
wavelength pair. For example, use 440 and 675 nm only when those are
the wavelengths of the matched observations and model diagnostics.

AERONET column AAOD includes all aerosol types and can be derived as
reported AOD times one minus reported single scattering albedo at the
same wavelength. Match its sampling and quality filters before comparing
with total model column AAOD. Model BrC-only AAOD can be reported as
an attribution. For PAS in situ absorption, divide layer AAOD by its
layer thickness in meters, then multiply by 1e6 for inverse megameters.
Sample at measurement wavelength, time, position and altitude and match
the instrument humidity and reporting conditions.

## Validation remaining

Add HISTORY output for the online species-resolved AAOD and test it in
a positive-mass fixture with qualified optical tables. Include multiple
wavelengths, bounds from zero to extinction, an independent dry-carrier
calculation, and a matched BrC-off run. Current archived extinction
fields do not contain the contemporaneous single scattering albedo
necessary to reconstruct this output exactly.
