# BrC bleaching and persistent primary emissions

Set `aerosols.carbon.bleach_scheme` in the GCClassic
`geoschem_config.yml` to select the bleaching method at runtime. The full
chemistry template sets scheme 4. The same executable can run all five schemes.
Use the same initial state, inventory, injection settings, and optics when
comparing bleaching schemes.

| Scheme | Behavior |
| --- | --- |
| 0 | Both bleaching transfers are disabled |
| 1 | Fixed 24-hour lifetime |
| 2 | Fixed 24-hour lifetime at or below 1 km AGL and zero transfer above 1 km |
| 3 | Viscosity-dependent lifetime with fixed ozone partial pressure of 3.5 × 10⁻⁸ atm |
| 4 | Viscosity-dependent lifetime with local modeled ozone |

The lifetime cap of 10⁸ seconds is finite. Schemes 3 and 4 can still
bleach at this cap. Scheme 2 explicitly bypasses the rate calculation
above its height threshold.

## Current parameter values

| Parameter | Current value | Basis and scope |
| --- | --- | --- |
| `TAU_MAX` | 10⁸ s | Numerical cap for the cold/dry regime in the viscosity formulation. Schnitzler et al. (2022) and Xie et al. (2025) motivate longer persistence in this regime. |
| `TAU_MIN` | 21,600 s | Development floor for the viscosity formulation, representing faster bleaching under warm, humid conditions. Global qualification remains open. |
| `A_PARTICLE` | 150 nm | Particle radius assumed in the Schnitzler et al. viscosity parameterization. This is a chemistry input. Aerosol optical radii are configured separately. |
| `OMOC_BBOA` | 1.8 kg OM per kg C | Fresh fire SOA conversion from FSOAS organic mass to BRCSOA carbon mass. See [BrC mass basis](brc-mass-basis.md) for the other mass conversions. |
| GFED persistent primary fraction | 0.25 | Regime-specific sensitivity setting. Long-transport absorption evidence motivates a persistent component and does not directly measure an emitted carbon fraction. |

For a GFED persistent-primary sensitivity, change the `NPBRCPOA persistent
fraction` option under extension 111 in `HEMCO_Config.rc`. The inventory
then repartitions the primary BrC source between NPBRCPOA and PBRCPOA.
Verify the combined carbon source is conserved and archive the two emission
columns. FINN, GFAS, and QFED use their own selected mappings.

Record the source and regime for each setting when comparing runs. The
minimum lifetime, particle radius, fresh OM:OC conversion, and persistent
fraction still require review before a global production prescription.
