# Claim-to-file map

This map identifies the primary numerical evidence for the manuscript's main claims. The GitHub repository contains the executable final computational core; the complete Zenodo archive retains all frozen tables, figures, manuscript binaries, source LaTeX, and development-stage assets.

| Manuscript claim | Primary numerical evidence | Availability |
|---|---|---|
| Physical cancellation survives exact acoustic recursion/internal multiples | `data/frozen/basalt_co2_v0_6_reflectivity_energy.csv` | GitHub + archive |
| PP/PS do not remove dissolution-cementation collinearity | `data/frozen/basalt_co2_v0_7_basis_similarity.csv` | GitHub + archive |
| DC response enters the final process inversion | `data/frozen/basalt_co2_v1_0_nominal_resistivity.csv` | compact GitHub input + full archive |
| SIP material response enters FEM and final inversion | `data/frozen/basalt_co2_v1_1_nominal_SIP_parameters.csv` | compact GitHub input + full archive |
| FEM-optimized crosswell SIP survives phase noise | `data/frozen/basalt_co2_v1_3_FEM_target_optimized.csv` and `code/reproduce_v1_3_fem.py` | GitHub + archive |
| Conventional property inversion can lose process separation | `basalt_co2_v1_3_FEM_Jacobian_inversion_MC.csv` | full archive |
| SIP reduces D+C null-direction variance by ~5x | `data/frozen/basalt_co2_v1_4_null_direction_variance.csv` and `code/reproduce_v1_4.py` | GitHub + archive |
| SIP approximately halves D/C posterior uncertainty | `data/frozen/basalt_co2_v1_4_headline_posterior.csv` | GitHub + archive |
| SIP benefit survives model discrepancy/noise stress tests | `basalt_co2_v1_4_stress_summary.csv`, `basalt_co2_v1_4_SIP_noise_sensitivity.csv` | full archive |

Binary manuscript/figure integrity is recorded in `docs/BINARY_ASSETS.sha256`.
