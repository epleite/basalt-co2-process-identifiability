# Numerical provenance

## Scope

All study outputs are synthetic. The model uses published literature to constrain plausible geological, petrophysical, reaction, acquisition, and electrical/SIP parameter ranges; no proprietary time-lapse Serra Geral dataset is used.

## Stage hierarchy

- **v0.1-v0.3:** reactive-state and facies-conditioned rock-physics construction.
- **v0.4-v0.6:** seismic monitorability, adversarial cancellation, exact recursive acoustic reflectivity.
- **v0.7-v0.8:** elastic PP/PS/AVA and HTI/AVAZ sensitivity.
- **v0.9:** geomechanical complement tests.
- **v1.0:** DC-resistivity complement and field-plausibility screening.
- **v1.1:** spectral induced-polarization material model.
- **v1.2:** acquisition-level crosswell screening.
- **v1.3:** 2-D complex-conductivity FEM and acquisition/inversion tests.
- **v1.4:** process-aware Bayesian joint inversion.

## Authoritative files

The manuscript's headline numbers are taken from the compact CSVs in `data/frozen/`. Large Monte Carlo scenario tables and the complete unfiltered frozen tables are retained in the archival package prepared for Zenodo.

Original executable scripts are retained when they exist as complete stand-alone programs. Several intermediate stages were developed interactively before the submission repository was frozen. For these, validated stage-level numerical tables are preserved rather than presenting later stubs as if they were the original complete executable workflows.

This distinction is particularly important for the v0.6-v0.7 wave-model development. The manuscript's exact/elastic headline tables are preserved, and the final Bayesian calculation consumes those frozen values explicitly. The repository therefore supports numerical verification of the paper without overstating code provenance.

## Public GitHub versus archival package

The public GitHub repository contains the lightweight executable computational core required to reproduce the final FEM and Bayesian headline results. The complete archival package prepared for Zenodo additionally contains manuscript binaries, publication figures, source LaTeX, retained development-stage scripts, additional frozen tables, and large numerical assets. Binary integrity is tracked in `docs/BINARY_ASSETS.sha256`.
