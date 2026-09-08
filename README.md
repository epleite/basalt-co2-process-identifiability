# Basalt CO2 reactive-process identifiability

Reproducibility repository for the manuscript:

> **Breaking Reactive-Process Null Directions in Multiphysics Monitoring of Basalt Carbon Mineralization**

Author: **Emilson Pereira Leite**  
Institute of Geosciences, University of Campinas (UNICAMP), Campinas, SP, Brazil  
ORCID: **0000-0003-1691-6243**

## Scientific result

The study distinguishes signal detectability from reactive-process identifiability during basalt CO2 mineralization. Dissolution and carbonate cementation can each be detectable while producing nearly compensating seismic and DC signatures. Spectral induced polarization (SIP) adds complementary interfacial-polarization information. Target-sensitive crosswell acquisition preserves that information in nonlinear FEM data, and process-aware joint inversion reduces variance along the remaining dissolution-plus-cementation null direction by approximately a factor of five.

All distributed numerical outputs are synthetic; no proprietary time-lapse Serra Geral field data are used.

## What this GitHub repository contains

This repository is the **lightweight executable computational core** of the study. It contains:

- final reproducibility scripts for the 2-D complex-conductivity FEM and process-aware Bayesian inversion;
- the compact frozen inputs required by those scripts;
- headline frozen outputs used in the manuscript;
- scientific/repository tests;
- environment specifications, licenses, provenance, and integrity manifests.

The complete archival package also contains the manuscript PDFs, seven publication figures, source LaTeX, retained development-stage scripts, additional frozen tables, and large numerical assets. Those binary/archival materials are deposited with the versioned Zenodo record; their SHA-256 hashes are listed in `docs/BINARY_ASSETS.sha256`.

The compact `v1_0_nominal_resistivity` and `v1_1_nominal_SIP_parameters` tables in this GitHub repository contain only the rows actually consumed by the final reproduction scripts. The full corresponding frozen tables are retained in the archival package.

## Quick start

The verified environment used Python 3.13.5.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python code/run_all.py
```

Windows activation:

```powershell
.venv\Scripts\activate
```

Generated outputs are written to `reproduced/`; frozen inputs are never overwritten.

## Reproduction checks

The final executable workflows were independently rerun before repository freeze:

- v1.3 FEM headline metrics reproduced to numerical precision;
- v1.4 Bayesian headline metrics reproduced to machine precision;
- release smoke tests passed.

Exact verification details are recorded in `docs/REPRODUCTION_VERIFICATION.md`.

## Repository layout

- `code/reproduce_v1_3_fem.py` — reproduces the target-optimized crosswell complex-conductivity FEM metrics.
- `code/reproduce_v1_4.py` — reproduces the final process-aware Bayesian inversion.
- `code/run_all.py` — runs the executable GitHub core end to end.
- `data/frozen/` — compact authoritative inputs and headline outputs.
- `tests/` — structural and scientific consistency tests.
- `docs/` — provenance, claim-to-file map, AI-assistance disclosure, reproduction verification, and binary-asset hashes.

## Provenance and scope

Published literature constrains geological, petrophysical, reaction, acquisition, and SIP parameter ranges. Numerical outputs are author-generated synthetic model results. Third-party publications and source datasets are not redistributed. See `docs/PROVENANCE.md`.

The study intentionally distinguishes fully executable final analyses from earlier development stages. Frozen validated products are retained when an intermediate stage was developed interactively and no complete stand-alone original script exists. This avoids overstating code provenance.

## Licenses

- Original code: **MIT License** (`LICENSE-CODE`).
- Author-generated numerical tables and figures: **CC BY 4.0** (`DATA_LICENSE.md`).
- Third-party publications, source data, and software remain under their original terms.

## Citation and archival release

GitHub repository: `https://github.com/epleite/basalt-co2-process-identifiability`

A versioned Zenodo DOI is used for the complete archival release. The manuscript Open Research statement should cite that version DOI together with this repository.
