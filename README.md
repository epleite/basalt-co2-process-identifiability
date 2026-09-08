# Basalt CO2 reactive-process identifiability

Reproducibility repository for **Breaking Reactive-Process Null Directions in Multiphysics Monitoring of Basalt Carbon Mineralization**.

**Author:** Emilson Pereira Leite, Institute of Geosciences, University of Campinas (UNICAMP), Brazil. ORCID: 0000-0003-1691-6243. Contact: emilson@unicamp.br.

## Publication status

Repository population is in progress. This is **not yet the complete v1.0.0 release**, and no Zenodo DOI has been issued for this project. The verified source package remains backed up in Google Drive. File availability in the repository must not be inferred from the planned release layout.

The initial published numerical core contains the authoritative v1.3 target-optimized FEM table, process-aware summary, and v1.4 posterior and null-direction tables. The remaining source files and binary assets are being transferred separately. The tests below check consistency of these frozen tables; they do not rerun the forward models or establish full release completeness.

## Scientific scope

The study distinguishes signal detectability from reactive-process identifiability during basalt CO2 mineralization. Dissolution and cementation can each be detectable while producing nearly compensating seismic and DC signatures. SIP supplies complementary interfacial-polarization information; its usefulness depends on target-sensitive acquisition and process-aware interpretation. All distributed numerical outputs are synthetic, not field observations.

## Check the published numerical core

```bash
python -m pip install pandas
python -m unittest discover -s tests -v
```

The complete pinned scientific environment is recorded in `requirements.txt` and `environment.yml` for the forthcoming full package.

## Licenses and provenance

Original code is MIT licensed (`LICENSE-CODE`). Author-generated numerical tables and figures are CC BY 4.0 (`DATA_LICENSE.md`). Third-party publications and source data are not redistributed or relicensed. See `docs/PROVENANCE.md`.

The `CITATION.cff` and `.zenodo.json` files prepare the software metadata; their presence does not mean that a release or DOI already exists.
