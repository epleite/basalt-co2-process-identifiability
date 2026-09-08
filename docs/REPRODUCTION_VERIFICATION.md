# Reproduction verification

The executable final workflows were rerun before freezing the lightweight GitHub repository.

## Verified environment

- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- SciPy 1.17.0
- Matplotlib 3.10.8

## Full archival package

Running the original archival `code/run_all.py` completed successfully, including:

1. target-optimized v1.3 complex-conductivity FEM;
2. v1.4 process-aware Bayesian inversion;
3. release-figure integrity generation.

With the full frozen inputs, the maximum difference across the headline FEM metrics was approximately `7.1e-15`. The maximum absolute difference across the frozen v1.4 posterior headline metrics was approximately `1.1e-16`.

## Lightweight GitHub computational core

The GitHub version uses compact resistivity and SIP input tables containing only the rows consumed by the final scripts. With those compact inputs:

- the v1.3 FEM maximum difference across headline metrics was approximately `1.1e-11`;
- the v1.4 posterior maximum absolute difference remained approximately `1.1e-16`.

These differences are numerical-roundoff scale and do not change any reported scientific result.

Release smoke tests passed before freeze.
