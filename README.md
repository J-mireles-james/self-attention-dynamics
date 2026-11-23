# From Low-Dimensional Chaos to High-Dimensional Dynamics in Self-Attention Systems

**J. D. Mireles James** · November 2026

This public research repository accompanies the manuscript *From Low-Dimensional Chaos to High-Dimensional Dynamics in Self-Attention Systems* and provides the MATLAB scripts and archived numerical data used for its figures and tables.

## Research website

The project webpage is in [`docs/`](docs/). Once GitHub Pages is enabled with **Deploy from a branch → `main` → `/docs`**, it will be available at **[Self-Attention Dynamics](https://j-mireles-james.github.io/self-attention-dynamics/)**.

## Manuscript

[Read the manuscript (PDF)](paper/chaosTransformers.pdf)

## Overview

We study continuous-time projected self-attention particle models with constant query, key, and value matrices. For three particles constrained to an ellipse, the equations define a vector field on a three-dimensional torus; analytical calculations establish synchronized periodic orbits, and numerical experiments illustrate periodic, torus-like, and chaotic behavior. We then couple three-particle subsystems through block-masked, group-label, and separately normalized two-head attention. In a 60-particle experiment, 26 of the largest 40 numerically estimated Lyapunov exponents exceed 0.01.

The numerical evidence for chaos and the observed invariant structures should not be read as computer-assisted proofs of their existence; see the paper for the precise scope of the results.

## Reproducible MATLAB experiments

All code and saved data are in [`matlab/`](matlab/). The package contains one driver per **Figure 1–9** (`figure_01.m` … `figure_09.m`) and per **Table 1–2** (`table_01.m`, `table_02.m`).

With MATLAB (R2023b recommended):

1. Open the `matlab` directory in MATLAB.
2. Run a script such as `figure_07` or `table_02`.
3. By default, `recompute = false` loads the supplied archived data. To repeat the numerical integrations, set `recompute = true` in the chosen script.

Scripts locate support functions and datasets relative to their own location. Default runs display results without saving new figures or data. Long parameter scans and the 60-particle variational integration can be computationally expensive. See the [MATLAB README](matlab/README.md) for the detailed script-to-experiment mapping and test scope.

The Figure 7 torus-like experiment uses **β = 0.59**. Figure 9 reports the corrected separately normalized two-head spectra, including leading exponents close to zero at **ε = 0.001**.

## Structure

- [`paper/chaosTransformers.pdf`](paper/chaosTransformers.pdf) — current manuscript.
- [`matlab/figure_01.m` … `figure_09.m`](matlab/) — figure drivers.
- [`matlab/table_01.m`, `table_02.m`](matlab/) — numerical table drivers.
- [`matlab/functions/`](matlab/functions/) — numerical routines and utilities.
- [`matlab/data/`](matlab/data/) — archived MATLAB results and CSV tables.

## License and acknowledgment

MATLAB code is distributed under the [MIT License](LICENSE). The manuscript is provided for scholarly reading and citation; it is **not** licensed under MIT. AI assistance in the development of the manuscript and numerical code is disclosed in the manuscript and the MATLAB README.

Please cite the accompanying manuscript when using these computations. A DOI can be added after archival release.
