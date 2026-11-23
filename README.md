# From Low-Dimensional Chaos to High-Dimensional Dynamics in Self-Attention Systems

**J. D. Mireles James** · November 2026

Public research repository for the manuscript *From Low-Dimensional Chaos to High-Dimensional Dynamics in Self-Attention Systems* and its reproducible MATLAB numerical experiments.

## Research overview

We study continuous-time projected self-attention particle models with constant query, key, and value matrices. For three particles constrained to an ellipse, the equations define a vector field on a three-dimensional torus. Analytical calculations establish synchronized periodic orbits; numerical explorations investigate periodic, torus-like, and chaotic behavior. We also couple three-particle subsystems through block-masked, group-label, and separately normalized two-head attention. In a 60-particle computation, 26 of the largest 40 numerically estimated Lyapunov exponents exceed 0.01.

Numerical evidence should not be interpreted as a computer-assisted proof of the existence of chaotic invariant sets.

## Manuscript and MATLAB code

**The PDF and the binary MATLAB datasets have not yet been deposited in this GitHub repository.** The accompanying release package has been prepared and awaits the final file upload. When uploaded, materials will be organized as follows:

- `paper/chaosTransformers.pdf` — manuscript PDF.
- `matlab/figure_01.m` through `figure_09.m` — one driver per figure.
- `matlab/table_01.m` and `table_02.m` — one driver per table.
- `matlab/functions/` — supporting implementations.
- `matlab/data/` — archived MATLAB results and CSV tables.
- `matlab/README.md` — reproduction instructions.

The MATLAB drivers default to reading archived data (`recompute = false`). Set `recompute = true` in a driver to repeat its experiment. MATLAB R2023b is recommended. The long parameter scans and 60-particle calculation can be computationally expensive.

## License and citation

The MATLAB software is licensed under the [MIT License](LICENSE). The research manuscript is **not** covered by the MIT license. Please cite the accompanying manuscript when using these results; an archival DOI may be added later.

AI assistance in the research workflow and software development is disclosed in the manuscript and MATLAB documentation.
