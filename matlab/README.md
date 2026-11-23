# From Low-Dimensional Chaos to High-Dimensional Dynamics in Self-Attention Systems

Numerical supplement by J.D. Mireles James. MIT license.

In MATLAB (R2023b recommended), open this folder and run a numbered script, e.g. `figure_07` or `table_02`. The default displays the archived numerical results in a figure window or in the command window. It writes no files. Set `recompute = true` in that script to repeat the integrations. Keep the data folder alongside the scripts: it supplies the original initial conditions, grids, and solver settings. Full scans and the 60-particle variational integration can take substantial time. Table recomputations write their checkpoints and reports to a fresh temporary directory.

| Script | Section | Experiment |
|---|---|---|
| figure_01 | 3.1 | Synchronized periodic orbit and sample trajectories |
| figure_02 | 3.2 | Forward parameter scan from zero |
| figure_03, figure_04 | 3.2 | Torus at 0.1 and chaos at 1.25 |
| figure_05 | 3.3 | Nonsynchronized periodic orbit at 0.55 |
| figure_06 | 3.3 | Two scans initialized on the periodic orbit |
| figure_07, figure_08 | 3.3 | Torus-like dynamics at 0.59 and chaos at 0.725 |
| figure_09 | 5.1 | Six-particle Lyapunov histories under two coupling strengths |
| table_01, table_02 | 5.2, 5.3 | Nine- and sixty-particle spectra and counts above 0.01 |

The upward and downward crossings in the fixed-parameter portraits use the sign of the first angular velocity. The scans display downward crossings. Figure 9 reuses the masked integration for the group-label panel by equation (15). For chaotic trajectories, compare the shape of the section and finite-time spectra across runs. Pointwise trajectories depend on solver rounding.

AI assistance was used in creating and testing this release. Package checks used Python to read the archived MATLAB arrays and verify dimensions, parameter settings, table counts and periodic closure. Execution of the rewritten MATLAB drivers and visual comparison in MATLAB remain for the author's prepublication check.
