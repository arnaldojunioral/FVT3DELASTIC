# FVT3DELASTIC

This repository provides a **free** MATLAB implementation of a three-dimensional finite-volume theory (FVT) for stress analysis in linear elastic structures. The formulation is particularly suited for cantilever beam configurations and now includes the capability to model C-shaped cross-sections, enabling the analysis of flexural, torsional, or combined flexural-torsional loading. The code offers:

* an efficient and flexible framework for conducting numerical investigations in solid mechanics,
* direct control over stress and displacement fields through a local equilibrium-based formulation, and
* easy extension to other problems involving three-dimensional structural components.

## Installation

Save the FVT3DELASTIC.m program (18 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

```diff
function in blue@@ **FVT3DELASTIC(nx, ny, nz)**

## Example: 3D Cantilever beam


