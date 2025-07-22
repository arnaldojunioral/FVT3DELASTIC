# FVT3DELASTIC

This repository provides a **free** MATLAB implementation of a three-dimensional finite-volume theory (FVT) for stress analysis in linear elastic structures. The formulation is particularly suited for cantilever beam configurations and now includes the capability to model C-shaped cross-sections, enabling the analysis of flexural, torsional, or combined flexural-torsional loading. The code offers:

* an efficient and flexible framework for conducting numerical investigations in solid mechanics,
* direct control over stress and displacement fields through a local equilibrium-based formulation, and
* easy extension to other problems involving three-dimensional structural components.

### Installation

Save the FVT3DELASTIC.m program (18 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

$\color{blue}{\textbf{\texttt{function}}}$ **FVT3DELASTIC(nx, ny, nz)**

where **nx**, **ny**, and **nz** define the number of subvolumes along the x, y, and z-directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below.

### Documentation


### Example: 3D Cantilever beam


#### Authors
* Arnaldo dos Santos Júnior  [a link]arnaldojunioral@gmail.com
* Marcelo Victor Oliveira Araujo
* Romildo Escarpini dos Santos Filho
* Eduardo Nobre Lages
* Márcio André Araujo Cavalcante
