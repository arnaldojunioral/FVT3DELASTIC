# FVT3DELASTIC

This repository provides a **free** MATLAB implementation of a three-dimensional finite-volume theory (FVT) for stress analysis in linear elastic structures. The formulation is particularly suited for cantilever beam configurations and now includes the capability to model C-shaped cross-sections, enabling the analysis of flexural, torsional, or combined flexural-torsional loading. The code offers:

* an efficient and flexible framework for conducting numerical investigations in solid mechanics,
* direct control over stress and displacement fields through a local equilibrium-based formulation, and
* easy extension to other problems involving three-dimensional structural components.

### Installation

Save the FVT3DELASTIC.m program (17 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

$\color{blue}{\textbf{\texttt{function}}}$ **FVT3DELASTIC(nx, ny, nz)**

where **nx**, **ny**, and **nz** define the number of subvolumes along the x, y, and z-directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below.

### Documentation

The journal article uses the FVT3DELASTIC to generate the examples presented.

### Example: 3D Cantilever beam

![Figure 7](https://github.com/user-attachments/assets/1cea754d-6dd5-469f-8e68-f7c2e77087ec)


<img width="2048" height="780" alt="image" src="https://github.com/user-attachments/assets/1a826ba6-0fc4-4749-b6d4-4bfef5d56534" />


#### Authors

Project developed by:

* Arnaldo dos Santos Júnior  arnaldo@ctec.ufal.br
* Marcelo Victor Oliveira Araujo marcelo.vitor.o.a@gmail.com
* Romildo dos Santos Escarpini Filho romildo.escarpini@penedo.ufal.br
* Eduardo Nobre Lages enl@ctec.ufal.br
* Márcio André Araujo Cavalcante marcio.cavalcante@ceca.ufal.br
