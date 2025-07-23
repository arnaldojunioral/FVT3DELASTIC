# FVT3DELASTIC

This repository provides a **free** MATLAB implementation of a three-dimensional finite-volume theory (FVT) for stress analysis in linear elastic structures. The formulation is particularly suited for cantilever beam configurations and now includes the capability to model C-shaped cross-sections, enabling the analysis of flexural, torsional, or combined flexural-torsional loading. The code offers:

* an efficient and flexible framework for conducting numerical investigations in solid mechanics;
* direct control over stress and displacement fields through a local equilibrium-based formulation; and
* easy extension to other problems involving three-dimensional structural components.

## Getting started

Save the FVT3DELASTIC.m program (17 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

**FVT3DELASTIC(nx, ny, nz)**

where **nx**, **ny**, and **nz** define the number of subvolumes along the x<sub>1</sub>, x<sub>2</sub>, and x<sub>3</sub> directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below.
<p align="center">
<img width="350" height="350" alt="image" src="https://github.com/user-attachments/assets/3d92838e-2fcb-40f7-b0da-80d891ec62d6" />
</p>

The table below summarizes the key input parameters used in the simulation, including beam geometry, material properties, loading conditions, and visualization settings.

#### Model Parameters

| Parameter       | Description                                          | Unit             |
|-----------------|------------------------------------------------------|------------------|
| `L`             | Beam length                                          | mm               |
| `H`             | Beam height                                          | mm               |
| `B`             | Beam width                                           | mm               |
| `E`             | Young's modulus (material stiffness)                 | GPa              |
| `nu`            | Poisson's ratio                                      | –                |
| `P`             | Applied load (negative indicates downward force)     | N                |
| `pb`            | `Problem: `'flexure'`, `'torsion'`, or `'torsion-flexure' | –          |
| `amp`           | Amplification factor for deformation visualization   | –                | 

<!-- ## Documentation -->

<!-- The journal article uses the FVT3DELASTIC to generate the examples presented. -->

## Example: 3D Cantilever beam

The example features a three-dimensional cantilever beam, with the analysis domain and boundary conditions illustrated in the figure below. This setup is a standard benchmark problem used to verify the functionality of the FVT3DELASTIC code.

<p align="center">
<img width="400" height="206" alt="image" src="https://github.com/user-attachments/assets/56962135-65c7-44d0-ba56-1cc5c18a9910" />
</p>

The beam dimensions are defined as:

- Length, \( L = 500 \) mm  
- Height, \( H = 100 \) mm  
- Width, \( B = 100 \) mm  

Material properties used in the analysis are:

- Young’s modulus, \( E = 150 \) GPa  
- Poisson’s ratio, \( nu = 0.3 \)  

An applied load of \( P = 2000 \) N is considered in the simulation.

Problem type: \(pb = 'flexure' \)

Deformation amplification: \(amp = 1e3 \)

#### Output

<table align="center">
  <tr>
    <td align="center" valign="top">
      <strong>Deformed structure</strong><br>
      <img width="300" height="166" alt="Deformed structure" src="https://github.com/user-attachments/assets/55c63898-61b3-42af-91ae-f9514570bcb1" />
    </td>
    <td align="center" valign="top">
      <strong>Stress fields</strong><br>
      <img width="900" height="300" alt="Stress fields" src="https://github.com/user-attachments/assets/1a826ba6-0fc4-4749-b6d4-4bfef5d56534" />
    </td>
  </tr>
</table>

## Authors

Project developed by:

* Arnaldo dos Santos Júnior  arnaldo@ctec.ufal.br
* Marcelo Victor Oliveira Araujo marcelo.vitor.o.a@gmail.com
* Romildo dos Santos Escarpini Filho romildo.escarpini@penedo.ufal.br
* Eduardo Nobre Lages enl@ctec.ufal.br
* Márcio André Araujo Cavalcante marcio.cavalcante@ceca.ufal.br
