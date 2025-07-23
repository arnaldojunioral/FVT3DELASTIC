# FVT3DELASTIC

This repository provides a **free** MATLAB implementation of a three-dimensional finite-volume theory (FVT) for stress analysis in linear elastic structures. The formulation is particularly suited for cantilever beam configurations and now includes the capability to model C-shaped cross-sections, enabling the analysis of flexural, torsional, or combined flexural-torsional loading. The code offers:

* an efficient and flexible framework for conducting numerical investigations in solid mechanics;
* direct control over stress and displacement fields through a local equilibrium-based formulation; and
* easy extension to other problems involving three-dimensional structural components.

## Getting started

Save the FVT3DELASTIC.m program (17 kB) and launch MATLAB in the same directory. The program can be executed with the following command:

**FVT3DELASTIC(nx, ny, nz)**

where **nx**, **ny**, and **nz** define the number of subvolumes along the x<sub>1</sub>, x<sub>2</sub>, and x<sub>3</sub> directions, respectively. These parameters specify the discretization of the three-dimensional domain, as illustrated in the figure below.
<!-- <p align="center">
<img width="350" height="350" alt="image" src="https://github.com/user-attachments/assets/3d92838e-2fcb-40f7-b0da-80d891ec62d6" />
</p> -->
<p align="center">
<img width="510" height="280" alt="image" src="https://github.com/user-attachments/assets/f6484d68-3d51-4040-b83f-3719f702d792" />
</p>

The table below summarizes the key input parameters used in the simulation, including beam geometry, material properties, loading conditions, and visualization settings.

#### Model Parameters

| Parameter       | Description                                          | Unit             |
|-----------------|------------------------------------------------------|------------------|
| `L`             | Beam length                                          | mm               |
| `H`             | Beam height                                          | mm               |
| `B`             | Beam width                                           | mm               |
| `E`             | Young's modulus (material stiffness)                 | MPa              |
| `nu`            | Poisson's ratio                                      | –                |
| `P`             | Applied load (negative indicates downward force)     | N                |
| `pb`            | `Problem: `'flexure'`, `'torsion'`, or `'torsion-flexure' | –          |
| `amp`           | Amplification factor for deformation visualization   | –                | 

<!-- ## Documentation -->

<!-- The journal article uses the FVT3DELASTIC to generate the examples presented. -->

## Example: 3D Cantilever beam

This example presents a standard benchmark problem involving a three-dimensional cantilever beam. The analysis domain and boundary conditions, shown in the figure below, are used to verify the functionality of the FVT3DELASTIC code.

<p align="center">
<img width="400" height="206" alt="image" src="https://github.com/user-attachments/assets/56962135-65c7-44d0-ba56-1cc5c18a9910" />
</p>

| Parameter               | Symbol | Value     | Unit    |
|-------------------------|--------|-----------|---------|
| **Beam Length**         | L      | 500       | mm      |
| **Beam Height**         | H      | 100       | mm      |
| **Beam Width**          | B      | 100       | mm      |
| **Young’s Modulus**     | E      | 150e3     | MPa     |
| **Poisson’s Ratio**     | ν      | 0.3       | —       |
| **Applied Load**        | P      | 2000      | N       |
| **Problem Type**        | pb     | 'flexure' | —       |
| **Deformation Amplification** | amp    | 1e3       | —       |

#### Output

<table align="center">
  <tr>
    <td align="center" valign="top">
      <strong>Deformed structure</strong><br>
      <img width="375" height="207" alt="Deformed" src="https://github.com/user-attachments/assets/55c63898-61b3-42af-91ae-f9514570bcb1" />
    </td>
    <td align="center" valign="top">
      <strong>Stress fields</strong><br>
      <img width="900" height="300" alt="Stress fields" src="https://github.com/user-attachments/assets/1a826ba6-0fc4-4749-b6d4-4bfef5d56534" />
    </td>
  </tr>
</table>

## Error Reporting

We strive to ensure that the implementation of the finite-volume theory is accurate, efficient, and well-documented. However, if you encounter unexpected behavior, inconsistencies, or potential bugs in the code, we welcome your feedback.

Please feel free to:

- **Open an issue** on the [GitHub repository](https://github.com/arnaldojunioral/FVT3DELASTIC/issues) describing the problem in detail.
- Or **contact us directly** via email at [arnaldo@ctec.ufal.br](mailto:arnaldo@ctec.ufal.br).

Your contributions help improve the reliability and usability of this project for the research community.

## Authors

Project developed by:

* Arnaldo dos Santos Júnior  arnaldo@ctec.ufal.br
* Marcelo Victor Oliveira Araujo marcelo.vitor.o.a@gmail.com
* Romildo dos Santos Escarpini Filho romildo.escarpini@penedo.ufal.br
* Eduardo Nobre Lages enl@ctec.ufal.br
* Márcio André Araujo Cavalcante marcio.cavalcante@ceca.ufal.br

## References

The following table summarizes the five relevant references supporting the development of the proposed three-dimensional finite-volume theory. These works were selected based on their conceptual alignment with the present formulation, their methodological contributions, and their scientific impact.

| Rank | Reference                                                                                          | Relevance to the Study      | Scientific Impact | Justification                                                                                                                                                  |
|------|----------------------------------------------------------------------------------------------------|-----------------------------|-------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1    | Cavalcante, M.A.A., Pindera, M.-J. (2012a,b). *Generalized finite-volume theory for elastic analysis in solid mechanics* (Parts I & II). J. Appl. Mech. | ⭐⭐⭐⭐⭐ | High              | Serves as the primary foundation for the present formulation. These papers established the generalized FVT framework and addressed key numerical challenges.     |
| 2    | Bansal, Y., Pindera, M.-J. (2003, 2005). Reformulated higher-order theory for periodic multiphase materials. | ⭐⭐⭐⭐⭐ | High              | Provides the theoretical basis for the discretization and stiffness matrix strategies adopted in this work, particularly for heterogeneous and periodic media.    |
| 3    | Zhong, Y., Bansal, Y., Pindera, M.-J. (2004). *Efficient reformulation of the thermal higher-order theory for FGMs.* | ⭐⭐⭐⭐☆ | Medium            | Offers important methodological insights relevant to extending the FVT to three-dimensional domains with variable material properties.                          |
| 4    | Cardiff, P., Demirdžić, I. (2021). *Thirty years of the finite volume method for solid mechanics.* Arch. Comput. Methods Eng. | ⭐⭐⭐⭐☆ | Very High         | A comprehensive review that situates the present study within the broader evolution of finite-volume methods in solid mechanics.                                 |
| 5    | Araujo, M.V.O., Lages, E.N., Cavalcante, M.A.A. (2020a, 2021). Applications of FVT in topology optimization and structural energy analysis. | ⭐⭐⭐⭐☆ | Medium            | Demonstrates the practical applicability of FVT to topology optimization and validates the versatility of the formulation for linear elastic structures.         |
