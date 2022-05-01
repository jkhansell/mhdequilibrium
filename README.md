# MHD Equilibrium

This is a repository of project scripts for the ITCR Plasma Laboratory. It considers some magnetic field calculation routines for the MEDUSA tokamak, particle trajectories using the SCR-1 magnetic field, and diagnostics for SCR-1 magnetic field errors.

## Curve_Opt

The jupyter notebooks in this project fit the following function to BS-SOLCTRA's data output, 

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\color{white}\large R(\theta, \phi) = \sum_{i=0}^{N}\sum_{j=0}^{M} A_{ij}\cos(j\theta-iN_{fp}\phi)">
</p>

<p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=\color{white}\large Z(\theta, \phi) = \sum_{i=0}^{N}\sum_{j=0}^{M} B_{ij}\sin(j\theta-iN_{fp}\phi)">
</p>

This allows for the analytical calculation of the area and volume of the magnetic field surfaces, which is used in the [magnetic field error minimization routine](Errors/README.md). 

## Errors 


<p align="center">
   <img src="https://render.githubusercontent.com/render/math?math=\color{white}\Large e^{i\pi} = -1" >
</p>

## Medusa
## old_code
## SCR_1
## SingleParticleKabreTest
## SPM
## StructuredPIC




These projects are still on-course and updated periodically.

Developer:

- Johansell Villalobos Cubillo -- IF



