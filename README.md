# MHD Equilibrium

This is a repository of project scripts for the ITCR Plasma Laboratory. It considers some magnetic field calculation routines for the MEDUSA tokamak, particle trajectories using the SCR-1 magnetic field, and diagnostics for SCR-1 magnetic field errors.

## Curve_Opt

The jupyter notebooks in this project fit the following function to BS-SOLCTRA's data output, 


$$R(\theta, \phi) = \sum_{i=0}^{N}\sum_{j=0}^{M} A_{ij}\cos(j\theta-iN_{fp}\phi)$$

$$Z(\theta, \phi) = \sum_{i=0}^{N}\sum_{j=0}^{M} B_{ij}\sin(j\theta-iN_{fp}\phi)$$

This allows for the analytical calculation of the area and volume of the magnetic field surfaces, which is used in the [magnetic field error minimization routine](Errors/README.md). 

## Errors 

In this module, a magnetic field coil optimization was carried out for the SCR-1 stellarator's modular coil design. The code was designed to optimize the fourier coefficients of the 12 modular coils using a genetic algorithm based on


## Medusa
## old_code
## SCR_1
## SingleParticleKabreTest
## SPM
## StructuredPIC

These projects are still on-course and updated periodically.

Developer:

- Johansell Villalobos Cubillo -- IF
