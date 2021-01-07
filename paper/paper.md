---
title: 'Coral: a parallel spectral solver for fluid dynamics and partial differential equations'
tags:
  - fluid dynamics
  - PDE
  - fortran
  - chebyshev
  - plane layer
authors:
  - name: Benjamin Miquel
    orcid: 0000-0001-6283-0382
    affiliation: 1
affiliations:
 - name: Université Paris-Saclay, CEA, CNRS, Service de Physique de l’Etat Condensé, 91191 Gif-sur-Yvette,France
   index: 1
date: 07 January 2021
bibliography: paper.bib
---

# Summary

`Coral` is a fast, flexible, and efficient time-stepper for solving a large class of partial differential equations, at the core of which are the Navier-Stokes equations that govern fluid motions. Written in Fortran and employing the MPI standard for parallelization, the scalability of `Coral` allows the code to leverage the resources of high-performance computing infrastructures (up to hundreds of thousands of core, see @decomp2d), while running efficiently on laptops and workstations. Equations are entered by the user in the form of a mere text file following a simple and legible syntax. No coding proficiency in Fortran is required. This flexibility makes `Coral` suitable for both students and researchers with no coding experience.


# Statement of need

Natural and industrial flows exist in numerous different flavours, including homogeneous incompressible flows, shear flow, stably or unstably stratified flows, rotating flows, and flows of an electrically conducting fluid. These flow, however, have in common that they can be modelled by sets of (quadratic) advection-diffusion equations for the velocity, and possibly for the density, the temperature, the salinity, the magnetic field, etc. Hard-coding the sets of equations corresponding to each of these flow configurations is both complex, time-consuming, and error-prone. These difficulties impede the development of new models. While `Coral` was initially motivated by the study of **Co**nvection in **Ra**pidly rotating **L**ayers, its scope has broadened and now encompasses solving homogeneous quadratic partial differential equations in a plane-layer geometry, i.e. a 3D domain with periodic boundary conditions along the two horizontal directions $x$ and $y$. Internally, `Coral` expands the variables along Fourier basis (horizontal directions) and Chebyshev polynomials (vertical direction). Transforms from physical to spectral space and domain decomposition are handled by the 2decomp&fft library [@decomp2d]. The quasi-inverse technique allows for employing an arbitrarily large numbers of Chebyshev polynomials, resulting in the ability to resolve thin boundary layers characteristic of turbulent flows without suffering from loss of accuracy. Early versions of coral have been used for studies concerning the turbulent motion of convective flows in presence of internal heat sources and sinks [@miquelPRF19, @miquelJFM20].


# Acknowledgements

The author warmly thanks Basile Gallet, Keith Julien and Nick Featherstone for discussions and encouragement during the genesis of this project. This work was granted access to the HPC resources of TGCC and CINES under the allocation 2020-A0082A10803 attributed by GENCI (Grand Equipement National de Calcul Intensif). 

# References
