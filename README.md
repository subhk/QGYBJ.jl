QG-YBJ+ model
=============

This is a numerical model for the two-way interaction of near-inertial waves with (Lagrangian-mean) balanced eddies. Wave evolution is governed by the YBJ+ equation (Asselin & Young 2019). The traditional quasigeostrophic equation dictates the evolution of potential vorticity, which includes the wave feedback term of Xie & Vanneste (2015). The model is pseudo-spectral in the horizontal and uses second-order finite differences to evaluate vertical and time derivatives. 

Code written by Olivier Asselin


Julia rewrite (WIP)
===================

This repository includes an in-progress Julia port using PencilArrays and PencilFFTs for distributed horizontal FFTs and pencil (slab) decompositions. The Julia package lives under `src/` as `QGYBJ.jl` and currently provides:

- Grid and parameter setup (`Grid`, `QGParams`, `init_grid`, `default_params`).
- FFT planning and wrappers using PencilFFTs when available, otherwise FFTW in serial (`plan_transforms!`, `fft_forward!`, `fft_backward!`).
- Spectral operators for basic diagnostics (`compute_velocities!`).
- Elliptic inversion along the vertical for each horizontal wavenumber (`invert_q_to_psi!`).

Getting started in Julia
------------------------

- Add dependencies in your Julia environment: `MPI`, `PencilArrays`, `PencilFFTs`, `FFTW`.
- Load the package and set up a model:

  ```julia
  using QGYBJ
  par = default_params(nx=128, ny=128, nz=64, Lx=2π, Ly=2π)
  G, S, plans = setup_model(; par)
  # supply a(z), b(z) for the elliptic operator and invert:
  a = ones(G.nz); b = zeros(G.nz)
  invert_q_to_psi!(S, G; a, b)
  compute_velocities!(S, G; plans)
  ```

Porting roadmap
---------------

- Map Fortran components (`parameters*.f90`, `init.f90`, `fft.f90`, `derivatives.f90`, `elliptic.f90`, `main_waqg.f90`) to Julia modules under `src/`.
- Flesh out vertical operators (nonuniform z, boundary conditions) to match the reference model.
- Implement YBJ+ and QG time stepping (projection + leapfrog) with de-aliasing and hyperdiffusion.
- Port I/O (NetCDF) as needed using `NCDatasets.jl`.


Brief overview of files
=======================

#Essentials


parametersXXX.f90: contains all the parameters determining the simulation.

init.f90:          initialization of all basic arrays, stratification profile, initial condition for eddies and waves.

IO_ncf.f90:        all things netCDF input/output.

lcometXXX          compiling and job launching script

main_waqg.f90:     main program performing the integration



#Under the hood


elliptic.f90:      routines pertaining to inversion of q for psi, and LA for A. 

derivatives.f90:   contains various subroutines computing derivatives and nonlinear terms via the transform method.

fft.f90            all things Fourier transforms

mpi.f90            all things parallelization via MPI



#Deprecated


diagnostic.f90:    contains a bunch of various old diagnostics (obsolete)

files.f90:         initialize all text files needed (obsolete)

special.f90:       contains a couple special functions for diagnostics (obsolete)
