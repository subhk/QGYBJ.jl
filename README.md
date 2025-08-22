QG-YBJ+ Model
==============

This is a numerical model for the two-way interaction of near-inertial waves with (Lagrangian-mean) balanced eddies. Wave evolution is governed by the YBJ+ equation (Asselin & Young 2019). The traditional quasigeostrophic equation dictates the evolution of potential vorticity, which includes the wave feedback term of Xie & Vanneste (2015). The model is pseudo-spectral in the horizontal and uses second-order finite differences to evaluate vertical and time derivatives.

**Original code by Olivier Asselin**

## Julia Implementation

This repository provides a comprehensive Julia implementation with modern features including distributed computing, advanced particle advection, and high-order interpolation schemes. The Julia package `QGYBJ.jl` offers:

### Core Numerical Methods
- **Grid and parameter setup**: `Grid`, `QGParams`, with automatic parallel domain decomposition
- **Distributed FFTs**: PencilFFTs for parallel transforms, FFTW fallback for serial
- **Spectral operators**: High-order differentiation and velocity computation
- **Elliptic solvers**: Vertical inversion with multiple boundary condition options
- **YBJ+ formulation**: Complete wave-mean flow interaction with configurable feedback
- **Time integration**: Forward Euler and Leapfrog with Robert filter and hyperdiffusion

### Advanced Particle Advection System 🚀 **NEW**
- **Unified serial/parallel**: Automatic MPI detection and domain decomposition
- **Multiple vertical levels**: 3D particle distributions with layered, random, and custom patterns
- **High-order interpolation**: Trilinear (O(h²)), Tricubic (O(h⁴)), and adaptive schemes
- **Cross-domain interpolation**: Halo exchange for accurate particle tracking across MPI boundaries
- **QG + YBJ vertical velocities**: Choose between omega equation and YBJ formulation
- **Multiple integration methods**: Euler, RK2, RK4 with automatic boundary conditions

### Wave-Mean Flow Interaction Controls
- **Configurable feedback**: Options to disable wave feedback or fix mean flow evolution
- **Physics validation**: Separate controls for wave-mean flow coupling terms
- **Parameter studies**: Systematic exploration of interaction strength effects

Getting started in Julia
------------------------

- Add dependencies in your Julia environment: `MPI`, `PencilArrays`, `PencilFFTs`, `FFTW`.
- Load the package and set up a model:

  ```julia
  using QGYBJ
  par = default_params(nx=128, ny=128, nz=64, Lx=2π, Ly=2π)
  G, S, plans, a = setup_model(; par)
  invert_q_to_psi!(S, G; a)            # ψ from q
  compute_velocities!(S, G; plans)
  S.B .= 0                             # set initial B spectrum as needed
  invert_B_to_A!(S, G, par, a)         # A and C=A_z from B (YBJ+)
  L = dealias_mask(G)
  first_projection_step!(S, G, par, plans; a, dealias_mask=L)
  Snp1 = deepcopy(S); Snm1 = deepcopy(S)
  leapfrog_step!(Snp1, S, Snm1, G, par, plans; a, dealias_mask=L)

Examples
--------

- Run `examples/demo_ybj_plus.jl` for a short YBJ+ run that writes NetCDF outputs.
- Run `examples/demo_ybj_normal.jl` to use the normal YBJ branch (`ybj_plus=false`).

MPI/PencilArrays
----------------

- Install `MPI`, `PencilArrays`, and `PencilFFTs` in your environment.
- Launch the MPI demo:
  - mpiexec -n 4 julia --project examples/demo_mpi.jl
- The code will initialize a pencil decomposition and prefer PencilFFTs for transforms if available.

Continuous Integration
----------------------

- GitHub Actions runs tests on Julia 1.9–1.11 (see `.github/workflows/ci.yml`).
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
