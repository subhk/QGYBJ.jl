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

## Quick Start

### Dependencies
Add to your Julia environment:
```julia
using Pkg
Pkg.add(["MPI", "PencilArrays", "PencilFFTs", "FFTW", "NCDatasets"])
```

### Basic Usage

#### Simple QG-YBJ Simulation
```julia
using QGYBJ

# Create simulation configuration
domain = create_domain_config(nx=64, ny=64, nz=32, Lx=2π, Ly=2π, Lz=π)
stratification = create_stratification_config(:constant_N, N0=1.0)
initial_conditions = create_initial_condition_config(:random, :random)
output = create_output_config(output_dir="./results", save_interval=0.1)

config = create_model_config(domain, stratification, initial_conditions, output,
                           total_time=2.0, dt=1e-3, Ro=0.1, Fr=0.1)

# Run simulation  
sim = setup_simulation(config)
run_simulation!(sim)
```

#### Advanced Particle Advection
```julia
# 3D particle distribution with multiple z-levels
particle_config = create_layered_distribution(
    π/2, 3π/2, π/2, 3π/2,              # Horizontal region
    [π/8, π/4, π/2, 3π/4, 7π/8],       # 5 depth layers
    10, 10,                            # 10×10 particles per layer
    use_ybj_w=true,                    # YBJ vertical velocity
    interpolation_method=TRICUBIC       # High-accuracy interpolation
)

# Initialize and run with particle tracking
tracker = ParticleTracker(particle_config, sim.grid)
initialize_particles!(tracker, particle_config)

# Simulation loop with particle advection
for step in 1:1000
    leapfrog_step!(sim.state, sim.grid, sim.params, sim.plans)
    advect_particles!(tracker, sim.state, sim.grid, sim.config.dt)
    
    if step % 100 == 0
        write_particle_snapshot("particles_$(step).nc", tracker, step * sim.config.dt)
    end
end
```

#### Parallel Execution
```bash
# Serial execution
julia examples/particle_advection_example.jl

# Parallel execution (automatic detection)
mpiexec -n 4 julia examples/particle_advection_example.jl
```

## Examples

### Core Model Examples
- **`examples/demo_ybj_plus.jl`**: YBJ+ simulation with NetCDF output
- **`examples/demo_ybj_normal.jl`**: Standard YBJ formulation
- **`examples/test_ybj_vertical_velocity.jl`**: Compare QG vs YBJ vertical velocities

### Particle Advection Examples 🆕
- **`examples/particle_advection_example.jl`**: Comprehensive particle tracking demonstration
- **`examples/3d_particle_distribution_example.jl`**: Multiple z-level and 3D distribution patterns
- **`examples/interpolation_comparison_example.jl`**: Performance comparison of interpolation methods

### Advanced Features
- **Wave-mean flow interaction controls**: Disable feedback or fix mean flow
- **High-order interpolation schemes**: Tricubic for improved accuracy
- **Parallel particle migration**: Seamless cross-domain particle tracking

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
