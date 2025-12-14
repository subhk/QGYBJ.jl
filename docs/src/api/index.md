# [API Reference](@id api-index)

```@meta
CurrentModule = QGYBJ
```

Complete API reference for QGYBJ.jl.

## Quick Links

- [Core Types](@ref api-types): `QGParams`, `Grid`, `State`
- [Physics Functions](@ref api-physics): Inversions, operators, diagnostics

## Core Types and Setup

```@docs
QGParams
default_params
Grid
State
setup_model
```

## Time Stepping

```@docs
first_projection_step!
leapfrog_step!
```

## Elliptic Solvers

```@docs
invert_q_to_psi!
invert_B_to_A!
invert_helmholtz!
```

## Nonlinear Operators

```@docs
jacobian_spectral!
convol_waqg!
refraction_waqg!
compute_qw!
dissipation_q_nv!
```

## Transforms

```@docs
fft_forward!
fft_backward!
plan_transforms!
dealias_mask
```

## Velocities

```@docs
compute_velocities!
compute_vertical_velocity!
compute_total_velocities!
```

## Diagnostics

```@docs
wave_energy
flow_kinetic_energy
omega_eqn_rhs!
```

## YBJ Normal Mode Functions

```@docs
sumB!
compute_sigma
compute_A!
```

## I/O

```@docs
OutputManager
ncdump_psi
ncdump_la
```

## MPI Functions

Available when MPI extension is loaded:

```@docs
setup_mpi_environment
init_mpi_grid
init_mpi_state
init_mpi_workspace
plan_mpi_transforms
gather_to_root
scatter_from_root
mpi_reduce_sum
mpi_barrier
transpose_to_z_pencil!
transpose_to_xy_pencil!
```

## Module Structure

```
QGYBJ
├── Core Types
│   ├── QGParams      # Model parameters
│   ├── Grid          # Spatial grid and wavenumbers
│   └── State         # Prognostic/diagnostic fields
├── Physics
│   ├── elliptic.jl   # Tridiagonal inversions
│   ├── nonlinear.jl  # Jacobians, refraction, qw
│   ├── operators.jl  # Velocities
│   └── transforms.jl # FFT wrappers
├── Time Stepping
│   └── timestep.jl   # Leapfrog with Robert-Asselin
├── YBJ Normal Mode
│   └── ybj_normal.jl # sumB!, compute_sigma, compute_A!
├── Diagnostics
│   └── diagnostics.jl # Energy, omega equation
├── Particles
│   └── particles/    # Lagrangian tracking
└── I/O
    └── netcdf_io.jl  # NetCDF output
```

## Naming Conventions

| Suffix | Meaning | Example |
|:-------|:--------|:--------|
| `!` | In-place modification | `compute_velocities!` |
| `_spectral` | Operates in spectral space | `jacobian_spectral!` |
| `_waqg` | Wave-related | `convol_waqg!` |
| `_mpi` | MPI-enabled version | `init_mpi_grid` |
