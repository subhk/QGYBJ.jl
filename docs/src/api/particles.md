# [Particle API](@id api-particles)

```@meta
CurrentModule = QGYBJ
```

This page documents the particle advection API.

## Particle Types

### Particles

```@docs
Particles
```

Main particle structure:

```julia
struct Particles{T<:AbstractFloat}
    x::Vector{T}      # x-coordinates
    y::Vector{T}      # y-coordinates
    z::Vector{T}      # z-coordinates
    active::BitVector # Active particle mask
    properties::Dict  # User-defined properties
end
```

### ParticleTrajectories

```@docs
ParticleTrajectories
```

For recording trajectories:

```julia
struct ParticleTrajectories{T}
    x::Matrix{T}    # (nparticles, nsteps)
    y::Matrix{T}
    z::Matrix{T}
    time::Vector{T} # Time at each step
end
```

## Particle Creation

### create_particles

```@docs
create_particles
```

**Signatures:**
```julia
# Random initialization
particles = create_particles(n, grid; init_type=:random)

# From coordinates
particles = create_particles(x0, y0, z0)

# In specific region
particles = create_particles(n, grid;
    x_range = (0, π),
    y_range = (0, 2π),
    z_range = (0.5, 1.0)
)
```

### create_particle_pairs

```@docs
create_particle_pairs
```

For FSLE analysis:
```julia
pairs = create_particle_pairs(grid;
    n_pairs = 1000,
    initial_separation = 0.01
)
```

## Advection Functions

### advect_particles!

```@docs
advect_particles!
```

**Signature:**
```julia
advect_particles!(particles, u, v, w, grid, dt;
    interp = :linear,
    time_scheme = :rk4,
    bc_horizontal = :periodic,
    bc_vertical = :reflecting
)
```

**Parameters:**
- `interp`: Interpolation method (`:nearest`, `:linear`, `:cubic`, `:spectral`)
- `time_scheme`: Time integration (`:euler`, `:rk2`, `:rk4`)
- `bc_horizontal`: Horizontal boundary (`:periodic`)
- `bc_vertical`: Vertical boundary (`:reflecting`, `:absorbing`)

### advect_particles_2d!

```@docs
advect_particles_2d!
```

Horizontal-only advection:
```julia
advect_particles_2d!(particles, u, v, grid, dt)
```

### advect_particle_pairs!

```@docs
advect_particle_pairs!
```

For pair separation tracking:
```julia
advect_particle_pairs!(pairs, state, grid, dt)
```

## Interpolation

### interpolate_to_particles!

```@docs
interpolate_to_particles!
```

Sample field at particle positions:
```julia
values = zeros(nparticles)
interpolate_to_particles!(values, field, particles, grid;
    method = :linear
)
```

### interpolate_velocity

```@docs
interpolate_velocity
```

Get velocity at single point:
```julia
u_p, v_p, w_p = interpolate_velocity(x, y, z, state, grid)
```

## Trajectory Recording

### record!

```@docs
record_trajectory!
```

Store current positions:
```julia
record_trajectory!(trajectories, particles, step)
```

### ParticleTrajectories constructors

```julia
# Pre-allocate
traj = ParticleTrajectories(nparticles, nsteps)

# From particles and expected steps
traj = ParticleTrajectories(particles, nsteps)
```

## Statistics Functions

### particle_dispersion

```@docs
particle_dispersion
```

Mean squared displacement:
```math
D(t) = \langle |x(t) - x(0)|^2 \rangle
```

```julia
D = particle_dispersion(trajectories)
```

### relative_dispersion

```@docs
relative_dispersion
```

Two-particle dispersion:
```julia
D_rel = relative_dispersion(pair_trajectories)
```

### lagrangian_diffusivity

```@docs
lagrangian_diffusivity
```

Effective diffusivity:
```julia
K = lagrangian_diffusivity(trajectories)
```

### velocity_autocorrelation

```@docs
velocity_autocorrelation
```

Lagrangian velocity autocorrelation:
```julia
R = velocity_autocorrelation(trajectories, dt)
```

## FSLE Functions

### compute_fsle

```@docs
compute_fsle
```

Finite-Size Lyapunov Exponents:
```julia
lambda = compute_fsle(pairs, separations)
```

### check_separations!

```@docs
check_separations!
```

Check if pairs reached target separations:
```julia
check_separations!(pairs, target_separations)
```

## Utility Functions

### count_active

```@docs
count_active
```

Count active particles:
```julia
n = count_active(particles)
```

### remove_inactive!

```@docs
remove_inactive!
```

Remove absorbed particles:
```julia
remove_inactive!(particles)
```

### add_particles!

```@docs
add_particles!
```

Add new particles:
```julia
add_particles!(particles, new_x, new_y, new_z)
```

## Properties

### set_property!

```@docs
set_property!
```

Attach data to particles:
```julia
set_property!(particles, :temperature, temp_values)
```

### get_property

```@docs
get_property
```

Retrieve particle data:
```julia
temp = get_property(particles, :temperature)
```

## MPI Functions

### create_mpi_particles

```@docs
create_mpi_particles
```

Create distributed particles:
```julia
particles = create_mpi_particles(n_total, grid, mpi_config)
```

### exchange_particles!

```@docs
exchange_particles!
```

Transfer particles between ranks:
```julia
exchange_particles!(particles, grid, mpi_config)
```

## Full API Reference

```@docs
Particles
ParticleTrajectories
create_particles
create_particle_pairs
advect_particles!
advect_particles_2d!
advect_particle_pairs!
interpolate_to_particles!
record_trajectory!
particle_dispersion
relative_dispersion
lagrangian_diffusivity
velocity_autocorrelation
compute_fsle
count_active
```
