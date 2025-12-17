# [Particle API](@id api-particles)

```@meta
CurrentModule = QGYBJ
```

This page documents the particle advection API for Lagrangian tracking.

## Core Types

```@docs
ParticleConfig
ParticleState
ParticleTracker
```

## Particle Initialization Constructors

Simple, intuitive functions for creating particle distributions:

```@docs
particles_in_box
particles_in_circle
particles_in_grid_3d
particles_in_layers
particles_random_3d
particles_custom
```

## Initialization and Advection

```@docs
initialize_particles!
advect_particles!
interpolate_velocity_at_position
```

## I/O Functions

```@docs
write_particle_trajectories
read_particle_trajectories
write_particle_snapshot
write_particle_trajectories_by_zlevel
```

## Interpolation Methods

```@docs
InterpolationMethod
```

## 3D Particle Types

```@docs
ParticleConfig3D
ParticleDistribution
initialize_particles_3d!
```

## Quick Reference

| Constructor | Description | Example |
|:------------|:------------|:--------|
| `particles_in_box(z; ...)` | 2D box at fixed z | `particles_in_box(π/2; nx=10, ny=10)` |
| `particles_in_circle(z; ...)` | Circular disk | `particles_in_circle(1.0; radius=0.5, n=100)` |
| `particles_in_grid_3d(; ...)` | 3D grid | `particles_in_grid_3d(; nx=10, ny=10, nz=5)` |
| `particles_in_layers(zs; ...)` | Multiple z-levels | `particles_in_layers([0.5, 1.0, 1.5]; nx=10, ny=10)` |
| `particles_random_3d(n; ...)` | Random 3D | `particles_random_3d(500)` |
| `particles_custom(pos; ...)` | Custom positions | `particles_custom([(1.0,1.0,0.5), ...])` |

## Usage Example

```julia
using QGYBJ

# Setup model
par = default_params(nx=64, ny=64, nz=32)
G, S, plans, a = setup_model(; par)

# Create particle configuration (100 particles in a box at z = π/2)
pconfig = particles_in_box(π/2;
    nx=10, ny=10,
    integration_method=:rk4,
    save_interval=0.1
)

# Or use a circular distribution
pconfig = particles_in_circle(π/2; radius=1.0, n=100)

# Or multiple z-levels
pconfig = particles_in_layers([π/4, π/2, 3π/4]; nx=10, ny=10)

# Create tracker and initialize
tracker = ParticleTracker(pconfig, G)
initialize_particles!(tracker, pconfig)

# Advection loop
dt = par.dt
for step in 1:1000
    compute_velocities!(S, G, plans)
    advect_particles!(tracker, S, G, dt, step * dt)
end

# Write trajectories
write_particle_trajectories("trajectories.nc", tracker)
```
