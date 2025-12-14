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

## Creating and Initializing Particles

```@docs
create_particle_config
initialize_particles!
```

## Advection and Interpolation

```@docs
advect_particles!
interpolate_velocity_at_position
```

## I/O Functions

```@docs
write_particle_trajectories
```

## Interpolation Methods

```@docs
InterpolationMethod
```

## 3D Particle Distributions

```@docs
ParticleConfig3D
create_particle_config_3d
initialize_particles_3d!
```

## Usage Example

```julia
using QGYBJ

# Setup model
par = default_params(nx=64, ny=64, nz=32)
G, S, plans, a = setup_model(; par)

# Create particle configuration
pconfig = create_particle_config(
    nparticles = 1000,
    domain = (G.Lx, G.Ly, 2π)
)

# Initialize particle state
particles = ParticleState(pconfig)
initialize_particles!(particles, G)

# Advection loop
dt = par.dt
for step in 1:1000
    compute_velocities!(S, G, plans)
    advect_particles!(particles, S, G, dt)
end

# Write trajectories
write_particle_trajectories("trajectories.nc", particles)
```
