# [Particle API](@id api-particles)

```@meta
CurrentModule = QGYBJ
```

This page documents the particle advection API for Lagrangian tracking.

## Core Types

### ParticleConfig

Configuration for particle advection.

```@docs
ParticleConfig
```

### ParticleState

State of particle positions.

```@docs
ParticleState
```

### ParticleTracker

Tracks particles over time.

```@docs
ParticleTracker
```

## Creating Particles

### Configuration

```@docs
create_particle_config
```

### Initialization

```@docs
initialize_particles!
```

## Advection

### Main Advection Function

```@docs
advect_particles!
```

## Interpolation

### Velocity Interpolation

```@docs
interpolate_velocity_at_position
```

## I/O Functions

### Writing Trajectories

```@docs
write_particle_trajectories
```

## Interpolation Methods

The following interpolation methods are available:

```@docs
InterpolationMethod
TRILINEAR
TRICUBIC
```

## 3D Particle Distributions

For more complex 3D particle setups:

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
    domain = (G.Lx, G.Ly, 2π),  # Domain size
    method = TRILINEAR           # Interpolation method
)

# Initialize particle state
particles = ParticleState(pconfig)
initialize_particles!(particles, G)

# Create tracker for trajectories
tracker = ParticleTracker(particles, nsteps=1000)

# Advection loop
dt = par.dt
for step in 1:1000
    # Update velocities
    compute_velocities!(S, G, plans)

    # Advect particles
    advect_particles!(particles, S, G, dt)

    # Record positions
    # tracker records automatically if configured
end

# Write trajectories
write_particle_trajectories("trajectories.nc", tracker)
```
