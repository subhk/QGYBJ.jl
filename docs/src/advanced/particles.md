# [Particle Advection](@id particles)

```@meta
CurrentModule = QGYBJ
```

This page describes Lagrangian particle tracking in QGYBJ.jl.

## Overview

Particle tracking allows you to:

- Follow **fluid parcels** as they move with the flow
- Compute **Lagrangian statistics** (dispersion, diffusivity)
- Track **tracer concentrations** along trajectories
- Study **mixing** and **transport**

## Quick Start

### Basic Setup

```julia
using QGYBJ

# Create grid and state (from simulation)
grid = Grid(nx=128, ny=128, nz=64)
state = create_state(grid)
# ... initialize flow ...

# Create particles
np = 1000  # Number of particles
particles = create_particles(np, grid;
    init_type = :random,  # Random initial positions
    vertical_range = (0.8, 1.0)  # Upper 20% of domain
)
```

### Advection

```julia
# In time loop
for step = 1:nsteps
    timestep!(state, ...)

    # Advect particles with flow
    advect_particles!(particles, state.u, state.v, state.w, grid, dt)
end
```

## Particle Initialization

### Random Positions

```julia
particles = create_particles(np, grid;
    init_type = :random,
    x_range = (0, 2π),
    y_range = (0, 2π),
    z_range = (0.5, 1.0)
)
```

### Regular Grid

```julia
particles = create_particles(np, grid;
    init_type = :regular,
    nx_part = 10,
    ny_part = 10,
    nz_part = 10
)
# Creates 1000 particles on regular 10×10×10 grid
```

### From Coordinates

```julia
# Custom initial positions
x0 = rand(np) .* 2π
y0 = rand(np) .* 2π
z0 = 0.9 .* ones(np)  # Near surface

particles = create_particles(x0, y0, z0)
```

### Within Eddies

```julia
# Initialize particles inside anticyclones
particles = create_particles(np, grid;
    init_type = :vortex,
    vorticity_sign = :negative,  # Anticyclones
    threshold = -0.5
)
```

## Advection Methods

### Interpolation Schemes

| Scheme | Order | Speed | Smoothness |
|:-------|:------|:------|:-----------|
| `:nearest` | 0 | Fastest | Discontinuous |
| `:linear` | 1 | Fast | C⁰ |
| `:cubic` | 3 | Medium | C¹ |
| `:spectral` | N | Slow | Spectral |

```julia
# Set interpolation method
advect_particles!(particles, u, v, w, grid, dt;
    interp = :cubic
)
```

### Time Integration

| Method | Order | Evaluations/step |
|:-------|:------|:-----------------|
| `:euler` | 1 | 1 |
| `:rk2` | 2 | 2 |
| `:rk4` | 4 | 4 |

```julia
advect_particles!(particles, u, v, w, grid, dt;
    time_scheme = :rk4
)
```

### Horizontal-Only Advection

For 2D studies (or when w is negligible):

```julia
# Only use horizontal velocities
advect_particles_2d!(particles, u, v, grid, dt)
```

## Particle Data

### Accessing Positions

```julia
# Current positions
x = particles.x  # Array of x-coordinates
y = particles.y  # Array of y-coordinates
z = particles.z  # Array of z-coordinates

# As matrix
positions = [particles.x particles.y particles.z]
```

### Particle Properties

You can attach properties to particles:

```julia
# Add temperature tracer
particles.T = zeros(np)

# Initialize from field
interpolate_to_particles!(particles.T, temperature_field, particles, grid)

# Update along trajectory (advection-diffusion)
advect_tracer!(particles.T, particles, state, grid, dt; kappa=1e-4)
```

### Trajectory History

```julia
# Record trajectories
traj = ParticleTrajectories(particles, nsteps)

for step = 1:nsteps
    timestep!(state, ...)
    advect_particles!(particles, ...)

    # Store current positions
    record!(traj, particles, step)
end

# Access history
x_history = traj.x  # (np, nsteps)
```

## Boundary Conditions

### Periodic

Default for horizontal directions:

```julia
# Particles wrap around domain edges
advect_particles!(particles, u, v, w, grid, dt;
    bc_horizontal = :periodic
)
```

### Reflecting Vertical

```julia
# Particles reflect off top/bottom
advect_particles!(particles, u, v, w, grid, dt;
    bc_vertical = :reflecting
)
```

### Absorbing

```julia
# Particles removed when hitting boundary
advect_particles!(particles, u, v, w, grid, dt;
    bc_vertical = :absorbing
)

# Check for removed particles
n_active = count(particles.active)
```

## Lagrangian Statistics

### Dispersion

```julia
# Single particle dispersion
D = particle_dispersion(traj)  # <|x(t) - x(0)|²>

# Plot dispersion
using Plots
plot(traj.time, D, xlabel="Time", ylabel="Dispersion")
```

### Two-Particle Statistics

```julia
# Relative dispersion
D_rel = relative_dispersion(traj)

# Richardson scaling: D_rel ~ t³
```

### Diffusivity

```julia
# Effective diffusivity
K_eff = lagrangian_diffusivity(traj)

# Or from velocity autocorrelation
K = velocity_autocorrelation_integral(traj)
```

### Velocity Statistics

```julia
# Lagrangian velocity
u_lag = lagrangian_velocity(particles, state.u, grid)

# Autocorrelation
R_u = velocity_autocorrelation(traj)

# Lagrangian integral time scale
T_L = integral_time_scale(R_u, dt)
```

## FSLE Analysis

Finite-Size Lyapunov Exponents for mixing:

```julia
# Initialize particle pairs
pairs = create_particle_pairs(grid;
    n_pairs = 500,
    initial_separation = 1e-3
)

# Track until specified separations
separations = [0.01, 0.1, 1.0]

for step = 1:nsteps
    timestep!(state, ...)
    advect_particle_pairs!(pairs, state, grid, dt)
    check_separations!(pairs, separations)
end

# Compute FSLE
lambda = compute_fsle(pairs, separations)
```

## Visualization

### Particle Positions

```julia
using Plots

# 2D horizontal slice
scatter(particles.x, particles.y,
    markersize=1, alpha=0.5,
    xlabel="x", ylabel="y",
    title="Particle Distribution"
)
```

### Trajectories

```julia
# Plot trajectories
p = plot()
for i in 1:min(100, np)  # First 100 particles
    plot!(p, traj.x[i, :], traj.y[i, :],
        alpha=0.3, legend=false)
end
display(p)
```

### Animation

```julia
anim = @animate for t in 1:10:size(traj.x, 2)
    scatter(traj.x[:, t], traj.y[:, t],
        markersize=1,
        xlim=(0, 2π), ylim=(0, 2π),
        title="t = $(t*dt)"
    )
end

gif(anim, "particles.gif", fps=20)
```

### Colored by Property

```julia
# Color by depth
scatter(particles.x, particles.y,
    zcolor=particles.z,
    colorbar=true,
    colorbar_title="Depth"
)
```

## MPI Parallelization

Particles work with MPI decomposition:

```julia
# Particles are distributed across ranks
particles = create_mpi_particles(np_total, grid, mpi_config)

# Advection handles communication
advect_particles_mpi!(particles, state, grid, mpi_config, dt)
```

Particles crossing domain boundaries are automatically transferred.

## Performance Tips

### Large Particle Counts

```julia
# Use GPU for many particles
using CUDA

particles_gpu = cu(particles)
advect_particles_gpu!(particles_gpu, state_gpu, grid_gpu, dt)
```

### Spatial Indexing

For efficient interpolation with many particles:

```julia
# Build spatial index
index = build_particle_index(particles, grid)

# Interpolation uses index
interpolate_to_particles!(field_p, field, particles, grid;
    index=index
)
```

### Chunked Processing

For very large particle counts:

```julia
# Process in chunks to limit memory
chunk_size = 10000
for i in 1:chunk_size:np
    chunk = i:min(i+chunk_size-1, np)
    advect_particles!(particles[chunk], state, grid, dt)
end
```

## API Reference

See the [Particle API Reference](../api/particles.md) for complete documentation of:
- `ParticleConfig`, `ParticleState`, `ParticleTracker`
- `create_particle_config`, `initialize_particles!`
- `advect_particles!`, `interpolate_velocity_at_position`
- `write_particle_trajectories`
