# [Initial Conditions](@id initial-conditions)

```@meta
CurrentModule = QGYBJ
```

This page describes how to set up initial conditions for QGYBJ.jl simulations.

## Overview

Initial conditions must be specified for:
- **Potential vorticity** ``q`` (or streamfunction ``\psi``)
- **Wave envelope** ``B`` (or wave amplitude ``A``)

## Random Initialization

### Random Flow

```julia
# Random streamfunction with specified energy level
initialize_random_flow!(state, grid;
    energy_level = 1.0,
    seed = 42
)
```

The random flow is band-limited to resolved wavenumbers with a specified spectral slope.

### Random Waves

```julia
# Random wave field
initialize_random_waves!(state, grid;
    amplitude = 0.1,
    vertical_mode = 1,  # Dominant vertical mode
    seed = 123
)
```

## Coherent Structures

### Single Vortex

```julia
# Gaussian vortex at domain center
initialize_vortex!(state, grid;
    x0 = π,
    y0 = π,
    radius = 0.5,
    intensity = 1.0,
    sign = :cyclone  # or :anticyclone
)
```

### Vortex Pair

```julia
# Dipole (vortex pair)
initialize_dipole!(state, grid;
    separation = 1.0,
    intensity = 1.0,
    angle = 0.0  # Orientation
)
```

### Jet

```julia
# Zonal jet with specified profile
initialize_jet!(state, grid;
    jet_width = 0.5,
    jet_intensity = 1.0,
    y_center = π
)
```

## Wave Packets

### Plane Wave

```julia
# Plane wave with specified wavenumber
initialize_plane_wave!(state, grid;
    kx = 4,
    ky = 0,
    amplitude = 0.1,
    vertical_mode = 1
)
```

### Gaussian Packet

```julia
# Localized wave packet
initialize_wave_packet!(state, grid;
    x0 = π,
    y0 = π,
    sigma_x = 0.5,
    sigma_y = 0.5,
    k0 = 4,     # Central wavenumber
    amplitude = 0.1
)
```

## From Spectra

### Prescribed Energy Spectrum

```julia
# Initialize with E(k) ~ k^(-3)
initialize_from_spectrum!(state, grid;
    spectrum_type = :power_law,
    exponent = -3,
    energy_level = 1.0
)
```

### Realistic Spectrum

```julia
# ECCO-like spectrum
initialize_from_spectrum!(state, grid;
    spectrum_type = :ECCO,
    energy_level = 1.0
)
```

## From Data

### From NetCDF

```julia
using NCDatasets

# Load from file
NCDataset("initial_conditions.nc") do ds
    psi_init = ds["psi"][:]
    B_init = ds["B"][:]
end

# Set state
state.psi .= psi_init
state.B .= B_init

# Compute derived quantities
invert_q_to_psi!(state, grid, params, a_ell)
compute_velocities!(state, grid, plans)
```

### Interpolation from Different Grid

```julia
# Interpolate from coarse to fine grid
psi_fine = interpolate_field(psi_coarse, grid_coarse, grid_fine)
state.psi .= psi_fine
```

## Spin-Up

For realistic simulations, start with random initialization and spin up:

```julia
# Initialize randomly
initialize_random_flow!(state, grid; energy_level=0.1)

# Spin-up phase (develop turbulence)
for step = 1:spinup_steps
    timestep!(state, grid, params, work, plans, a_ell, dt)
end

# Now add waves
initialize_random_waves!(state, grid; amplitude=0.1)

# Production run
for step = 1:nsteps
    timestep!(...)
end
```

## Verification

### Check Initial Energy

```julia
KE = flow_kinetic_energy(state.u, state.v, grid)
WE = wave_energy(state.B, state.A, grid)[2]

println("Initial KE: $KE")
println("Initial WE: $WE")
```

### Visualize

```julia
using Plots

# Surface vorticity
zeta = compute_vorticity(state.psi, grid, plans)
heatmap(zeta[:, :, end], title="Initial Surface Vorticity")

# Wave amplitude
A2 = abs2.(state.A)
heatmap(A2[:, :, end], title="Initial Wave Intensity")
```

## API Reference

Initial conditions can be set using:
- `init_random_psi!` - Initialize random streamfunction field
- Direct assignment to `state.q`, `state.B` in spectral space
- Reading from NetCDF files using `ncread_psi!`, `ncread_la!`

See the [Grid & State API](../api/grid_state.md) for state initialization functions.
