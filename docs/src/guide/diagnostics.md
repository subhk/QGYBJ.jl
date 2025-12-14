# [Diagnostics](@id diagnostics)

```@meta
CurrentModule = QGYBJ
```

This page describes diagnostic quantities and analysis tools in QGYBJ.jl.

## Energy Diagnostics

### Flow Kinetic Energy

The kinetic energy of the balanced flow:

```math
KE = \frac{1}{2}\int (u^2 + v^2) \, dV
```

```julia
KE = flow_kinetic_energy(state.u, state.v, grid)
```

### Flow Potential Energy

Available potential energy:

```math
PE = \frac{1}{2}\int \frac{f_0^2}{N^2}\left(\frac{\partial\psi}{\partial z}\right)^2 dV
```

```julia
PE = flow_potential_energy(state.psi, grid, params)
```

### Total Flow Energy

```julia
E_flow = flow_total_energy(state, grid, params)
# Equivalent to: KE + PE
```

### Wave Energy

```math
E_{wave} = \frac{1}{2}\int |A|^2 \, dV
```

```julia
E_B, E_A = wave_energy(state.B, state.A, grid)
```

!!! note
    ``E_B`` and ``E_A`` differ because the ``L^+`` operator is not unitary.
    ``E_A`` is the physical wave energy.

## Enstrophy

### Relative Enstrophy

```math
Z_{\zeta} = \frac{1}{2}\int \zeta^2 \, dV
```

where ``\zeta = \nabla^2\psi`` is relative vorticity.

```julia
Z_rel = relative_enstrophy(state.psi, grid)
```

### Potential Enstrophy

```math
Z_q = \frac{1}{2}\int q^2 \, dV
```

```julia
Z_pot = potential_enstrophy(state.q, grid)
```

## Spectral Diagnostics

### Energy Spectra

Compute 1D energy spectrum:

```julia
# Horizontal wavenumber spectrum
k, E_k = horizontal_energy_spectrum(state.psi, grid)

# Plot
using Plots
plot(k, E_k, xscale=:log10, yscale=:log10,
    xlabel="k", ylabel="E(k)", label="Energy")
plot!(k, 0.1 .* k.^(-3), label="k⁻³", linestyle=:dash)
```

### Vertical Spectra

```julia
# Vertical mode decomposition
m, E_m = vertical_energy_spectrum(state.psi, grid, params)
```

### 2D Spectrum

```julia
# 2D spectrum E(kx, ky)
E_2d = compute_2d_spectrum(state.psi, grid)
heatmap(log10.(E_2d), xlabel="kx", ylabel="ky")
```

## Vorticity Statistics

### Vorticity Field

```julia
# Compute vorticity
zeta = compute_vorticity(state.psi, grid, plans)
```

### Statistics

```julia
# Mean, variance, skewness, kurtosis
stats = vorticity_statistics(zeta)

println("Mean: $(stats.mean)")
println("Variance: $(stats.variance)")
println("Skewness: $(stats.skewness)")
println("Kurtosis: $(stats.kurtosis)")
```

### PDF

```julia
# Probability density function
bins, pdf = vorticity_pdf(zeta; nbins=100)

plot(bins, pdf, xlabel="ζ", ylabel="P(ζ)")
```

## Wave Diagnostics

### Wave Amplitude

```julia
# |A|² field
A2 = abs2.(state.A)

# Maximum amplitude
A_max = maximum(sqrt.(A2))

# Volume-averaged amplitude
A_rms = sqrt(mean(A2))
```

### Wave Fluxes

Energy flux in vertical:

```julia
# Vertical wave energy flux
Fz = compute_wave_energy_flux_z(state.A, grid, params)
```

### Wave Polarization

Horizontal wave velocity components:

```julia
u_wave, v_wave = compute_wave_velocities(state.A, grid, plans)

# Polarization (should be ~1 for pure inertial oscillation)
polarization = mean(abs2.(u_wave) ./ abs2.(v_wave))
```

## Omega Equation

Compute ageostrophic vertical velocity:

```julia
# Solve omega equation
w = compute_omega(state.psi, grid, params, plans)

# Maximum vertical velocity
w_max = maximum(abs.(w))
```

The omega equation:
```math
\nabla^2 w + \frac{N^2}{f_0^2}\frac{\partial^2 w}{\partial z^2} = 2J\left(\frac{\partial\psi}{\partial z}, \nabla^2\psi\right)
```

## Eddy Identification

### Vortex Cores

Identify coherent vortices:

```julia
# Okubo-Weiss parameter
OW = compute_okubo_weiss(state.psi, grid, plans)

# Vortex cores: OW < 0
vortex_mask = OW .< -threshold
```

### Eddy Census

```julia
# Count and characterize eddies
eddies = identify_eddies(state.psi, grid;
    threshold = 0.1,
    min_radius = 3  # grid points
)

println("Found $(length(eddies)) eddies")
for e in eddies
    println("  Radius: $(e.radius), Intensity: $(e.intensity)")
end
```

## Time Series Analysis

### Recording Diagnostics

```julia
# Initialize time series storage
ts = DiagnosticsTimeSeries()

for step = 1:nsteps
    timestep!(state, ...)

    # Record diagnostics
    push!(ts.time, step * dt)
    push!(ts.KE, flow_kinetic_energy(state.u, state.v, grid))
    push!(ts.PE, flow_potential_energy(state.psi, grid, params))
    push!(ts.WE, wave_energy(state.B, state.A, grid)[2])
end
```

### Energy Conservation Check

```julia
# Total energy should be conserved (inviscid)
E_total = ts.KE .+ ts.PE .+ ts.WE

# Check conservation
dE = (E_total[end] - E_total[1]) / E_total[1]
println("Relative energy change: $dE")

if abs(dE) > 1e-6
    @warn "Energy not well conserved!"
end
```

### Growth Rates

```julia
# Compute growth rate from time series
using Statistics

# Linear fit in log space
log_E = log.(ts.KE)
t = ts.time

growth_rate = (log_E[end] - log_E[1]) / (t[end] - t[1])
println("KE growth rate: $growth_rate")
```

## Budget Analysis

### Energy Budget

```julia
# Compute all energy budget terms
budget = energy_budget(state, state_old, grid, params, dt)

println("dKE/dt = $(budget.dKE_dt)")
println("Advection: $(budget.advection)")
println("Dissipation: $(budget.dissipation)")
println("Wave feedback: $(budget.wave_feedback)")
println("Residual: $(budget.residual)")
```

### Enstrophy Budget

```julia
budget_Z = enstrophy_budget(state, grid, params)
```

## Spatial Averaging

### Horizontal Mean

```julia
# Average over horizontal domain
psi_z = horizontal_mean(state.psi, grid)  # Function of z only
```

### Vertical Mean

```julia
# Depth-averaged field
psi_xy = vertical_mean(state.psi, grid)  # Function of x, y only
```

### Zonal Mean

```julia
# Average over x (zonal direction)
psi_yz = zonal_mean(state.psi, grid)
```

## Correlation Functions

### Autocorrelation

```julia
# Spatial autocorrelation of vorticity
lag, R = spatial_autocorrelation(zeta, grid; direction=:x)

# Integral length scale
L_int = sum(R .* diff([0; lag])) / R[1]
```

### Cross-Correlation

```julia
# Correlation between wave amplitude and vorticity
R_wave_zeta = cross_correlation(abs.(state.A), zeta, grid)
```

## Diagnostic Output

### Comprehensive Diagnostics

```julia
function compute_all_diagnostics(state, grid, params, plans)
    diag = Dict{String, Any}()

    # Energy
    diag["KE"] = flow_kinetic_energy(state.u, state.v, grid)
    diag["PE"] = flow_potential_energy(state.psi, grid, params)
    diag["WE_B"], diag["WE_A"] = wave_energy(state.B, state.A, grid)

    # Enstrophy
    diag["Z_rel"] = relative_enstrophy(state.psi, grid)
    diag["Z_pot"] = potential_enstrophy(state.q, grid)

    # Extrema
    zeta = compute_vorticity(state.psi, grid, plans)
    diag["zeta_max"] = maximum(zeta)
    diag["zeta_min"] = minimum(zeta)
    diag["A_max"] = maximum(abs.(state.A))

    return diag
end
```

### Formatted Output

```julia
function print_diagnostics(diag, step, time)
    println("=" ^ 60)
    println("Step: $step, Time: $(round(time, digits=4))")
    println("-" ^ 60)
    println("  KE = $(round(diag["KE"], sigdigits=6))")
    println("  PE = $(round(diag["PE"], sigdigits=6))")
    println("  WE = $(round(diag["WE_A"], sigdigits=6))")
    println("  ζ_max = $(round(diag["zeta_max"], sigdigits=4))")
    println("  |A|_max = $(round(diag["A_max"], sigdigits=4))")
    println("=" ^ 60)
end
```

## API Reference

See the [Physics API Reference](../api/physics.md) for complete documentation of diagnostic functions:
- `wave_energy` - Wave energy diagnostics
- `omega_eqn_rhs!` - Omega equation RHS computation

Additional diagnostic utilities are available through the model interface.
