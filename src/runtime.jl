"""
Convenience setup helpers to bootstrap a simulation.
"""

using SpecialFunctions: erf

"""
    setup_model(par::QGParams) -> (G, S, plans, a)

Initialize grid, state, FFT plans, and elliptic coefficient for a basic run.

# Arguments
- `par`: Model parameters (REQUIRED - use `default_params(Lx=..., Ly=..., Lz=...)`)

# Returns
Tuple of (Grid, State, Plans, a_ell)

# Note on Stratification
This function uses constant N² (from par.N²) for the elliptic coefficient.
If `par.stratification != :constant_N`, use `setup_model_with_profile()` instead
to ensure consistent physics across all operators.

# Example
```julia
par = default_params(Lx=500e3, Ly=500e3, Lz=4000.0)
G, S, plans, a = setup_model(par)
```
"""
function setup_model(par::QGParams)
    # Warn if stratification is non-constant
    if hasfield(typeof(par), :stratification) && par.stratification != :constant_N
        @warn "setup_model: Using constant N² from par.N² but par.stratification=$(par.stratification). " *
              "For non-constant stratification, use setup_model_with_profile() or the high-level " *
              "QGYBJSimulation API to ensure consistent physics in wave dispersion and vertical velocity." maxlog=1
    end

    G = init_grid(par)
    S = init_state(G)
    plans = plan_transforms!(G)
    a = a_ell_ut(par, G)
    return G, S, plans, a
end

"""
    setup_model_with_profile(par::QGParams) -> (G, S, plans, a, N2_profile)

Initialize grid, state, FFT plans, elliptic coefficient, and N² profile.

This function properly handles both constant and non-constant stratification,
returning the N² profile needed for consistent physics in wave dispersion
and vertical velocity calculations.

# Returns
Tuple of (Grid, State, Plans, a_ell, N2_profile)

# Example
```julia
par = default_params(Lx=500e3, Ly=500e3, Lz=4000.0, stratification=:skewed_gaussian)
G, S, plans, a, N2_profile = setup_model_with_profile(par)
# Use N2_profile in compute_vertical_velocity!, etc.
```
"""
function setup_model_with_profile(par::QGParams)
    G = init_grid(par)
    S = init_state(G)
    plans = plan_transforms!(G)

    # Compute N² profile based on stratification type
    if hasfield(typeof(par), :stratification) && par.stratification == :skewed_gaussian
        N2_profile = compute_N2_skewed_gaussian(G.z, par)
        a = a_ell_from_N2(N2_profile, par)
    else
        # Constant stratification
        N2_profile = fill(par.N², G.nz)
        a = a_ell_ut(par, G)
    end

    return G, S, plans, a, N2_profile
end

"""
    compute_N2_skewed_gaussian(z::Vector, par::QGParams) -> Vector

Compute N² profile using skewed Gaussian parameters from par.
"""
function compute_N2_skewed_gaussian(z::Vector{T}, par::QGParams) where T
    nz = length(z)
    N2_profile = Vector{T}(undef, nz)

    # Get skewed Gaussian parameters
    N02 = par.N02_sg
    N12 = par.N12_sg
    σ = par.sigma_sg
    z0 = par.z0_sg
    α = par.alpha_sg

    # Compute N² at each level
    for k in 1:nz
        zk = z[k] / par.Lz  # Normalize by domain depth
        N2_profile[k] = N02 + N12 * exp(-((zk - z0)/σ)^2 * (1 + erf(α * (zk - z0))))
    end

    return N2_profile
end
