# [Time Stepping](@id api-timestepping)

```@meta
CurrentModule = QGYBJ
```

This page documents the time integration functions.

## Main Time Stepping Scheme

QGYBJ.jl uses a **Leapfrog scheme with Robert-Asselin filter** for time integration. This provides second-order accuracy while maintaining stability through computational mode damping.

### Overview

The time stepping consists of two functions:
1. `first_projection_step!` - Forward Euler initialization
2. `leapfrog_step!` - Main leapfrog integration with Robert-Asselin filter

## Forward Euler Projection Step

```@docs
first_projection_step!
```

**Purpose:** Initialize the leapfrog scheme by providing values at times n and n-1.

**Algorithm:**
1. Compute tendencies at time n (advection, refraction, diffusion)
2. Apply physics switches (linear, inviscid, etc.)
3. Forward Euler update with integrating factors
4. Wave feedback (optional)
5. Diagnostic inversions (q → ψ → u, v)

**Usage:**
```julia
# Serial mode
first_projection_step!(state, grid, params, plans; a=a_ell, dealias_mask=L)

# Parallel mode (2D decomposition)
first_projection_step!(state, grid, params, plans; a=a_ell, dealias_mask=L, workspace=workspace)
```

## Leapfrog Step with Robert-Asselin Filter

```@docs
leapfrog_step!
```

**The Leapfrog scheme:**
```math
\phi^{n+1} = \phi^{n-1} + 2\Delta t \times F^n
```

**Robert-Asselin filter (damps computational mode):**
```math
\tilde{\phi}^n = \phi^n + \gamma(\phi^{n-1} - 2\phi^n + \phi^{n+1})
```

**With integrating factor for hyperdiffusion:**
```math
\phi^{n+1} = \phi^{n-1} \times e^{-2\lambda\Delta t} + 2\Delta t \times F^n \times e^{-\lambda\Delta t}
```

**Usage:**
```julia
# Serial mode
leapfrog_step!(Snp1, Sn, Snm1, grid, params, plans; a=a_ell, dealias_mask=L)

# Parallel mode (2D decomposition)
leapfrog_step!(Snp1, Sn, Snm1, grid, params, plans; a=a_ell, dealias_mask=L, workspace=workspace)

# Time level rotation after each step
Snm1, Sn, Snp1 = Sn, Snp1, Snm1
```

## Forward Euler Update

Used for the first (projection) step:
```math
\phi^{n+1} = \left[\phi^n - \Delta t \times F\right] \times e^{-\lambda\Delta t}
```

The integrating factor `e^{-λΔt}` handles hyperdiffusion exactly.

## Tendency Computation

### compute_rhs!

```@docs
compute_rhs_qg!
compute_rhs_wave!
```

Computes the right-hand side tendencies:

**QG:**
```math
F_q = -J(\psi, q) - J(\psi, q^w) + \text{dissipation}
```

**Wave:**
```math
F_B = -J(\psi, B) - B\frac{\partial\zeta}{\partial t} - i\frac{N^2}{2f_0}\nabla^2 A + \text{dissipation}
```

## Integrating Factors

### Purpose

For stiff diffusion terms, we transform:
```math
\tilde{q} = q \cdot e^{\nu k^{2p} t}
```

### Functions

```@docs
compute_integrating_factors
apply_integrating_factor!
remove_integrating_factor!
```

**Usage:**
```julia
# Setup (once)
IF_q, IF_B = compute_integrating_factors(grid, params, dt)

# Each step
apply_integrating_factor!(q, IF_q)
# ... time step ...
remove_integrating_factor!(q, IF_q)
```

## Startup Procedure

### First Steps

```@docs
startup_ab3!
```

AB3 requires history. For the first two steps:

1. **Step 1**: Forward Euler
2. **Step 2**: AB2

```julia
# Automatic handling
for step = 1:nsteps
    timestep!(state, grid, params, work, plans, a_ell, dt)
    # First 2 steps use Euler/AB2 automatically
end
```

## CFL Condition

### Stability Constraint

```julia
function compute_cfl(state, grid, dt)
    u_max = maximum(abs.(state.u))
    v_max = maximum(abs.(state.v))
    return dt * max(u_max/grid.dx, v_max/grid.dy)
end
```

For stability, CFL < 1 is required. Recommended: CFL ≈ 0.5.

### Adaptive Time Stepping

```julia
function adaptive_dt(state, grid; cfl_target=0.5, dt_max=0.01)
    u_max = maximum(abs.(state.u)) + 1e-10  # Avoid division by zero
    v_max = maximum(abs.(state.v)) + 1e-10

    dt = cfl_target * min(grid.dx/u_max, grid.dy/v_max)
    return min(dt, dt_max)
end
```

## Sub-stepping

For very stiff problems:

```@docs
substep!
```

```julia
# Take N substeps within one outer step
n_sub = 4
dt_sub = dt / n_sub

for _ in 1:n_sub
    substep!(state, grid, params, work, plans, a_ell, dt_sub)
end
```

## State History

### Managing History Arrays

The state contains history for AB3:

```julia
# Current tendency
state.rq_new

# Previous tendencies
state.rq_old   # n-1
state.rq_old2  # n-2
```

### Rotating History

```@docs
rotate_history!
```

Called automatically at end of timestep:
```julia
state.rq_old2 .= state.rq_old
state.rq_old .= state.rq_new
```

## Performance

### Timing Breakdown

Typical distribution:
| Component | Fraction |
|:----------|:---------|
| FFTs | 40-50% |
| Elliptic solves | 20-30% |
| Array operations | 15-25% |
| History management | 5% |

### Optimization

```julia
# Pre-compute integrating factors
IF_q, IF_B = compute_integrating_factors(grid, params, dt)

# Reuse for all steps (if dt is constant)
for step = 1:nsteps
    timestep_with_IF!(state, ..., IF_q, IF_B)
end
```

## API Reference

```@docs
timestep!
ab3_step!
ab2_step!
euler_step!
compute_rhs_qg!
compute_rhs_wave!
compute_integrating_factors
apply_integrating_factor!
startup_ab3!
```
