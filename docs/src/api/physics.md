# [Physics Functions](@id api-physics)

```@meta
CurrentModule = QGYBJ
```

This page documents the physics functions in QGYBJ.jl.

## Elliptic Inversions

### Streamfunction Inversion

```@docs
invert_q_to_psi!
```

**Solves:** ``\nabla^2\psi + \frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\frac{\partial\psi}{\partial z}\right) = q``

**Usage:**
```julia
# Serial mode
invert_q_to_psi!(state, grid; a=a_ell)

# Parallel mode (with workspace for 2D decomposition)
invert_q_to_psi!(state, grid; a=a_ell, workspace=workspace)
# Updates state.psi from state.q
```

### Wave Amplitude Inversion

```@docs
invert_B_to_A!
```

**Solves:** ``\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\frac{\partial A}{\partial z}\right) - \frac{k_h^2}{4}A = B``

**Usage:**
```julia
# Serial mode
invert_B_to_A!(state, grid, params, a_ell)

# Parallel mode (with workspace for 2D decomposition)
invert_B_to_A!(state, grid, params, a_ell; workspace=workspace)
# Updates state.A from state.B
```

### Helmholtz Solver

```@docs
invert_helmholtz!
```

## Nonlinear Terms

### Jacobian

```@docs
jacobian_spectral!
```

**Computes:** ``J(a, b) = \frac{\partial a}{\partial x}\frac{\partial b}{\partial y} - \frac{\partial a}{\partial y}\frac{\partial b}{\partial x}``

### Wave Advection

```@docs
convol_waqg!
```

**Wave Advection:** ``J(\psi, B)`` - Advection of wave envelope by streamfunction

### Refraction

```@docs
refraction_waqg!
```

**Computes:** ``B \zeta`` - Wave refraction by vorticity

### Wave Feedback

```@docs
compute_qw!
```

**Computes:** ``q^w`` - Wave feedback on mean flow (Xie & Vanneste 2015)

## Velocity Computation

```@docs
compute_velocities!
```

**Computes:**
- ``u = -\partial\psi/\partial y``
- ``v = \partial\psi/\partial x``

### Vertical Velocity

```@docs
compute_vertical_velocity!
```

### Total Velocities

```@docs
compute_total_velocities!
```

## Dissipation

### Vertical Diffusion

```@docs
dissipation_q_nv!
```

**Applies:** Vertical diffusion term ``\nu_z \partial^2 q / \partial z^2``

### Integrating Factor

```@docs
int_factor
```

For handling stiff diffusion terms efficiently.

## Diagnostics Functions

### Energy

```@docs
wave_energy
flow_kinetic_energy
```

### Omega Equation

```@docs
omega_eqn_rhs!
```

## Transform Functions

### Forward Transforms

```@docs
fft_forward!
```

Real space → Spectral space

### Backward Transforms

```@docs
fft_backward!
```

Spectral space → Real space

### Dealiasing

```@docs
dealias_mask
```

Creates dealiasing mask using 2/3 rule.

## YBJ Normal Mode Functions

```@docs
sumB!
compute_sigma
compute_A!
```

## Function Signatures Summary

| Function | Input | Output | In-place |
|:---------|:------|:-------|:---------|
| `invert_q_to_psi!` | q | psi | Yes |
| `invert_B_to_A!` | B | A, C | Yes |
| `jacobian_spectral!` | a, b | J(a,b) | Yes |
| `compute_velocities!` | psi | u, v | Yes |
| `flow_kinetic_energy` | u, v | scalar | No |
| `wave_energy` | B, A | scalar | No |

## Performance Notes

- All physics functions are **in-place** to avoid allocations
- FFT plans are **pre-computed** for efficiency
- Tridiagonal systems use **Thomas algorithm** (O(n))
- Functions are **type-stable** for optimal JIT compilation

## 2D Decomposition Notes

Functions requiring vertical operations automatically detect 2D decomposition and use the appropriate method:

| Function | Serial | Parallel (2D) |
|:---------|:-------|:--------------|
| `invert_q_to_psi!` | Direct solve | Transpose → solve → transpose |
| `invert_B_to_A!` | Direct solve | Transpose → solve → transpose |
| `compute_vertical_velocity!` | Direct solve | Transpose → solve → transpose |
| `dissipation_q_nv!` | Direct | Transpose if needed |

**Pattern:**
```julia
need_transpose = G.decomp !== nothing && hasfield(typeof(G.decomp), :pencil_z)
if need_transpose
    _function_2d!(...)   # Uses transpose operations
else
    _function_direct!(...)  # Direct vertical access
end
```

**Workspace requirement:** Pass `workspace` argument for parallel mode to avoid repeated allocation of z-pencil arrays.
