# [Stratification](@id stratification)

```@meta
CurrentModule = QGYBJ
```

This page explains how to configure ocean stratification profiles in QGYBJ.jl.

## Why Stratification Matters

The buoyancy frequency ``N(z)`` affects:

- **Wave propagation**: Dispersion depends on ``N^2``
- **Vertical structure**: Mode shapes vary with ``N(z)``
- **Refraction**: Waves bend toward regions of lower ``N``
- **Energy flux**: Vertical group velocity scales with ``N``

## Built-in Profiles

### Constant N

Uniform stratification throughout the water column:

```julia
setup_stratification!(grid, params, :constant_N)
```

Profile:
```math
N^2(z) = N_0^2 = \text{const}
```

Best for:
- Idealized studies
- Analytical comparisons
- Simple mode structure

### Exponential

Stratification decreasing exponentially with depth:

```julia
setup_stratification!(grid, params, :exponential;
    depth_scale = 0.1  # e-folding scale
)
```

Profile:
```math
N^2(z) = N_0^2 \, e^{z/d}
```

where ``d`` is the depth scale and ``z`` is negative (depth).

Best for:
- Realistic ocean upper thermocline
- Simple continuous variation

### Skewed Gaussian (Pycnocline)

Sharp pycnocline with gradual decrease below:

```julia
setup_stratification!(grid, params, :skewed_gaussian;
    pycnocline_depth = 0.05,   # Depth of N² maximum
    pycnocline_width = 0.02,   # Width of pycnocline
    deep_N2_ratio = 0.01       # N²(deep) / N²(surface)
)
```

Profile (schematic):
```
N² →
│
│    ╱╲
│   ╱  ╲
│  ╱    ╲____
│ ╱           ╲____
│╱                  ──────
└─────────────────────────── z (depth)
     ↑
  pycnocline
```

Best for:
- Realistic subtropical ocean
- Strong near-surface trapping
- Wave focusing studies

### Two-Layer

Step-like stratification:

```julia
setup_stratification!(grid, params, :two_layer;
    interface_depth = 0.2,     # Depth of interface
    upper_N2 = 1.0,            # N² in upper layer
    lower_N2 = 0.01            # N² in lower layer
)
```

Profile:
```
N² →
│
│ ────────────
│             │
│             │
│             ─────────────
│
└────────────────────────── z
```

Best for:
- Simple layered models
- Mode-1 dominated dynamics
- Pedagogical examples

## Custom Profiles

### From Function

```julia
# Define custom N² profile
function my_N2(z; N0=1.0, z_pyc=-0.1, sigma=0.02)
    return N0^2 * exp(-((z - z_pyc)/sigma)^2)
end

# Evaluate on grid
N2_profile = [my_N2(z) for z in grid.z]

# Apply to grid
set_stratification!(grid, N2_profile)
```

### From Data

```julia
using DelimitedFiles

# Load from file (depth, N2 columns)
data = readdlm("N2_profile.txt")
z_data = data[:, 1]
N2_data = data[:, 2]

# Interpolate to grid
using Interpolations
itp = linear_interpolation(z_data, N2_data)
N2_profile = [itp(z) for z in grid.z]

set_stratification!(grid, N2_profile)
```

### From Observations

```julia
using NCDatasets

# Load from NetCDF (e.g., Argo profile)
ds = NCDataset("argo_profile.nc")
depth = ds["depth"][:]
temp = ds["temperature"][:]
salt = ds["salinity"][:]
close(ds)

# Compute N² from T, S profiles
N2_profile = compute_N2_from_TS(depth, temp, salt)

# Interpolate and apply
# ... (similar to above)
```

## Stratification Structure

### Internal Storage

The grid stores stratification at cell centers and faces:

```julia
grid.N2      # N² at cell centers [nz]
grid.N2_face # N² at cell faces [nz+1]
grid.a       # f₀²/N² at cell centers [nz]
grid.a_face  # f₀²/N² at cell faces [nz+1]
```

### Accessing Values

```julia
# Get N² profile
N2 = get_stratification(grid)

# Get at specific depth
z_target = -0.5
N2_at_z = interpolate_N2(grid, z_target)
```

## Effects on Dynamics

### Vertical Modes

Stratification determines the vertical mode structure:

```julia
# Compute vertical modes
modes, eigenvalues = compute_vertical_modes(grid, params)

# modes[:, n] is the n-th mode shape
# eigenvalues[n] is the corresponding Rossby radius⁻²
```

Mode-1 dominates when:
- Strong pycnocline exists
- ``N^2`` decreases smoothly with depth

### Wave Trapping

Strong surface stratification traps waves near the surface:

```
Strong N²    Weak N²
near surface everywhere

    │           │
  ──┴──       ──┴──
  Wave        Wave
  trapped     penetrates
  above       to depth
  pycnocline
```

### Deformation Radius

The first baroclinic deformation radius:

```math
R_d = \frac{1}{f_0\pi}\int_0^H N(z) \, dz
```

Computed as:
```julia
Rd = compute_deformation_radius(grid, params)
```

## Visualization

### Profile Plot

```julia
using Plots

z = grid.z
N2 = grid.N2

plot(sqrt.(N2), z,
    xlabel = "N (s⁻¹)",
    ylabel = "Depth",
    title = "Stratification Profile",
    legend = false,
    flip = true  # Depth increases downward
)
```

### With Modes

```julia
modes, _ = compute_vertical_modes(grid, params)

p = plot(layout=(1, 2))

# Left: N² profile
plot!(p[1], N2, z, xlabel="N²", ylabel="z")

# Right: First 3 modes
for n in 1:3
    plot!(p[2], modes[:, n], z, label="Mode $n")
end
xlabel!(p[2], "Mode amplitude")
```

## Best Practices

### Resolution Guidelines

| Profile Type | Recommended ``n_z`` |
|:-------------|:--------------------|
| Constant N | 16-32 |
| Smooth exponential | 32-64 |
| Sharp pycnocline | 64-128 |
| Two-layer | 32-64 |

Sharp features need higher resolution to avoid Gibbs phenomena.

### Smoothness

Ensure ``N^2(z)`` is smooth for numerical stability:

```julia
# Bad: Discontinuous
N2[z .> -0.1] .= 1.0
N2[z .<= -0.1] .= 0.01

# Good: Smooth transition
N2 = 0.505 .+ 0.495 .* tanh.((z .+ 0.1) ./ 0.02)
```

### Physical Constraints

- ``N^2 > 0`` everywhere (stable stratification)
- ``N^2`` should decrease with depth (typically)
- Avoid very small ``N^2`` (causes numerical issues in inversions)

```julia
# Ensure minimum N²
N2_min = 1e-6
N2 = max.(N2, N2_min)
```

## Stratification Types Reference

| Type | Parameters | Best For |
|:-----|:-----------|:---------|
| `:constant_N` | `N0` | Idealized |
| `:exponential` | `depth_scale` | Thermocline |
| `:skewed_gaussian` | `pycnocline_depth`, `pycnocline_width` | Realistic ocean |
| `:two_layer` | `interface_depth`, `upper_N2`, `lower_N2` | Layered models |
| `:custom` | `N2_profile` | Observations |

## Related Topics

- [QG Equations](@ref qg-equations): How N² enters PV
- [YBJ+ Wave Model](@ref ybj-plus): Wave dispersion with N²
- [Configuration](@ref configuration): Setting up simulations
