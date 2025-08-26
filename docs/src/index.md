# QGYBJ.jl

```@meta
CurrentModule = QGYBJ
```

QGYBJ.jl is a Julia implementation of a quasi‑geostrophic (QG) model
coupled with the Young–Ben Jelloul (YBJ) wave envelope dynamics (YBJ+).
It mirrors the structure of a reference Fortran code, providing:

- Spectral 2D FFTs per vertical level and vertical tridiagonal solvers
- QG streamfunction inversion (q → ψ) and YBJ+ recovery (B=L⁺A → A)
- Nonlinear advection, refraction, and wave feedback operators
- Diagnostics, vertical velocity (QG omega or YBJ form)
- A high‑level configuration API + examples
- Optional NetCDF I/O for state and diagnostics
- Hooks for parallel execution and particle advection

Use the left‑hand navigation to dive into setup, configuration, I/O,
and a complete API reference.

If you’re new, start with Getting Started.

