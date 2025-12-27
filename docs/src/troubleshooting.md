## Troubleshooting

### Package Instantiation Errors

- Symptom: `Package <X> is a direct dependency, but does not appear in the manifest`
- Fix: run `Pkg.resolve(); Pkg.instantiate(); Pkg.precompile()` in the project.

### Missing NCDatasets (NetCDF I/O)

- Symptom: `NCDatasets not available. Install NCDatasets.jl or skip NetCDF I/O.`
- Fix: `julia --project=. -e 'using Pkg; Pkg.add("NCDatasets")'`
  - Alternatively, disable NetCDF I/O by setting `save_*` to `false` in `OutputConfig`.

### FFT Issues on HPC/Clusters

- Install system FFTW or let FFTW.jl download binaries automatically.
- Ensure `using FFTW` works and that shared libraries are available in your environment.

### MPI/Pencil Setup

- Install MPI.jl and PencilArrays/PencilFFTs.
- Launch with `mpiexec -n <N> julia --project=. your_script.jl` and pass `use_mpi=true` to `setup_simulation`.

### Stability/Time Step

- If you see blow‑ups or NaNs:
  - Reduce `dt`
  - Increase resolution
  - Enable viscosity/hyperdiffusion (see parameters in `QGParams`)
  - Start with `linear=true` to isolate linear dynamics

### 2D MPI Decomposition Issues

#### Segmentation Faults or Memory Corruption

- **Symptom**: Crashes with "malloc(): invalid next size" or segfaults when using many MPI ranks
- **Cause**: Spectral and physical arrays have different local dimensions in 2D decomposition

Always get dimensions from the array you're iterating over:

```julia
# Get physical array dimensions after backward FFT
phys_arr = allocate_fft_backward_dst(spectral_arr, plans)
fft_backward!(phys_arr, spectral_arr, plans)
nz_phys, nx_phys, ny_phys = size(parent(phys_arr))
for k in 1:nz_phys, j in 1:ny_phys, i in 1:nx_phys
    phys_arr[k, i, j] = ...
end
```

#### Pencil Topology Mismatch Error

- **Symptom**: `ArgumentError: pencil topologies must be the same`
- **Cause**: Using `deepcopy(state)` instead of `copy_state(state)` for leapfrog time stepping

Use `copy_state(S)` instead of `deepcopy(S)` to preserve pencil topology.

### Wave Dispersion CFL Constraint

For the YBJ+ model, the wave dispersion term imposes a CFL-like constraint:

```
dt ≤ 2f₀/N²
```

For typical ocean parameters (f₀ = 10⁻⁴ s⁻¹, N² = 10⁻⁵ s⁻²), this gives dt ≤ 20s.
If your simulation blows up rapidly, check this constraint.

### Hyperdiffusion for Stability

Add 4th order (biharmonic) hyperdiffusion to damp grid-scale noise:

```julia
par = default_params(
    Lx=70e3, Ly=70e3, Lz=3000.0,
    nx=128, ny=128, nz=64,
    νₕ₁ʷ = 1.0e7,  # Wave hyperdiffusion [m⁴/s]
    ilap1w = 2,    # 4th order (biharmonic)
    γ = 0.01       # Robert-Asselin filter (default: 1e-3)
)
```

