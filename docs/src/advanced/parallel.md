# [MPI Parallelization](@id parallel)

```@meta
CurrentModule = QGYBJ
```

This page describes how to run QGYBJ.jl on distributed memory systems using MPI.

## Overview

QGYBJ.jl uses **pencil decomposition** for MPI parallelization:

```
              Full Domain                    Pencil Decomposition
        ┌─────────────────────┐         ┌───────┬───────┬───────┐
        │                     │         │ P0    │ P1    │ P2    │
        │                     │         │───────┼───────┼───────│
        │    nx × ny × nz     │   →     │ P3    │ P4    │ P5    │
        │                     │         │───────┼───────┼───────│
        │                     │         │ P6    │ P7    │ P8    │
        └─────────────────────┘         └───────┴───────┴───────┘
```

Each MPI rank owns a "pencil" of data spanning one full dimension.

## Requirements

Install MPI-related packages:

```julia
using Pkg
Pkg.add(["MPI", "PencilArrays", "PencilFFTs"])
```

Ensure you have a working MPI installation (OpenMPI, MPICH, etc.).

## Quick Start

### Basic MPI Script

```julia
# parallel_run.jl
using QGYBJ
using MPI

# Initialize MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
nprocs = MPI.Comm_size(comm)

# Create configuration
config = create_simple_config(
    nx = 256,
    ny = 256,
    nz = 128,
    dt = 0.001,
    total_time = 10.0,
    parallel = true,
    comm = comm
)

# Run simulation (automatically parallelized)
result = run_simple_simulation(config)

# Output from rank 0 only
if rank == 0
    println("Simulation complete!")
    println("Final KE: ", result.diagnostics.KE[end])
end

MPI.Finalize()
```

### Running

```bash
# Run on 16 processes
mpiexec -n 16 julia --project parallel_run.jl
```

## Pencil Decomposition

### How It Works

The domain is decomposed into 2D slabs ("pencils"):

| Configuration | Layout | Best For |
|:--------------|:-------|:---------|
| Pencil-X | Full in X, split in Y,Z | X-direction FFTs |
| Pencil-Y | Full in Y, split in X,Z | Y-direction FFTs |
| Pencil-Z | Full in Z, split in X,Y | Vertical operations |

### Automatic Transposition

FFTs require data along the transform dimension. QGYBJ.jl automatically:

1. Starts in Pencil-X for X-FFT
2. Transposes to Pencil-Y for Y-FFT
3. Transposes to Pencil-Z for vertical operations
4. Transposes back for next time step

### Decomposition Setup

```julia
# Manual decomposition control
decomp = PencilDecomposition(
    comm,
    (nx, ny, nz),
    proc_dims = (4, 4)  # 4×4 = 16 processes
)

grid = create_mpi_grid(decomp, Lx=2π, Ly=2π, H=1.0)
```

## Full MPI Example

```julia
using QGYBJ
using MPI

function main()
    MPI.Init()

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    # Global parameters
    nx, ny, nz = 256, 256, 128
    dt = 0.001
    nsteps = 10000

    # Set up MPI environment
    mpi_config = setup_mpi_environment(comm; nx=nx, ny=ny, nz=nz)

    # Create grid with decomposition
    grid = init_mpi_grid(mpi_config, Lx=2π, Ly=2π, H=1.0)

    # Parameters
    params = QGParams(; ybj_plus=true, nu_h2=1e-10)

    # Create distributed state
    state = init_mpi_state(grid, mpi_config)

    # Initialize (each rank initializes its local portion)
    initialize_random_flow!(state, grid; seed=rank+42)
    initialize_random_waves!(state, grid; seed=rank+123)

    # FFT plans (MPI-aware)
    plans = plan_mpi_transforms!(grid, mpi_config)

    # Work arrays
    work = create_mpi_work_arrays(grid, mpi_config)
    a_ell = setup_elliptic_matrices(grid, params)

    # Time stepping
    for step = 1:nsteps
        timestep!(state, grid, params, work, plans, a_ell, dt)

        if step % 100 == 0 && rank == 0
            KE = mpi_reduce_sum(flow_kinetic_energy(state.u, state.v), comm)
            println("Step $step: KE = $KE")
        end

        # Checkpointing
        if step % 1000 == 0
            save_mpi_checkpoint(state, grid, step, mpi_config)
        end
    end

    MPI.Finalize()
end

main()
```

## Domain Decomposition Details

### Process Topology

```julia
# Automatic topology (recommended)
decomp = PencilDecomposition(comm, (nx, ny, nz))

# Custom topology (for specific architectures)
decomp = PencilDecomposition(comm, (nx, ny, nz);
    proc_dims = (Px, Py)
)

# Constraint: Px × Py = nprocs
```

### Local Array Size

Each rank holds a portion of the global array:

```julia
# Global size
global_size = (nx, ny, nz)

# Local size (varies by rank)
local_size = size(state.psi)  # e.g., (64, 64, 128) for 16 procs
```

### Index Mapping

```julia
# Get global indices for local data
i_global, j_global, k_global = local_to_global_indices(grid, mpi_config)

# Example: local point (i, j, k) maps to global (i_global[i], j_global[j], k)
```

## Communication Patterns

### Halo Exchange

For stencil operations (e.g., finite differences):

```julia
# Exchange halo regions between neighbors
halo_exchange!(field, grid, mpi_config)
```

### Global Reductions

```julia
# Sum across all ranks
global_sum = mpi_reduce_sum(local_value, comm)

# Maximum across all ranks
global_max = MPI.Allreduce(local_max, MPI.MAX, comm)
```

### Gather/Scatter

```julia
# Gather full field to rank 0
if rank == 0
    global_field = zeros(ComplexF64, nx, ny, nz)
end
gather_to_root!(global_field, local_field, mpi_config)

# Scatter from rank 0 to all
scatter_from_root!(local_field, global_field, mpi_config)
```

## I/O in Parallel

### Parallel NetCDF

```julia
using NCDatasets

# Collective file creation
ds = NCDataset("output.nc", "c";
    format = :netcdf4,
    comm = comm
)

# Each rank writes its portion
local_range = get_local_range(grid, mpi_config)
ds["psi"][local_range..., t] = state.psi

close(ds)
```

### Gather-and-Write

Alternative for smaller datasets:

```julia
if rank == 0
    global_psi = zeros(ComplexF64, nx, ny, nz)
end

gather_to_root!(global_psi, state.psi, mpi_config)

if rank == 0
    # Write from rank 0 only
    NCDataset("output.nc", "a") do ds
        ds["psi"][:, :, :, t] = global_psi
    end
end
```

## Performance Considerations

### Scaling

| Regime | Characteristic |
|:-------|:---------------|
| Strong scaling | Fixed problem size, increasing procs |
| Weak scaling | Fixed work per proc, increasing total size |

QGYBJ.jl scales well up to ~1000 cores for typical configurations.

### Bottlenecks

1. **Transposition**: All-to-all communication for pencil swaps
2. **I/O**: Parallel file writes can be slow
3. **Load imbalance**: Uneven decomposition

### Optimization Tips

```julia
# 1. Use power-of-2 process counts
mpiexec -n 16 julia script.jl  # Good
mpiexec -n 15 julia script.jl  # May be slower

# 2. Match decomposition to architecture
# For 256-core nodes with 4 sockets:
decomp = PencilDecomposition(comm, dims; proc_dims=(16, 16))

# 3. Minimize I/O frequency
output_interval = 1000  # Not too frequent

# 4. Use asynchronous I/O if available
enable_async_io!(output_config)
```

### Profiling

```julia
# Built-in timing
using MPI
using QGYBJ

MPI.Init()

# Enable timing
set_timing!(true)

# Run simulation
result = run_simulation(config)

# Print timing summary
if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    print_timing_summary()
end
```

## Job Submission

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=qgybj
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --time=24:00:00

module load julia
module load openmpi

mpiexec -n 64 julia --project run_simulation.jl
```

### PBS Example

```bash
#!/bin/bash
#PBS -N qgybj
#PBS -l nodes=4:ppn=16
#PBS -l walltime=24:00:00

module load julia openmpi

cd $PBS_O_WORKDIR
mpiexec -n 64 julia --project run_simulation.jl
```

## Troubleshooting

### Common Issues

**MPI not found**
```julia
# Install system MPI package first, then:
using MPIPreferences
MPIPreferences.use_system_binary()
```

**Deadlock in communication**
- Ensure all ranks call collective operations
- Check for mismatched send/receive counts

**Memory errors**
- Local arrays might be too large
- Increase number of processes

### Debugging

```julia
# Debug output with rank info
function mpi_println(msg)
    rank = MPI.Comm_rank(MPI.COMM_WORLD)
    println("[Rank $rank] $msg")
    flush(stdout)
end

# Synchronize before critical sections
MPI.Barrier(comm)
mpi_println("Reached checkpoint")
```

## API Reference

```@docs
setup_mpi_environment
init_mpi_grid
init_mpi_state
plan_mpi_transforms!
gather_to_root!
scatter_from_root!
mpi_reduce_sum
```
