#=
================================================================================
            MPI Parallel Example for QGYBJ.jl
================================================================================

This example demonstrates how to run QGYBJ.jl with MPI parallelization using
PencilArrays and PencilFFTs for distributed arrays and FFT transforms.

PREREQUISITES:
--------------
Install the required packages:

    using Pkg
    Pkg.add(["MPI", "PencilArrays", "PencilFFTs"])

For MPI, you also need an MPI implementation installed on your system
(OpenMPI, MPICH, or Intel MPI).

RUNNING:
--------
    mpiexecjl -n 4 julia examples/mpi_parallel_example.jl

Or with OpenMPI:
    mpirun -np 4 julia examples/mpi_parallel_example.jl

================================================================================
=#

# Load MPI packages FIRST (before QGYBJ)
using MPI
using PencilArrays
using PencilFFTs

# Now load QGYBJ - the extension will be automatically activated
using QGYBJ

function main()
    # Initialize MPI
    MPI.Init()

    try
        run_parallel_simulation()
    finally
        MPI.Finalize()
    end
end

function run_parallel_simulation()
    #==========================================================================
    STEP 1: Set up MPI environment
    ==========================================================================#
    mpi_config = QGYBJ.setup_mpi_environment()

    # Only print from rank 0
    is_root = mpi_config.is_root
    if is_root
        println("="^60)
        println("QGYBJ.jl MPI Parallel Simulation")
        println("="^60)
        println("Number of processes: ", mpi_config.nprocs)
        println("Process topology: ", mpi_config.topology)
        println()
    end

    #==========================================================================
    STEP 2: Create model parameters
    ==========================================================================#
    # Use larger domain for parallel efficiency
    params = QGYBJ.default_params(
        nx = 128,
        ny = 128,
        nz = 64,
        dt = 0.001,
        nt = 100,
        nu_h = 1e-4,
        nu_v = 1e-4
    )

    if is_root
        println("Grid size: $(params.nx) x $(params.ny) x $(params.nz)")
        println("Time step: $(params.dt)")
        println("Number of steps: $(params.nt)")
        println()
    end

    #==========================================================================
    STEP 3: Initialize parallel grid and state
    ==========================================================================#
    if is_root
        println("Initializing parallel grid...")
    end

    grid = QGYBJ.init_mpi_grid(params, mpi_config)

    if is_root
        println("Initializing parallel state...")
    end

    state = QGYBJ.init_mpi_state(grid, mpi_config)
    state_old = QGYBJ.init_mpi_state(grid, mpi_config)

    # Synchronize all processes
    QGYBJ.mpi_barrier(mpi_config)

    #==========================================================================
    STEP 4: Set up FFT transforms
    ==========================================================================#
    if is_root
        println("Setting up parallel FFT plans...")
    end

    plans = QGYBJ.plan_mpi_transforms(grid, mpi_config)

    #==========================================================================
    STEP 5: Initialize fields with deterministic random values
    ==========================================================================#
    if is_root
        println("Initializing fields...")
    end

    # Initialize streamfunction with random phases (deterministic across ranks)
    QGYBJ.init_mpi_random_field!(state.psi, grid, 0.1, 0)

    # Initialize wave envelope
    QGYBJ.init_mpi_random_field!(state.B, grid, 0.01, 12345)

    QGYBJ.mpi_barrier(mpi_config)

    #==========================================================================
    STEP 6: Compute initial diagnostics
    ==========================================================================#
    # Get local portion of the state for diagnostics
    local_range = QGYBJ.local_indices(grid)

    # Compute local energy contributions
    local_psi_energy = 0.0
    local_B_energy = 0.0

    psi_parent = parent(state.psi)
    B_parent = parent(state.B)

    for idx in eachindex(psi_parent)
        local_psi_energy += abs2(psi_parent[idx])
        local_B_energy += abs2(B_parent[idx])
    end

    # Sum across all processes
    global_psi_energy = QGYBJ.mpi_reduce_sum(local_psi_energy, mpi_config)
    global_B_energy = QGYBJ.mpi_reduce_sum(local_B_energy, mpi_config)

    if is_root
        println()
        println("Initial diagnostics:")
        println("  Total |ψ|² = ", global_psi_energy)
        println("  Total |B|² = ", global_B_energy)
        println()
    end

    #==========================================================================
    STEP 7: Time integration loop
    ==========================================================================#
    if is_root
        println("Starting time integration...")
        println()
    end

    # Physics setup
    a_ell = QGYBJ.a_ell_ut(params, grid)
    L = QGYBJ.dealias_mask(grid)

    # Time integration
    start_time = time()

    for n in 1:params.nt
        if n == 1
            # First step: Forward Euler
            QGYBJ.first_projection_step!(state, state_old, grid, params, plans, L, a_ell)
        else
            # Subsequent steps: Leapfrog with Robert-Asselin filter
            QGYBJ.leapfrog_step!(state, state_old, grid, params, plans, L, a_ell)
        end

        # Synchronize after each step
        QGYBJ.mpi_barrier(mpi_config)

        # Print progress every 10 steps
        if is_root && n % 10 == 0
            elapsed = time() - start_time
            steps_per_sec = n / elapsed

            # Compute current energy
            local_psi = 0.0
            for idx in eachindex(psi_parent)
                local_psi += abs2(psi_parent[idx])
            end
            global_psi = QGYBJ.mpi_reduce_sum(local_psi, mpi_config)

            println("Step $n / $(params.nt): |ψ|² = $(round(global_psi, sigdigits=4)), ",
                    "$(round(steps_per_sec, digits=2)) steps/sec")
        end
    end

    QGYBJ.mpi_barrier(mpi_config)
    total_time = time() - start_time

    #==========================================================================
    STEP 8: Final diagnostics and output
    ==========================================================================#
    # Compute final energy
    local_psi_final = 0.0
    local_B_final = 0.0
    for idx in eachindex(psi_parent)
        local_psi_final += abs2(psi_parent[idx])
        local_B_final += abs2(B_parent[idx])
    end

    global_psi_final = QGYBJ.mpi_reduce_sum(local_psi_final, mpi_config)
    global_B_final = QGYBJ.mpi_reduce_sum(local_B_final, mpi_config)

    if is_root
        println()
        println("="^60)
        println("Simulation complete!")
        println("="^60)
        println("Total time: $(round(total_time, digits=2)) seconds")
        println("Average: $(round(params.nt / total_time, digits=2)) steps/second")
        println()
        println("Final diagnostics:")
        println("  Total |ψ|² = ", global_psi_final)
        println("  Total |B|² = ", global_B_final)
        println()
        println("Energy change:")
        println("  Δ|ψ|² = ", global_psi_final - global_psi_energy)
        println("  Δ|B|² = ", global_B_final - global_B_energy)
    end

    #==========================================================================
    STEP 9: (Optional) Gather and save output
    ==========================================================================#
    if is_root
        println()
        println("Gathering final state to root...")
    end

    # Gather final state to rank 0 for output
    psi_gathered = QGYBJ.gather_to_root(state.psi, grid, mpi_config)

    if is_root && psi_gathered !== nothing
        println("Gathered array size: ", size(psi_gathered))
        # Here you could save to NetCDF or other format
        # NCDatasets.write("output.nc", psi_gathered)
    end

    QGYBJ.mpi_barrier(mpi_config)

    if is_root
        println()
        println("Done!")
    end
end

# Run the main function
main()
