#=
================================================================================
    Asselin et al. (2020) JPO Dipole Example - MPI Parallel Version
================================================================================

MPI-parallel version of the barotropic dipole simulation from:

    Asselin, O., L. N. Thomas, W. R. Young, and L. Rainville (2020)
    "Refraction and Straining of Near-Inertial Waves by Barotropic Eddies"
    Journal of Physical Oceanography, 50, 3439-3454

USAGE:
------
    mpirun -n 4 julia --project examples/asselin_jpo2020_dipole_mpi.jl
    mpirun -n 16 julia --project examples/asselin_jpo2020_dipole_mpi.jl

For optimal performance, use power-of-2 process counts (4, 16, 64, etc.)

================================================================================
=#

using MPI
using PencilArrays
using PencilFFTs
using QGYBJ
using Printf

# ============================================================================
#                       MODEL PARAMETERS
# ============================================================================

# Grid resolution
const nx = 128
const ny = 128
const nz = 64

# Time stepping
const n_inertial_periods = 15
const dt = 0.001
const T_inertial = 2π  # Inertial period when f=1
const nt = round(Int, n_inertial_periods * T_inertial / dt)

# Flow and wave parameters (nondimensional)
const U_flow = 1.0     # Max flow velocity
const u0_wave = 0.3    # Wave amplitude (relative to flow)
const sigma_z = 0.01 * 2π  # Surface layer depth

# ============================================================================
#                       MAIN FUNCTION
# ============================================================================

function main()
    MPI.Init()
    mpi_config = QGYBJ.setup_mpi_environment()
    is_root = mpi_config.is_root

    if is_root
        println("="^70)
        println("Asselin et al. (2020) Dipole Simulation - MPI Parallel")
        println("="^70)
        println("\nRunning on $(mpi_config.nprocs) processes")
        println("Topology: $(mpi_config.topology)")
        println("\nSimulation Setup:")
        @printf("  Resolution:    %d × %d × %d\n", nx, ny, nz)
        @printf("  Time step dt:  %.4f\n", dt)
        @printf("  Total steps:   %d\n", nt)
        @printf("  Duration:      %.1f inertial periods\n", n_inertial_periods)
    end

    # Parameters with Ro=Bu=1 (dimensional-like)
    par = QGYBJ.QGParams{Float64}(
        nx = nx, ny = ny, nz = nz,
        Lx = 2π, Ly = 2π,
        dt = dt, nt = nt,
        f0 = 1.0,
        Ro = 1.0,  # Dimensional-like
        Bu = 1.0,  # Dimensional-like
        W2F = u0_wave^2,
        gamma = 1e-3,
        nuh1 = 1e-4, ilap1 = 2,
        nuh2 = 1e-2, ilap2 = 6,
        nuh1w = 0.0, ilap1w = 2,
        nuh2w = 1e-2, ilap2w = 6,
        nuz = 0.0,
        nu_h = 0.0, nu_v = 0.0,
        linear_vert_structure = 0,
        stratification = :constant_N,
        inviscid = false,
        linear = false,
        no_dispersion = false,
        passive_scalar = false,
        ybj_plus = true,
        no_feedback = true,
        fixed_flow = true,
        no_wave_feedback = true,
        N02_sg = 1.0, N12_sg = 0.0, sigma_sg = 1.0, z0_sg = π, alpha_sg = 0.0
    )

    # Initialize distributed grid and state
    if is_root; println("\nInitializing distributed grid and state..."); end

    G = QGYBJ.init_mpi_grid(par, mpi_config)
    S = QGYBJ.init_mpi_state(G, mpi_config)
    workspace = QGYBJ.init_mpi_workspace(G, mpi_config)
    plans = QGYBJ.plan_mpi_transforms(G, mpi_config)

    local_range = QGYBJ.get_local_range_xy(G)
    z = G.z
    dx, dy = G.dx, G.dy

    # ========================================================================
    #                       SET UP DIPOLE STREAMFUNCTION
    # ========================================================================

    if is_root; println("Setting up dipole streamfunction..."); end

    const k_dipole = 1.0
    const psi_amp = U_flow / k_dipole

    psi_local = parent(S.psi)

    for k_local in axes(psi_local, 3)
        k_global = local_range[3][k_local]
        for j_local in axes(psi_local, 2)
            j_global = local_range[2][j_local]
            y = (j_global - 1) * dy
            for i_local in axes(psi_local, 1)
                i_global = local_range[1][i_local]
                x = (i_global - 1) * dx
                psi_val = psi_amp * sin(k_dipole * (x - π/2)) * cos(k_dipole * y)
                psi_local[i_local, j_local, k_local] = complex(psi_val)
            end
        end
    end

    QGYBJ.fft_forward!(S.psi, S.psi, plans)

    # Compute q = -kh² × ψ
    q_local = parent(S.q)
    psi_local = parent(S.psi)

    for k_local in axes(q_local, 3)
        k_global = local_range[3][k_local]
        for j_local in axes(q_local, 2)
            j_global = local_range[2][j_local]
            ky_val = G.ky[j_global]
            for i_local in axes(q_local, 1)
                i_global = local_range[1][i_local]
                kx_val = G.kx[i_global]
                kh2 = kx_val^2 + ky_val^2
                q_local[i_local, j_local, k_local] = -kh2 * psi_local[i_local, j_local, k_local]
            end
        end
    end

    # ========================================================================
    #                       SET UP WAVE INITIAL CONDITION
    # ========================================================================

    if is_root; println("Setting up wave initial condition..."); end

    z_surface = 2π
    B_local = parent(S.B)

    for k_local in axes(B_local, 3)
        k_global = local_range[3][k_local]
        z_k = z[k_global]
        depth = z_surface - z_k
        wave_profile = exp(-(depth^2) / (sigma_z^2))

        for j_local in axes(B_local, 2)
            j_global = local_range[2][j_local]
            for i_local in axes(B_local, 1)
                i_global = local_range[1][i_local]
                # Only (0,0) mode
                if i_global == 1 && j_global == 1
                    B_local[i_local, j_local, k_local] = u0_wave * wave_profile * (nx * ny)
                else
                    B_local[i_local, j_local, k_local] = 0.0
                end
            end
        end
    end

    # ========================================================================
    #                       TIME INTEGRATION
    # ========================================================================

    a_ell = QGYBJ.a_ell_ut(par, G)
    L_mask = QGYBJ.dealias_mask(G)

    QGYBJ.compute_velocities!(S, G; plans=plans, params=par, workspace=workspace)

    if is_root
        println("\n" * "="^70)
        println("Starting time integration...")
        println("="^70)
    end

    output_interval_steps = round(Int, T_inertial / dt)

    # Initial diagnostics
    local_EB = sum(abs2.(parent(S.B)))
    global_EB = QGYBJ.mpi_reduce_sum(local_EB, mpi_config)

    if is_root
        @printf("\nt = 0.0 IP: E_B = %.4e\n", global_EB / (nx * ny * nz))
    end

    # First step
    if is_root; println("\nRunning projection step..."); end
    QGYBJ.first_projection_step!(S, G, par, plans; a=a_ell, dealias_mask=L_mask, workspace=workspace)

    Sn = deepcopy(S)
    Snm1 = deepcopy(S)
    Snp1 = deepcopy(S)

    # Main loop
    for step in 1:nt
        QGYBJ.leapfrog_step!(Snp1, Sn, Snm1, G, par, plans;
                             a=a_ell, dealias_mask=L_mask, workspace=workspace)
        Snm1, Sn, Snp1 = Sn, Snp1, Snm1

        if step % output_interval_steps == 0
            t_IP = step * dt / T_inertial
            local_EB = sum(abs2.(parent(Sn.B)))
            global_EB = QGYBJ.mpi_reduce_sum(local_EB, mpi_config)

            if is_root
                @printf("t = %.1f IP: E_B = %.4e\n", t_IP, global_EB / (nx * ny * nz))
            end
        end
    end

    if is_root
        println("\n" * "="^70)
        println("Simulation complete!")
        println("="^70)
    end

    MPI.Finalize()
    return nothing
end

main()
