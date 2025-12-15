#=
================================================================================
    Asselin et al. (2020) JPO Dipole Example
================================================================================

This example reproduces the barotropic dipole simulation from:

    Asselin, O., L. N. Thomas, W. R. Young, and L. Rainville (2020)
    "Refraction and Straining of Near-Inertial Waves by Barotropic Eddies"
    Journal of Physical Oceanography, 50, 3439-3454
    DOI: 10.1175/JPO-D-20-0109.1

PHYSICAL SETUP:
---------------
- Location: NISKINe region (~58.5°N, 500 km south of Iceland)
- Domain: 70 km × 70 km horizontally, 3 km depth
- Flow: Steady barotropic vortex dipole
- Waves: Surface-confined near-inertial oscillations (k=0 initially)
- Stratification: Uniform N² = 10⁻⁵ s⁻²

NONDIMENSIONALIZATION:
----------------------
We use Ro=1, Bu=1 by default so that:
- The dispersion coefficient 1/(2*Bu*Ro) = 0.5 = N²/(2f) when N²=f=1
- Physical parameters enter through the flow amplitude and wave parameters
- Domain is [0, 2π]³ in nondimensional coordinates

================================================================================
=#

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
const n_inertial_periods = 15  # Run for 15 inertial periods
const dt = 0.001               # Time step (nondimensional)

# Inertial period in nondimensional time: T_inertial = 2π/f = 2π (when f=1)
const T_inertial = 2π
const nt = round(Int, n_inertial_periods * T_inertial / dt)

# Flow parameters (nondimensional)
# The dipole has max velocity U = 1 (nondimensional)
const U_flow = 1.0

# Wave parameters (nondimensional)
# From paper: u_wave/U_flow ≈ 0.3 (10 cm/s wave vs 33.5 cm/s jet)
const u0_wave = 0.3

# Surface layer depth (as fraction of domain)
# Paper: σ = 30 m out of 3000 m depth = 1% of depth
# In nondim: σ_nd = 0.01 × 2π ≈ 0.063
const sigma_z = 0.01 * 2π

println("="^70)
println("Asselin et al. (2020) Dipole Simulation")
println("="^70)

println("\nSimulation Setup:")
@printf("  Resolution:    %d × %d × %d\n", nx, ny, nz)
@printf("  Time step dt:  %.4f\n", dt)
@printf("  Total steps:   %d\n", nt)
@printf("  Duration:      %.1f inertial periods\n", n_inertial_periods)
@printf("  Wave amplitude: %.2f (relative to flow)\n", u0_wave)
@printf("  Surface layer:  %.3f (nondim)\n", sigma_z)

# Create parameter struct with Ro=Bu=1 (dimensional-like)
par = QGYBJ.QGParams{Float64}(
    # Domain
    nx = nx, ny = ny, nz = nz,
    Lx = 2π, Ly = 2π,

    # Time stepping
    dt = dt,
    nt = nt,

    # Physical parameters - use Ro=Bu=1 for simplicity
    f0 = 1.0,
    Ro = 1.0,  # Rossby number = 1 (dimensional-like)
    Bu = 1.0,  # Burger number = 1 (dimensional-like)
    W2F = u0_wave^2,  # Wave-to-flow ratio squared
    gamma = 1e-3,  # Robert-Asselin filter

    # Hyperdiffusion (for numerical stability)
    nuh1 = 1e-4, ilap1 = 2,    # Biharmonic for flow
    nuh2 = 1e-2, ilap2 = 6,    # Hyper-6 for flow
    nuh1w = 0.0, ilap1w = 2,   # Wave hyperdiffusion
    nuh2w = 1e-2, ilap2w = 6,
    nuz = 0.0,                  # No vertical diffusion

    # Legacy viscosity (unused)
    nu_h = 0.0, nu_v = 0.0,

    # Physics switches
    linear_vert_structure = 0,
    stratification = :constant_N,  # Uniform N² = 1 (nondimensional)
    inviscid = false,
    linear = false,                # Include nonlinear advection
    no_dispersion = false,         # Include wave dispersion
    passive_scalar = false,        # Waves are dynamically active
    ybj_plus = true,               # Use YBJ+ formulation
    no_feedback = true,            # No wave feedback on mean flow
    fixed_flow = true,             # STEADY flow (key assumption!)
    no_wave_feedback = true,

    # Stratification parameters (not used for constant_N)
    N02_sg = 1.0, N12_sg = 0.0, sigma_sg = 1.0, z0_sg = π, alpha_sg = 0.0
)

# ============================================================================
#                       INITIALIZE GRID AND STATE
# ============================================================================

println("\nInitializing grid and state...")

G = QGYBJ.init_grid(par)
S = QGYBJ.init_state(G)

z = G.z
dz = G.nz > 1 ? z[2] - z[1] : 1.0

println("  Grid spacing: dx = $(G.dx), dz = $dz")
println("  Vertical levels: z ∈ [$(minimum(z)), $(maximum(z))]")

# ============================================================================
#                       SET UP DIPOLE STREAMFUNCTION
# ============================================================================
#
# From paper Eq. (2): ψ = U/κ × sin(κx) × cos(κy)
#
# For a 2π periodic domain, use k = 1 (fundamental mode)
# This gives a dipole with:
#   - Anticyclone centered at (x,y) = (π/2, 0)
#   - Cyclone centered at (x,y) = (3π/2, 0)
#   - Jet between them along y = 0

println("\nSetting up dipole streamfunction...")

const k_dipole = 1.0
const psi_amp = U_flow / k_dipole  # ψ = U/κ sin(κx) cos(κy)

# Initialize ψ in real space, then FFT
psi_real = zeros(Float64, nx, ny, nz)

for k in 1:nz, j in 1:ny, i in 1:nx
    x = (i-1) * G.dx
    y = (j-1) * G.dy

    # Dipole streamfunction: ψ = (U/κ) sin(κx) cos(κy)
    # Shift x by π/2 to center anticyclone at origin
    psi_real[i,j,k] = psi_amp * sin(k_dipole * (x - π/2)) * cos(k_dipole * y)
end

# FFT to spectral space
plans = QGYBJ.plan_transforms!(G)
QGYBJ.fft_forward!(S.psi, complex.(psi_real), plans)

# Compute q = ∇²ψ (in spectral space, q = -kh² × ψ)
for k in 1:nz, j in 1:ny, i in 1:nx
    kh2 = G.kx[i]^2 + G.ky[j]^2
    S.q[i,j,k] = -kh2 * S.psi[i,j,k]
end

# Compute max vorticity for verification
zeta_spec = similar(S.q)
for k in 1:nz, j in 1:ny, i in 1:nx
    kh2 = G.kx[i]^2 + G.ky[j]^2
    zeta_spec[i,j,k] = -kh2 * S.psi[i,j,k]
end
zeta_real = similar(S.psi)
QGYBJ.fft_backward!(zeta_real, zeta_spec, plans)
zeta_max = maximum(abs.(real.(zeta_real))) / (nx * ny)

@printf("  Dipole wavenumber κ = %.3f\n", k_dipole)
@printf("  Max vorticity |ζ|/f = %.3f\n", zeta_max)

# ============================================================================
#                       SET UP WAVE INITIAL CONDITION
# ============================================================================
#
# From paper Eq. (4): u(t=0) = u₀ exp(-z²/σ²), v(t=0) = 0
#
# For horizontally uniform initial condition, only (kx,ky)=(0,0) mode is nonzero.
# In QGYBJ.jl, z goes from 0 (bottom) to 2π (top/surface)

println("\nSetting up wave initial condition...")

z_surface = 2π  # Surface location in nondim coordinates

for k in 1:nz
    z_k = z[k]

    # Distance from surface
    depth = z_surface - z_k

    # Gaussian profile peaking at surface: exp(-depth²/σ²)
    wave_profile = exp(-(depth^2) / (sigma_z^2))

    # Set only the (0,0) mode (i=1, j=1)
    S.B[1, 1, k] = u0_wave * wave_profile * (nx * ny)  # FFT normalization factor
end

# Verify wave energy distribution
wave_profile_check = [real(S.B[1,1,k]) / (nx*ny) for k in 1:nz]
@printf("  Wave amplitude at surface: %.3f\n", maximum(wave_profile_check))
@printf("  Wave amplitude at mid-depth: %.4f\n", wave_profile_check[nz÷2])

# ============================================================================
#                       DIAGNOSTIC SETUP
# ============================================================================

# Compute elliptic coefficient a = 1/N² = 1 for constant_N
a_ell = QGYBJ.a_ell_ut(par, G)

# Dealiasing mask
L_mask = QGYBJ.dealias_mask(G)

# Compute velocities from ψ
QGYBJ.compute_velocities!(S, G; plans=plans, params=par)

# Verify flow setup
u_max = maximum(abs.(S.u))
v_max = maximum(abs.(S.v))
@printf("\nFlow verification:\n")
@printf("  Max |u| = %.3f\n", u_max)
@printf("  Max |v| = %.3f\n", v_max)

# ============================================================================
#                       TIME INTEGRATION
# ============================================================================

println("\n" * "="^70)
println("Starting time integration...")
println("="^70)

# Output interval (every inertial period)
output_interval_steps = round(Int, T_inertial / dt)

# Storage for diagnostics
times = Float64[]
wave_energies = Float64[]

# Initial diagnostics
EB = sum(abs2.(S.B)) / (nx * ny * nz)
push!(times, 0.0)
push!(wave_energies, EB)

@printf("\nt = 0.0 IP: E_B = %.4e\n", EB)

# Run projection step (Forward Euler for first step)
println("\nRunning projection step...")
QGYBJ.first_projection_step!(S, G, par, plans; a=a_ell, dealias_mask=L_mask)

# Create states for leapfrog
Sn = deepcopy(S)     # State at time n
Snm1 = deepcopy(S)   # State at time n-1
Snp1 = deepcopy(S)   # State at time n+1

# Main time loop
for step in 1:nt
    # Leapfrog step
    QGYBJ.leapfrog_step!(Snp1, Sn, Snm1, G, par, plans; a=a_ell, dealias_mask=L_mask)

    # Rotate states
    Snm1, Sn, Snp1 = Sn, Snp1, Snm1

    # Output diagnostics
    if step % output_interval_steps == 0
        t_IP = step * dt / T_inertial
        EB = sum(abs2.(Sn.B)) / (nx * ny * nz)

        push!(times, t_IP)
        push!(wave_energies, EB)

        @printf("t = %.1f IP: E_B = %.4e\n", t_IP, EB)
    end
end

println("\n" * "="^70)
println("Simulation complete!")
println("="^70)

# ============================================================================
#                       FINAL ANALYSIS
# ============================================================================

S_final = Sn

# Wave energy at different depths
println("\nWave Energy Distribution (Final State):")
println("-"^50)

# Surface layer (top 10%)
k_surface = round(Int, 0.9 * nz):nz
EB_surface = sum(abs2.(S_final.B[:,:,k_surface]))

# Interior (10-50% depth)
k_interior = round(Int, 0.5 * nz):round(Int, 0.9 * nz)
EB_interior = sum(abs2.(S_final.B[:,:,k_interior]))

# Deep (below 50%)
k_deep = 1:round(Int, 0.5 * nz)
EB_deep = sum(abs2.(S_final.B[:,:,k_deep]))

EB_total = EB_surface + EB_interior + EB_deep

@printf("  Surface layer (top 10%%):    %.1f%%\n", 100 * EB_surface / EB_total)
@printf("  Interior (10-50%% depth):    %.1f%%\n", 100 * EB_interior / EB_total)
@printf("  Deep (below 50%%):           %.1f%%\n", 100 * EB_deep / EB_total)

println("\nKey Physics (from Asselin et al. 2020):")
println("-"^50)
println("1. Wave energy concentrates in anticyclone (ζ < 0)")
println("2. Shear bands form with linear growth of k_h and m")
println("3. For steady barotropic flow, straining is INEFFECTIVE")
println("   because k-advection exactly cancels it")

# Return key objects for further analysis
(S_final, G, par, times, wave_energies)
