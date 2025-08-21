"""
Time stepping for QG–YBJ (projection + leapfrog). Implements:
- Nonlinear advection J(psi, q) and J(psi, B)
- YBJ+ inversion B -> A and C=A_z
- Hyperdiffusion via integrating factors (horizontal)
- Robert–Asselin filter
- 2/3 dealiasing mask on updates
Refraction and feedback terms can be extended to match the full Fortran model.
"""

"""
    first_projection_step!(S, G, par, plans; a)

Projection step (Forward Euler) example: given q^n, invert to ψ^n, compute
diagnostics. Extend to compute B and apply LA->A inversion similarly to Fortran.
"""
function first_projection_step!(S::State, G::Grid, par::QGParams, plans; a, dealias_mask=nothing)
    invert_q_to_psi!(S, G; a)
    invert_B_to_A!(S, G, par, a)
    compute_velocities!(S, G; plans)
    return S
end

"""
    leapfrog_step!(Snp1, Sn, Snm1, G, par, plans; a)

Skeleton for leapfrog advance; fill with actual tendencies from Fortran model.
Here just carries ψ diagnostic refresh.
"""
function leapfrog_step!(Snp1::State, Sn::State, Snm1::State,
                        G::Grid, par::QGParams, plans; a, dealias_mask=nothing)
    nx, ny, nz = G.nx, G.ny, G.nz
    # Dealiasing mask
    Lmask = isnothing(dealias_mask) ? trues(nx,ny) : dealias_mask
    # Nonlinear terms
    nqk = similar(Sn.q)
    nBk = similar(Sn.B)
    jacobian_spectral!(nqk, Sn.psi, Sn.q, G, plans)      # J(psi, q)
    jacobian_spectral!(nBk, Sn.psi, Sn.B, G, plans)      # J(psi, B)
    # Refraction term placeholder (rBk): set zero; user to extend
    rBk = fill!(similar(Sn.B), 0)
    # Hyperdiffusion integrating factors
    int_q   = zeros(Float64, nx, ny)
    int_B   = zeros(Float64, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        int_q[i,j] = par.dt * ( par.nu_h * (abs(G.kx[i])^(2) + abs(G.ky[j])^(2)) )
        int_B[i,j] = par.dt * ( par.nu_h * (abs(G.kx[i])^(2) + abs(G.ky[j])^(2)) )
    end
    # Leapfrog formula with integrating factor and 2/3 dealiasing on horizontal modes
    qtemp = similar(Sn.q)
    Btemp = similar(Sn.B)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        if Lmask[i,j]
            efq = exp(-int_q[i,j])
            efB = exp(-int_B[i,j])
            qtemp[i,j,k] = Snm1.q[i,j,k]*exp(-2*int_q[i,j]) - 2*par.dt*nqk[i,j,k]*efq
            Btemp[i,j,k] = Snm1.B[i,j,k]*exp(-2*int_B[i,j]) - 2*par.dt*( nBk[i,j,k] )*efB
        else
            qtemp[i,j,k] = 0
            Btemp[i,j,k] = 0
        end
    end
    # Robert–Asselin filter
    γ = 1e-3
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        if Lmask[i,j]
            Snm1.q[i,j,k] = Sn.q[i,j,k] + γ*( Snm1.q[i,j,k] - 2Sn.q[i,j,k] + qtemp[i,j,k] )
            Snm1.B[i,j,k] = Sn.B[i,j,k] + γ*( Snm1.B[i,j,k] - 2Sn.B[i,j,k] + Btemp[i,j,k] )
        else
            Snm1.q[i,j,k] = 0
            Snm1.B[i,j,k] = 0
        end
    end
    # Accept the new fields
    Snp1.q .= qtemp
    Snp1.B .= Btemp
    # Recover ψ and A, velocities
    invert_q_to_psi!(Snp1, G; a)
    invert_B_to_A!(Snp1, G, par, a)
    compute_velocities!(Snp1, G; plans)
    return Snp1
end
