"""
Elliptic inversion routines, adapted from ellliptic.f90:
Solve for ψ in
    a(z) ∂²ψ/∂z² + b(z) ∂ψ/∂z - k_h² ψ = q,
for each horizontal wavenumber (kx,ky). We apply a tridiagonal solver along z
for each (kx,ky) independently.
"""
module Elliptic

using ..QGYBJ: Grid, State

"""
    invert_q_to_psi!(S, G; a, b)

Invert spectral PV `q(kx,ky,z)` to obtain `psi(kx,ky,z)`. The vertical operator
coefficients `a(z)` and `b(z)` should be provided as vectors of length `nz` and
`nz` respectively (interpreted with second-order finite differences). For a
basic QG model with constant stratification, set `a .= 1`, `b .= 0`.
"""
function invert_q_to_psi!(S::State, G::Grid; a::AbstractVector, b::AbstractVector)
    nx, ny, nz = G.nx, G.ny, G.nz
    @assert length(a) == nz
    @assert length(b) == nz
    ψ = S.psi
    q = S.q

    # Build constant vertical stencil coefficients (Thomas algorithm) per (i,j)
    dl = zeros(eltype(a), nz)   # lower diag
    d  = zeros(eltype(a), nz)   # main diag
    du = zeros(eltype(a), nz)   # upper diag

    # Interior points: a d2/dz2 + b d/dz - kh2 ψ
    # Using second-order central differences on uniform z for now
    # Handle simple boundary conditions ψ_z = 0 at top/bottom (Neumann)
    dz = G.dz
    # For now assume uniform dz for stencil estimates (extend later)
    Δ = nz > 1 ? (G.z[2]-G.z[1]) : 1.0

    for j in 1:ny, i in 1:nx
        kh2 = G.kh2[i,j]
        # Fill diagonals
        fill!(dl, 0); fill!(d, 0); fill!(du, 0)
        # Top boundary (Neumann): approximate ∂ψ/∂z = 0 -> ψ₀ = ψ₁
        d[1]  = a[1]/Δ^2 + kh2
        du[1] = -a[1]/Δ^2
        # Interior
        @inbounds for k in 2:nz-1
            dzz = a[k]/Δ^2
            dz1 = b[k]/(2Δ)
            dl[k] = -dzz - dz1
            d[k]  =  2dzz + kh2
            du[k] = -dzz + dz1
        end
        # Bottom boundary (Neumann): ψ_{nz} = ψ_{nz-1}
        dl[nz] = -a[nz]/Δ^2
        d[nz]  = a[nz]/Δ^2 + kh2

        # RHS is q(i,j,:)
        rhs = view(q, i, j, :)
        sol = view(ψ, i, j, :)
        thomas_solve!(sol, dl, d, du, rhs)
    end
    return S
end

"""
    thomas_solve!(x, dl, d, du, b)

In-place Thomas algorithm for tridiagonal systems with diagonals `dl,d,du`.
Accepts vector views for `x` and `b`.
"""
function thomas_solve!(x, dl, d, du, b)
    n = length(d)
    c = copy(du)
    bb = copy(d)
    x .= b
    # Forward sweep
    c[1] /= bb[1]
    x[1] /= bb[1]
    @inbounds for i in 2:n
        denom = bb[i] - dl[i]*c[i-1]
        c[i] /= denom
        x[i] = (x[i] - dl[i]*x[i-1]) / denom
    end
    # Back substitution
    @inbounds for i in n-1:-1:1
        x[i] -= c[i]*x[i+1]
    end
    return x
end

end # module

using .Elliptic: invert_q_to_psi!

