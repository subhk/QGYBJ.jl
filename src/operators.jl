"""
Operators and utilities mapping to key Fortran routines in derivatives.f90:
 - compute_streamfunction (q -> psi inversion via elliptic solver; in elliptic.jl)
 - compute_velo (psi -> u,v,b,w diagnostics)
This file wires spectral-to-real conversions and simple diagnostic operators.
"""
module Operators

using ..QGYBJ: Grid, State
using ..QGYBJ: fft_forward!, fft_backward!, plan_transforms!, compute_wavenumbers!

"""
    compute_velocities!(S, G)

Given spectral streamfunction `psi(kx,ky,z)`, compute velocities:
- Horizontal: `u = -∂ψ/∂y`, `v = ∂ψ/∂x` using spectral differentiation
- Vertical: `w` from QG ageostrophic vertical velocity (omega equation)

The vertical velocity computation requires solving the omega equation:
∇²w + (N²/f²)(∂²w/∂z²) = 2 J(ψ_z, ∇²ψ)
"""
function compute_velocities!(S::State, G::Grid; plans=nothing, params=nothing, compute_w=true)
    # Spectral differentiation: û = -i ky ψ̂, v̂ = i kx ψ̂
    ψk = S.psi
    uk = similar(ψk)
    vk = similar(ψk)
    @inbounds for k in axes(ψk,3), j in 1:G.ny, i in 1:G.nx
        ikx = im * G.kx[i]
        iky = im * G.ky[j]
        uk[i,j,k] = -iky * ψk[i,j,k]
        vk[i,j,k] =  ikx * ψk[i,j,k]
    end
    # Inverse FFT to real space
    if plans === nothing
        plans = plan_transforms!(G)
    end
    tmpu = similar(ψk)
    tmpv = similar(ψk)
    fft_backward!(tmpu, uk, plans)
    fft_backward!(tmpv, vk, plans)
    # Normalization: FFTW.ifft returns unnormalized; divide by (nx*ny)
    norm = (G.nx * G.ny)
    @inbounds for k in axes(S.u,3)
        # Real part
        S.u[:,:,k] .= real.(tmpu[:,:,k]) ./ norm
        S.v[:,:,k] .= real.(tmpv[:,:,k]) ./ norm
    end
    
    # Compute vertical velocity if requested
    if compute_w
        compute_vertical_velocity!(S, G, plans, params)
    else
        # Set w to zero (leading-order QG approximation)
        fill!(S.w, 0.0)
    end
    
    return S
end

"""
    compute_vertical_velocity!(S, G, plans, params)

Compute QG ageostrophic vertical velocity by solving the omega equation:
∇²w + (N²/f²)(∂²w/∂z²) = 2 J(ψ_z, ∇²ψ)

This is the diagnostic vertical velocity from quasi-geostrophic theory.
"""
function compute_vertical_velocity!(S::State, G::Grid, plans, params)
    nx, ny, nz = G.nx, G.ny, G.nz
    
    # For simple implementation, use existing omega_eqn_rhs! to get RHS
    # Then solve the elliptic equation for w
    using ..QGYBJ: omega_eqn_rhs!
    
    rhsk = similar(S.psi)
    omega_eqn_rhs!(rhsk, S.psi, G, plans)
    
    # Solve ∇²w + (N²/f²)(∂²w/∂z²) = RHS
    # For now, implement a simplified version focusing on horizontal structure
    wk = similar(S.psi)
    
    # Get stratification parameters
    if params !== nothing && hasfield(typeof(params), :Bu)
        Bu = params.Bu  # Burger number = N²H²/(f²L²)
        Ro = params.Ro  # Rossby number
        Fr = params.Fr  # Froude number
    else
        # Default values if params not provided
        Bu = 1.0
        Ro = 0.1
        Fr = 0.1
    end
    
    dz = nz > 1 ? (G.z[2] - G.z[1]) : 1.0
    
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        kh2 = G.kh2[i,j]
        
        if kh2 > 0  # Avoid division by zero at k=0
            # Simplified omega equation solver
            # ∇²w ≈ RHS (neglecting vertical structure for now)
            wk[i,j,k] = -rhsk[i,j,k] / kh2
        else
            wk[i,j,k] = 0.0
        end
    end
    
    # Apply boundary conditions: w = 0 at top and bottom
    if nz > 1
        wk[:,:,1] .= 0.0    # Bottom
        wk[:,:,nz] .= 0.0   # Top
    end
    
    # Transform to real space
    tmpw = similar(wk)
    fft_backward!(tmpw, wk, plans)
    
    # Store in state (real part, normalized)
    norm = nx * ny
    @inbounds for k in 1:nz
        S.w[:,:,k] .= real.(tmpw[:,:,k]) ./ norm
    end
    
    return S
end

end # module

using .Operators: compute_velocities!

