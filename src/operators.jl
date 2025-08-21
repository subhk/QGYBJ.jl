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

Given spectral streamfunction `psi(kx,ky,z)`, compute horizontal velocities
`u = -∂ψ/∂y`, `v =  ∂ψ/∂x` in real space. Uses spectral differentiation.
"""
function compute_velocities!(S::State, G::Grid; plans=nothing)
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
    return S
end

end # module

using .Operators: compute_velocities!

