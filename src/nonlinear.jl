"""
Nonlinear and linear tendency computations: Jacobians, refraction placeholder,
and hyperdiffusion factors with 2/3 dealiasing support.
"""

module Nonlinear

using ..QGYBJ: Grid
using ..QGYBJ: plan_transforms!, fft_forward!, fft_backward!

"""
    jacobian_spectral!(dstk, phik, chik, G, plans)

Compute J(phi, chi) = phi_x chi_y - phi_y chi_x using spectral derivatives and
2D transforms per z-slab. Writes result in spectral space dstk.
"""
function jacobian_spectral!(dstk, phik, chik, G::Grid, plans)
    nx, ny, nz = G.nx, G.ny, G.nz
    # spectral derivatives
    phixk = similar(phik); phiyk = similar(phik)
    chixk = similar(chik); chiyk = similar(chik)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        phixk[i,j,k] = im*G.kx[i]*phik[i,j,k]
        phiyk[i,j,k] = im*G.ky[j]*phik[i,j,k]
        chixk[i,j,k] = im*G.kx[i]*chik[i,j,k]
        chiyk[i,j,k] = im*G.ky[j]*chik[i,j,k]
    end
    # inverse to real space
    phix = similar(phik); phiy = similar(phik)
    chix = similar(chik); chiy = similar(chik)
    fft_backward!(phix, phixk, plans)
    fft_backward!(phiy, phiyk, plans)
    fft_backward!(chix, chixk, plans)
    fft_backward!(chiy, chiyk, plans)
    # form J in real space
    J = similar(phik)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        J[i,j,k] = (real(phix[i,j,k])*real(chiy[i,j,k]) - real(phiy[i,j,k])*real(chix[i,j,k]))
    end
    # forward to spectral
    fft_forward!(dstk, J, plans)
    # normalization of FFTs (inverse wasn’t normalized)
    norm = nx*ny
    @inbounds dstk .= dstk ./ norm
    return dstk
end

"""
    hyperdiffusion_factor(kx, ky, par; waves=false)

Compute scalar integrating factor argument used in leapfrog updates.
"""
function hyperdiffusion_factor(kx::Real, ky::Real, par; waves::Bool=false)
    if waves
        return par.dt * ( par.nu_h * ((abs(kx))^(2) + (abs(ky))^(2)) )
    else
        return par.dt * ( par.nu_h * ((abs(kx))^(2) + (abs(ky))^(2)) )
    end
end

end # module

using .Nonlinear: jacobian_spectral!, hyperdiffusion_factor

