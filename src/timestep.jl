"""
Time stepping skeleton for QG–YBJ (projection + leapfrog). This is a placeholder
focused on wiring MPI/Pencil-compatible data flow; the physics terms (YBJ+ and
QG nonlinearities, diffusion, and dealiasing) should be ported next.
"""

"""
    first_projection_step!(S, G, par, plans; a)

Projection step (Forward Euler) example: given q^n, invert to ψ^n, compute
diagnostics. Extend to compute B and apply LA->A inversion similarly to Fortran.
"""
function first_projection_step!(S::State, G::Grid, par::QGParams, plans; a)
    invert_q_to_psi!(S, G; a)
    compute_velocities!(S, G; plans)
    return S
end

"""
    leapfrog_step!(Snp1, Sn, Snm1, G, par, plans; a)

Skeleton for leapfrog advance; fill with actual tendencies from Fortran model.
Here just carries ψ diagnostic refresh.
"""
function leapfrog_step!(Snp1::State, Sn::State, Snm1::State,
                        G::Grid, par::QGParams, plans; a)
    # Placeholder: copy spectral q forward (no dynamics) and recompute ψ
    Snp1.q .= Sn.q
    invert_q_to_psi!(Snp1, G; a)
    compute_velocities!(Snp1, G; plans)
    return Snp1
end

