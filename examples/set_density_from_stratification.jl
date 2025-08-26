"""
Example: derive vertical density-like profiles from a stratification and run
one ψ inversion using density-weighted operators.

This serves as a template — the mapping from N² to `rho_ut`/`rho_st` here is a
placeholder (unity). Replace with your preferred model-based mapping.
"""

using QGYBJ

function main()
    # Build params/grid
    par = default_params(nx=32, ny=32, nz=16, stratification=:constant_N)
    G, S, plans, a = setup_model(; par)

    # Compute N² profile and derive density-like profiles
    using .QGYBJ: N2_ut, derive_density_profiles
    N2 = N2_ut(par, G)
    rho_ut_profile, rho_st_profile = derive_density_profiles(par, G; N2_profile=N2)
    par2 = with_density_profiles(par; rho_ut=rho_ut_profile, rho_st=rho_st_profile)

    # Fill a simple q mode and invert to ψ with density-weighted operator
    S.q[3,3,8] = 1 + 0im
    invert_q_to_psi!(S, G; a=a_ell_ut(par2, G), par=par2)
    println("L2(||ψ||) = ", sqrt(sum(abs2, S.psi)))
end

abspath(PROGRAM_FILE) == @__FILE__ && main()
