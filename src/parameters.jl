"""
    QGParams

Container for physical and numerical parameters of the QG–YBJ model.
This is a simplified starting point; extend as needed while porting.
"""
Base.@kwdef struct QGParams{T}
    # Domain sizes
    nx::Int
    ny::Int
    nz::Int
    Lx::T
    Ly::T

    # Time stepping
    dt::T
    nt::Int

    # Physical parameters (nondimensionalize as appropriate)
    Ro::T              # Rossby number
    Fr::T              # Froude number
    f0::T              # Coriolis parameter

    # Viscosity/hyperviscosity
    nu_h::T            # horizontal viscosity
    nu_v::T            # vertical viscosity

    # Flags
    linear_vert_structure::Int  # mapping from Fortran param
end

"""
    default_params(; kwargs...)

Construct a reasonable default parameter set for experimentation.
"""
function default_params(; nx=64, ny=64, nz=64, Lx=1.0, Ly=1.0,
                           dt=1e-3, nt=10_000, Ro=0.1, Fr=0.1, f0=1.0,
                           nu_h=0.0, nu_v=0.0, linear_vert_structure=0)
    return QGParams{Float64}(; nx, ny, nz, Lx, Ly, dt, nt, Ro, Fr, f0, nu_h, nu_v, linear_vert_structure)
end

