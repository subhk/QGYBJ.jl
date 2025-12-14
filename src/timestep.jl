#=
================================================================================
                    timestep.jl - Time Integration
================================================================================

This file implements the time stepping schemes for the QG-YBJ+ model:

1. FORWARD EULER (Projection Step):
   Used for the first time step to initialize the leapfrog scheme.
   Simple but only first-order accurate.

2. LEAPFROG with ROBERT-ASSELIN FILTER:
   The main time integration scheme. Second-order accurate in time.
   The Robert-Asselin filter damps the computational mode that can
   grow with leapfrog schemes.

TIME INTEGRATION ALGORITHM:
---------------------------
For each time step from n to n+1:

1. Compute tendencies at time n:
   - Advection: J(ψ, q), J(ψ, B)
   - Refraction: B × ζ
   - Vertical diffusion: νz ∂²q/∂z²

2. Apply physics switches:
   - linear: zero nonlinear terms
   - inviscid: zero dissipation
   - passive_scalar: zero dispersion and refraction
   - fixed_flow: zero q tendency

3. Time step with hyperdiffusion integrating factors:
   - Leapfrog: φ^(n+1) = φ^(n-1) × e^(-2λdt) - 2dt × tendency × e^(-λdt)
   - Forward Euler: φ^(n+1) = φ^n × e^(-λdt) - dt × tendency

4. Robert-Asselin filter (leapfrog only):
   φ̃^n = φ^n + γ(φ^(n-1) - 2φ^n + φ^(n+1))
   where γ ~ 10⁻³ is small to minimize damping

5. Wave feedback on mean flow:
   q* = q - qʷ (if wave feedback is enabled)

6. Diagnostic updates:
   - Invert q → ψ
   - Invert B → A (YBJ+) or compute A directly (normal YBJ)
   - Compute velocities from ψ

FORTRAN CORRESPONDENCE:
----------------------
The time stepping matches main_waqg.f90 in the Fortran QG_YBJp code.
The integrating factor approach for hyperdiffusion is from the Fortran.

STABILITY:
----------
- Leapfrog CFL condition: dt < min(dx/|u|, dy/|v|)
- Hyperdiffusion CFL: dt × ν × k^(2n) < 1 for largest k
- Robert-Asselin γ too large → excessive damping
- Robert-Asselin γ too small → computational mode growth

================================================================================
=#

#=
================================================================================
                    FORWARD EULER (Projection Step)
================================================================================
The projection step initializes the leapfrog scheme by providing values
at times n and n-1 from a single initial condition.

After this step:
- S contains fields at time n+1 (after one Euler step)
- The original S becomes the n-1 state for leapfrog
================================================================================
=#

"""
    first_projection_step!(S, G, par, plans; a, dealias_mask=nothing)

Forward Euler initialization step for the leapfrog time stepper.

# Purpose
The leapfrog scheme requires values at two time levels (n and n-1).
This function takes the initial state and advances it by one Forward Euler
step, providing the needed second time level.

# Algorithm
1. **Compute tendencies at time n:**
   - Advection of q and B by geostrophic flow
   - Wave refraction by vorticity
   - Vertical diffusion

2. **Apply physics switches:**
   - `linear`: Zero nonlinear advection
   - `inviscid`: Zero dissipation
   - `passive_scalar`: Zero dispersion and refraction
   - `fixed_flow`: Mean flow doesn't evolve

3. **Forward Euler update:**
   For each spectral mode:
   ```
   q^(n+1) = [q^n - dt × tendency_q + dt × diffusion] × exp(-λ_q × dt)
   B^(n+1) = [B^n - dt × tendency_B] × exp(-λ_B × dt)
   ```
   where λ is the hyperdiffusion factor.

4. **Wave feedback (optional):**
   ```
   q* = q - qʷ
   ```

5. **Diagnostic inversions:**
   - q → ψ (elliptic inversion)
   - B → A, C (YBJ+ inversion)
   - ψ → u, v (velocity computation)

# Arguments
- `S::State`: State to advance (modified in place)
- `G::Grid`: Grid struct
- `par::QGParams`: Model parameters
- `plans`: FFT plans
- `a`: Elliptic coefficient array a_ell(z) = Bu/N²
- `dealias_mask`: Optional 2/3 dealiasing mask (nx × ny)

# Returns
Modified state S at time n+1.

# Fortran Correspondence
This matches the projection step in main_waqg.f90.

# Example
```julia
# Initialize and run projection step
state = init_state(grid, params)
init_random_psi!(state, grid, params, plans)
a = a_ell_ut(params, grid)
L = dealias_mask(params, grid)
first_projection_step!(state, grid, params, plans; a=a, dealias_mask=L)
```
"""
function first_projection_step!(S::State, G::Grid, par::QGParams, plans; a, dealias_mask=nothing)
    #= Setup =#
    L = isnothing(dealias_mask) ? trues(G.nx,G.ny) : dealias_mask

    # Allocate tendency arrays
    nqk  = similar(S.q)    # J(ψ, q) advection of PV
    nBRk = similar(S.B)    # J(ψ, BR) advection of wave real part
    nBIk = similar(S.B)    # J(ψ, BI) advection of wave imaginary part
    rBRk = similar(S.B)    # BR × ζ refraction real part
    rBIk = similar(S.B)    # BI × ζ refraction imaginary part
    dqk  = similar(S.B)    # Vertical diffusion of q

    #= Split B into real and imaginary parts for computation
    The wave field B is complex; we work with BR = Re(B), BI = Im(B) =#
    BRk = similar(S.B); BIk = similar(S.B)
    @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
        BRk[i,j,k] = Complex(real(S.B[i,j,k]), 0)
        BIk[i,j,k] = Complex(imag(S.B[i,j,k]), 0)
    end

    #= Step 1: Compute diagnostic fields ψ and velocities =#
    invert_q_to_psi!(S, G; a, par=par)           # q → ψ
    compute_velocities!(S, G; plans, params=par) # ψ → u, v

    #= Step 2: Compute nonlinear tendencies =#

    # Advection: J(ψ, q), J(ψ, BR), J(ψ, BI)
    convol_waqg!(nqk, nBRk, nBIk, S.u, S.v, S.q, BRk, BIk, G, plans; Lmask=L)

    # Wave refraction: B × ζ where ζ = ∇²ψ
    refraction_waqg!(rBRk, rBIk, BRk, BIk, S.psi, G, plans; Lmask=L)

    # Vertical diffusion: νz ∂²q/∂z²
    dissipation_q_nv!(dqk, S.q, par, G)

    #= Step 3: Apply physics switches =#

    # inviscid: No dissipation
    if par.inviscid; dqk .= 0; end

    # linear: No nonlinear advection
    if par.linear
        nqk .= 0; nBRk .= 0; nBIk .= 0
    end

    # no_dispersion: Waves don't disperse (A = 0)
    if par.no_dispersion
        S.A .= 0; S.C .= 0
    end

    # passive_scalar: Waves are passive tracers (no dispersion, no refraction)
    if par.passive_scalar
        S.A .= 0; S.C .= 0; rBRk .= 0; rBIk .= 0
    end

    # fixed_flow: Mean flow doesn't evolve
    if par.fixed_flow; nqk .= 0; end

    #= Store old fields for time stepping =#
    qok  = copy(S.q)
    BRok = similar(S.B); BIok = similar(S.B)
    @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
        BRok[i,j,k] = Complex(real(S.B[i,j,k]), 0)
        BIok[i,j,k] = Complex(imag(S.B[i,j,k]), 0)
    end

    #= Step 4: Forward Euler with integrating factors =#
    # The integrating factor handles hyperdiffusion exactly:
    # φ^(n+1) = [φ^n - dt × F] × exp(-λ×dt)

    @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
        if L[i,j]
            kx = G.kx[i]; ky = G.ky[j]; kh2 = G.kh2[i,j]

            # Integrating factors for hyperdiffusion
            If = int_factor(kx, ky, par; waves=false)   # For mean flow
            Ifw = int_factor(kx, ky, par; waves=true)   # For waves

            #= Update q (QGPV) =#
            if par.fixed_flow
                # Keep q unchanged when mean flow is fixed
                S.q[i,j,k] = qok[i,j,k]
            else
                # q^(n+1) = [q^n - dt×J(ψ,q) + dt×diffusion] × exp(-λ×dt)
                S.q[i,j,k] = ( qok[i,j,k] - par.dt*nqk[i,j,k] + par.dt*dqk[i,j,k] ) * exp(-If)
            end

            #= Update B (wave envelope)
            The YBJ+ equation for B is:
                ∂B/∂t + J(ψ,B) = -(ikₕ²/2)A + (1/2)B×ζ
            In terms of real/imaginary parts:
                ∂BR/∂t = -J(ψ,BR) - (kₕ²/2)AI + (1/2)BI×ζ
                ∂BI/∂t = -J(ψ,BI) + (kₕ²/2)AR - (1/2)BR×ζ =#
            BRnew = ( BRok[i,j,k] - par.dt*nBRk[i,j,k]
                      - par.dt*0.5*kh2*Complex(imag(S.A[i,j,k]),0)
                      + par.dt*0.5*rBIk[i,j,k] ) * exp(-Ifw)
            BInew = ( BIok[i,j,k] - par.dt*nBIk[i,j,k]
                      + par.dt*0.5*kh2*Complex(real(S.A[i,j,k]),0)
                      - par.dt*0.5*rBRk[i,j,k] ) * exp(-Ifw)

            # Recombine into complex B
            S.B[i,j,k] = Complex(real(BRnew), 0) + im*Complex(real(BInew), 0)
        else
            # Zero out dealiased modes
            S.q[i,j,k] = 0
            S.B[i,j,k] = 0
        end
    end

    #= Step 5: Wave feedback on mean flow =#
    # q* = q - qʷ where qʷ is the wave feedback term
    wave_feedback_enabled = !par.no_feedback && !par.no_wave_feedback
    if wave_feedback_enabled
        qwk = similar(S.q)

        # Rebuild BR/BI from updated B
        @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
            BRk[i,j,k] = Complex(real(S.B[i,j,k]), 0)
            BIk[i,j,k] = Complex(imag(S.B[i,j,k]), 0)
        end

        # Compute qʷ from B
        compute_qw!(qwk, BRk, BIk, par, G, plans; Lmask=L)

        # Subtract from q
        @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
            if L[i,j]
                S.q[i,j,k] -= qwk[i,j,k]
            else
                S.q[i,j,k] = 0
            end
        end
    end

    #= Step 6: Update diagnostic fields =#

    # Invert q → ψ (only if mean flow evolves)
    if !par.fixed_flow
        invert_q_to_psi!(S, G; a, par=par)
    end

    # Recover A from B
    if par.ybj_plus
        # YBJ+: Solve elliptic problem B → A, C
        invert_B_to_A!(S, G, par, a)
    else
        # Normal YBJ: Different procedure
        BRk2 = similar(S.B); BIk2 = similar(S.B)
        @inbounds for k in 1:G.nz, j in 1:G.ny, i in 1:G.nx
            BRk2[i,j,k] = Complex(real(S.B[i,j,k]), 0)
            BIk2[i,j,k] = Complex(imag(S.B[i,j,k]), 0)
        end
        sumB!(S.B, G; Lmask=L)  # Remove vertical mean
        sigma = compute_sigma(par, G, nBRk, nBIk, rBRk, rBIk; Lmask=L)
        compute_A!(S.A, S.C, BRk2, BIk2, sigma, par, G; Lmask=L)
    end

    # Compute velocities from ψ
    compute_velocities!(S, G; plans, params=par)

    return S
end

#=
================================================================================
                    LEAPFROG with ROBERT-ASSELIN FILTER
================================================================================
The leapfrog scheme is:
    φ^(n+1) = φ^(n-1) + 2dt × F^n

This is second-order accurate but has a computational mode that can grow.
The Robert-Asselin filter damps this mode:
    φ̃^n = φ^n + γ(φ^(n-1) - 2φ^n + φ^(n+1))

With the integrating factor for hyperdiffusion:
    φ^(n+1) = φ^(n-1) × e^(-2λdt) + 2dt × F^n × e^(-λdt)

This ensures exact treatment of the linear diffusion terms.
================================================================================
=#

"""
    leapfrog_step!(Snp1, Sn, Snm1, G, par, plans; a, dealias_mask=nothing)

Advance the solution by one leapfrog time step with Robert-Asselin filtering.

# Algorithm

**1. Compute tendencies at time n:**
```
F_q^n = J(ψ^n, q^n) - νz∂²q^(n-1)/∂z²
F_B^n = J(ψ^n, B^n) + dispersion + refraction
```

**2. Leapfrog update with integrating factors:**
For each spectral mode (k):
```
q^(n+1) = q^(n-1) × e^(-2λdt) - 2dt × F_q^n × e^(-λdt) + 2dt × diff^(n-1) × e^(-2λdt)
B^(n+1) = B^(n-1) × e^(-2λdt) - 2dt × F_B^n × e^(-λdt)
```

**3. Robert-Asselin filter:**
```
q̃^n = q^n + γ(q^(n-1) - 2q^n + q^(n+1))
B̃^n = B^n + γ(B^(n-1) - 2B^n + B^(n+1))
```
The filtered values are stored in Snm1 for the next step.

**4. Wave feedback (if enabled):**
```
q*^(n+1) = q^(n+1) - qʷ^(n+1)
```

**5. Diagnostic inversions:**
- q^(n+1) → ψ^(n+1)
- B^(n+1) → A^(n+1), C^(n+1)
- ψ^(n+1) → u^(n+1), v^(n+1)

# Arguments
- `Snp1::State`: State at time n+1 (output)
- `Sn::State`: State at time n (input, filter applied to Snm1)
- `Snm1::State`: State at time n-1 (input, receives filtered values)
- `G::Grid`: Grid struct
- `par::QGParams`: Model parameters
- `plans`: FFT plans
- `a`: Elliptic coefficient array
- `dealias_mask`: Optional dealiasing mask

# Returns
Modified Snp1 with solution at time n+1.

# Time Level Management
After this call:
- Snp1 contains fields at n+1
- Snm1 contains **filtered** fields at n (for next step's n-1)
- Sn is unchanged (use Snm1 as new Sn in next call)

Typical loop structure:
```julia
for iter in 1:nsteps
    leapfrog_step!(Snp1, Sn, Snm1, G, par, plans; a=a)
    # Rotate: Snm1 ← Sn, Sn ← Snp1
    Snm1, Sn, Snp1 = Sn, Snp1, Snm1
end
```

# Fortran Correspondence
This matches the main leapfrog loop in main_waqg.f90.

# Example
```julia
# After projection step, run leapfrog
for iter in 1:1000
    leapfrog_step!(Snp1, Sn, Snm1, grid, params, plans; a=a, dealias_mask=L)
    Snm1, Sn, Snp1 = Sn, Snp1, Snm1
end
```
"""
function leapfrog_step!(Snp1::State, Sn::State, Snm1::State,
                        G::Grid, par::QGParams, plans; a, dealias_mask=nothing)
    nx, ny, nz = G.nx, G.ny, G.nz
    L = isnothing(dealias_mask) ? trues(nx,ny) : dealias_mask

    #= Step 1: Update diagnostics for current state =#
    if !par.fixed_flow
        invert_q_to_psi!(Sn, G; a, par=par)
    end
    compute_velocities!(Sn, G; plans, params=par)

    #= Step 2: Allocate and compute tendencies =#
    nqk  = similar(Sn.q)    # Advection of q
    nBRk = similar(Sn.B)    # Advection of BR
    nBIk = similar(Sn.B)    # Advection of BI
    rBRk = similar(Sn.B)    # Refraction of BR
    rBIk = similar(Sn.B)    # Refraction of BI
    dqk  = similar(Sn.B)    # Vertical diffusion

    # Split B into real/imaginary
    BRk = similar(Sn.B); BIk = similar(Sn.B)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        BRk[i,j,k] = Complex(real(Sn.B[i,j,k]), 0)
        BIk[i,j,k] = Complex(imag(Sn.B[i,j,k]), 0)
    end

    # Compute tendencies
    convol_waqg!(nqk, nBRk, nBIk, Sn.u, Sn.v, Sn.q, BRk, BIk, G, plans; Lmask=L)
    refraction_waqg!(rBRk, rBIk, BRk, BIk, Sn.psi, G, plans; Lmask=L)

    # Vertical diffusion uses q at n-1 (for leapfrog stability)
    dissipation_q_nv!(dqk, Snm1.q, par, G)

    #= Step 3: Apply physics switches =#
    if par.inviscid; dqk .= 0; end
    if par.linear; nqk .= 0; nBRk .= 0; nBIk .= 0; end
    if par.no_dispersion; Sn.A .= 0; end
    if par.passive_scalar; Sn.A .= 0; rBRk .= 0; rBIk .= 0; end
    if par.fixed_flow; nqk .= 0; end

    #= Step 4: Leapfrog update with integrating factors =#
    qtemp = similar(Sn.q)
    BRtemp = similar(Sn.B); BItemp = similar(Sn.B)

    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        if L[i,j]
            kx = G.kx[i]; ky = G.ky[j]; kh2 = G.kh2[i,j]
            If  = int_factor(kx, ky, par; waves=false)
            Ifw = int_factor(kx, ky, par; waves=true)

            #= Update q
            q^(n+1) = q^(n-1)×e^(-2λdt) - 2dt×J(ψ,q)^n×e^(-λdt) + 2dt×diff^(n-1)×e^(-2λdt) =#
            if par.fixed_flow
                qtemp[i,j,k] = Sn.q[i,j,k]  # Keep unchanged
            else
                qtemp[i,j,k] = Snm1.q[i,j,k]*exp(-2If) -
                               2*par.dt*nqk[i,j,k]*exp(-If) +
                               2*par.dt*dqk[i,j,k]*exp(-2If)
            end

            #= Update B (real and imaginary parts)
            BR^(n+1) = BR^(n-1)×e^(-2λdt) - 2dt×[J(ψ,BR) + (kₕ²/2)AI - (1/2)BI×ζ]×e^(-λdt)
            BI^(n+1) = BI^(n-1)×e^(-2λdt) - 2dt×[J(ψ,BI) - (kₕ²/2)AR + (1/2)BR×ζ]×e^(-λdt) =#
            BRtemp[i,j,k] = Complex(real(Snm1.B[i,j,k]),0)*exp(-2Ifw) -
                           2*par.dt*( nBRk[i,j,k] +
                                     0.5*kh2*Complex(imag(Sn.A[i,j,k]),0) -
                                     0.5*rBIk[i,j,k] )*exp(-Ifw)
            BItemp[i,j,k] = Complex(imag(Snm1.B[i,j,k]),0)*exp(-2Ifw) -
                           2*par.dt*( nBIk[i,j,k] -
                                     0.5*kh2*Complex(real(Sn.A[i,j,k]),0) +
                                     0.5*rBRk[i,j,k] )*exp(-Ifw)
        else
            qtemp[i,j,k] = 0; BRtemp[i,j,k] = 0; BItemp[i,j,k] = 0
        end
    end

    #= Step 5: Robert-Asselin filter
    Damps the computational mode: φ̃^n = φ^n + γ(φ^(n-1) - 2φ^n + φ^(n+1))
    Store filtered values in Snm1 (they become the "n-1" for next step) =#
    γ = par.gamma
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        if L[i,j]
            # Filter q
            Snm1.q[i,j,k] = Sn.q[i,j,k] + γ*( Snm1.q[i,j,k] - 2Sn.q[i,j,k] + qtemp[i,j,k] )

            # Filter B
            Bnp1 = Complex(real(BRtemp[i,j,k]),0) + im*Complex(real(BItemp[i,j,k]),0)
            Snm1.B[i,j,k] = Sn.B[i,j,k] + γ*( Snm1.B[i,j,k] - 2Sn.B[i,j,k] + Bnp1 )
        else
            Snm1.q[i,j,k] = 0; Snm1.B[i,j,k] = 0
        end
    end

    #= Step 6: Accept the new solution =#
    Snp1.q .= qtemp
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        Snp1.B[i,j,k] = Complex(real(BRtemp[i,j,k]),0) + im*Complex(real(BItemp[i,j,k]),0)
    end

    #= Step 7: Wave feedback on mean flow =#
    wave_feedback_enabled = !par.no_feedback && !par.no_wave_feedback
    if wave_feedback_enabled
        qwk = similar(Snp1.q)

        # Rebuild BR/BI from updated B
        BRk2 = similar(Snp1.B); BIk2 = similar(Snp1.B)
        @inbounds for kk in 1:nz, jj in 1:ny, ii in 1:nx
            BRk2[ii,jj,kk] = Complex(real(Snp1.B[ii,jj,kk]),0)
            BIk2[ii,jj,kk] = Complex(imag(Snp1.B[ii,jj,kk]),0)
        end

        compute_qw!(qwk, BRk2, BIk2, par, G, plans; Lmask=L)

        @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
            if L[i,j]
                Snp1.q[i,j,k] -= qwk[i,j,k]
            else
                Snp1.q[i,j,k] = 0
            end
        end
    end

    #= Step 8: Update diagnostics for new state =#

    # Invert q → ψ
    if !par.fixed_flow
        invert_q_to_psi!(Snp1, G; a, par=par)
    end

    # Recover A from B
    if par.ybj_plus
        invert_B_to_A!(Snp1, G, par, a)
    else
        # Normal YBJ path
        BRk3 = similar(Snp1.B); BIk3 = similar(Snp1.B)
        @inbounds for kk in 1:nz, jj in 1:ny, ii in 1:nx
            BRk3[ii,jj,kk] = Complex(real(Snp1.B[ii,jj,kk]), 0)
            BIk3[ii,jj,kk] = Complex(imag(Snp1.B[ii,jj,kk]), 0)
        end
        sumB!(Snp1.B, G; Lmask=L)
        sigma2 = compute_sigma(par, G, nBRk, nBIk, rBRk, rBIk; Lmask=L)
        compute_A!(Snp1.A, Snp1.C, BRk3, BIk3, sigma2, par, G; Lmask=L)
    end

    # Compute velocities
    compute_velocities!(Snp1, G; plans, params=par)

    return Snp1
end
