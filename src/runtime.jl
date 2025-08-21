"""
Convenience setup helpers to bootstrap a simulation.
"""

"""
    setup_model(; par=default_params()) -> (G, S, plans)

Initialize grid, state and FFT plans for a basic run.
"""
function setup_model(; par=default_params())
    G = init_grid(par)
    S = init_state(G)
    plans = plan_transforms!(G)
    a = a_ell_ut(par, G)
    return G, S, plans, a
end
