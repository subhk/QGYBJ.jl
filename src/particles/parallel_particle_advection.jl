"""
Parallel-aware particle advection for distributed QG-YBJ simulations.

This module extends the basic particle advection to work with MPI domain decomposition
using PencilArrays, handling particle migration between domains and distributed I/O.
"""

module ParallelParticleAdvection

using ..QGYBJ: Grid, State, QGParams
using ..ParticleAdvection: ParticleConfig, ParticleState, ParticleTracker, 
                          initialize_particles!, interpolate_velocity,
                          apply_boundary_conditions!, save_particle_state!

export ParallelParticleTracker, setup_parallel_particles!, 
       advect_particles_parallel!, migrate_particles!,
       gather_particles_for_io, write_parallel_particle_trajectories

"""
    ParallelParticleTracker

Particle tracker that handles MPI domain decomposition.
"""
mutable struct ParallelParticleTracker{T<:AbstractFloat}
    # Basic particle tracker
    base_tracker::ParticleTracker{T}
    
    # Parallel information
    comm::Any  # MPI communicator
    rank::Int
    nprocs::Int
    
    # Domain decomposition info
    local_domain::NamedTuple  # (x_start, x_end, y_start, y_end, z_start, z_end)
    global_grid::Grid        # Full global grid
    
    # Particle migration buffers
    send_particles::Vector{Vector{T}}  # Particles to send to each rank
    recv_particles::Vector{Vector{T}}  # Particles received from each rank
    
    # I/O configuration
    is_io_rank::Bool         # True if this rank handles I/O
    gather_for_io::Bool      # True if gathering particles for I/O
end

"""
    setup_parallel_particles!(config, grid, parallel_config)

Set up particle tracking for parallel execution.
"""
function setup_parallel_particles!(particle_config::ParticleConfig{T}, 
                                   grid::Grid, 
                                   parallel_config) where T
    
    # Try to get MPI info
    comm = nothing
    rank = 0
    nprocs = 1
    
    try
        import MPI
        if MPI.Initialized()
            comm = MPI.COMM_WORLD
            rank = MPI.Comm_rank(comm)
            nprocs = MPI.Comm_size(comm)
        end
    catch
        @warn "MPI not available, using serial particle advection"
    end
    
    # Create base tracker
    base_tracker = ParticleTracker(particle_config, grid)
    
    # Determine local domain bounds (simplified - would use actual PencilArrays info)
    local_domain = compute_local_domain(grid, rank, nprocs)
    
    # Initialize particles only in local domain
    initialize_local_particles!(base_tracker, particle_config, local_domain)
    
    # Create parallel tracker
    parallel_tracker = ParallelParticleTracker{T}(
        base_tracker,
        comm, rank, nprocs,
        local_domain, grid,
        [T[] for _ in 1:nprocs],  # send buffers
        [T[] for _ in 1:nprocs],  # recv buffers
        rank == 0,               # rank 0 is I/O rank
        true                     # gather for I/O by default
    )
    
    return parallel_tracker
end

"""
    compute_local_domain(grid, rank, nprocs)

Compute local domain bounds for this MPI rank (simplified 1D decomposition).
"""
function compute_local_domain(grid::Grid, rank::Int, nprocs::Int)
    # Simple 1D decomposition in x-direction
    nx_local = grid.nx ÷ nprocs
    remainder = grid.nx % nprocs
    
    if rank < remainder
        nx_local += 1
        x_start = rank * nx_local
    else
        x_start = remainder * (nx_local + 1) + (rank - remainder) * nx_local
    end
    
    x_end = x_start + nx_local - 1
    
    # Convert to physical coordinates
    dx = grid.Lx / grid.nx
    x_start_phys = x_start * dx
    x_end_phys = (x_end + 1) * dx
    
    return (x_start=x_start_phys, x_end=x_end_phys,
            y_start=0.0, y_end=grid.Ly,
            z_start=0.0, z_end=grid.Lz)
end

"""
    initialize_local_particles!(tracker, config, local_domain)

Initialize particles only within the local domain.
"""
function initialize_local_particles!(tracker::ParticleTracker{T}, 
                                    config::ParticleConfig{T},
                                    local_domain) where T
    # Find intersection of particle region with local domain
    x_min = max(config.x_min, local_domain.x_start)
    x_max = min(config.x_max, local_domain.x_end)
    
    if x_min >= x_max
        # No particles in this domain
        resize!(tracker.particles.x, 0)
        resize!(tracker.particles.y, 0)
        resize!(tracker.particles.z, 0)
        resize!(tracker.particles.u, 0)
        resize!(tracker.particles.v, 0)
        resize!(tracker.particles.w, 0)
        tracker.particles.np = 0
        return tracker
    end
    
    # Calculate local particle counts
    x_frac = (x_max - x_min) / (config.x_max - config.x_min)
    nx_local = max(1, round(Int, config.nx_particles * x_frac))
    
    # Initialize particles in local region
    local_config = ParticleConfig{T}(
        x_min=x_min, x_max=x_max,
        y_min=config.y_min, y_max=config.y_max,
        z_level=config.z_level,
        nx_particles=nx_local, ny_particles=config.ny_particles,
        use_ybj_w=config.use_ybj_w,
        use_3d_advection=config.use_3d_advection,
        integration_method=config.integration_method,
        periodic_x=config.periodic_x,
        periodic_y=config.periodic_y,
        reflect_z=config.reflect_z
    )
    
    initialize_particles!(tracker, local_config)
    
    return tracker
end

"""
    advect_particles_parallel!(parallel_tracker, state, grid, dt)

Advect particles in parallel with domain decomposition.
"""
function advect_particles_parallel!(parallel_tracker::ParallelParticleTracker{T},
                                   state::State, grid::Grid, dt::T) where T
    
    # Update local velocity fields
    update_velocity_fields!(parallel_tracker.base_tracker, state, grid, nothing, nothing)
    
    # Advect local particles
    advect_particles!(parallel_tracker.base_tracker, state, grid, dt)
    
    # Handle particle migration between domains
    if parallel_tracker.nprocs > 1
        migrate_particles!(parallel_tracker)
    end
    
    # Apply boundary conditions
    apply_boundary_conditions!(parallel_tracker.base_tracker)
    
    return parallel_tracker
end

"""
    migrate_particles!(parallel_tracker)

Handle particles that have moved outside local domain.
"""
function migrate_particles!(parallel_tracker::ParallelParticleTracker{T}) where T
    
    if parallel_tracker.comm === nothing
        return  # No MPI available
    end
    
    try
        import MPI
        
        particles = parallel_tracker.base_tracker.particles
        local_domain = parallel_tracker.local_domain
        
        # Clear send buffers
        for i in 1:parallel_tracker.nprocs
            empty!(parallel_tracker.send_particles[i])
        end
        
        # Find particles that need to migrate
        keep_indices = Int[]
        
        for i in 1:particles.np
            x = particles.x[i]
            
            # Determine which rank this particle belongs to
            target_rank = find_target_rank(x, parallel_tracker.global_grid, parallel_tracker.nprocs)
            
            if target_rank == parallel_tracker.rank
                # Keep this particle
                push!(keep_indices, i)
            else
                # Send to target rank
                particle_data = [particles.x[i], particles.y[i], particles.z[i],
                               particles.u[i], particles.v[i], particles.w[i]]
                append!(parallel_tracker.send_particles[target_rank + 1], particle_data)
            end
        end
        
        # Keep only local particles
        particles.x = particles.x[keep_indices]
        particles.y = particles.y[keep_indices]
        particles.z = particles.z[keep_indices]
        particles.u = particles.u[keep_indices]
        particles.v = particles.v[keep_indices]
        particles.w = particles.w[keep_indices]
        particles.np = length(keep_indices)
        
        # Exchange particles using MPI_Alltoall-like communication
        exchange_particles!(parallel_tracker)
        
    catch e
        @warn "Particle migration failed: $e"
    end
    
    return parallel_tracker
end

"""
    find_target_rank(x, grid, nprocs)

Find which MPI rank should own a particle at position x.
"""
function find_target_rank(x::T, grid::Grid, nprocs::Int) where T
    # Handle periodic boundaries
    x_periodic = mod(x, grid.Lx)
    
    # Simple 1D decomposition
    dx = grid.Lx / nprocs
    rank = min(nprocs - 1, floor(Int, x_periodic / dx))
    
    return rank
end

"""
    exchange_particles!(parallel_tracker)

Exchange particles between MPI ranks.
"""
function exchange_particles!(parallel_tracker::ParallelParticleTracker{T}) where T
    
    try
        import MPI
        
        comm = parallel_tracker.comm
        nprocs = parallel_tracker.nprocs
        
        # Clear receive buffers
        for i in 1:nprocs
            empty!(parallel_tracker.recv_particles[i])
        end
        
        # Send particle counts first
        send_counts = [length(parallel_tracker.send_particles[i]) ÷ 6 for i in 1:nprocs]
        recv_counts = MPI.Alltoall(send_counts, comm)
        
        # Send actual particle data
        for src_rank in 0:nprocs-1
            if src_rank == parallel_tracker.rank
                continue
            end
            
            # Send to src_rank
            if !isempty(parallel_tracker.send_particles[src_rank + 1])
                MPI.Send(parallel_tracker.send_particles[src_rank + 1], src_rank, 0, comm)
            end
            
            # Receive from src_rank
            if recv_counts[src_rank + 1] > 0
                recv_data = Vector{T}(undef, recv_counts[src_rank + 1] * 6)
                MPI.Recv!(recv_data, src_rank, 0, comm)
                parallel_tracker.recv_particles[src_rank + 1] = recv_data
            end
        end
        
        # Add received particles to local collection
        add_received_particles!(parallel_tracker)
        
    catch e
        @warn "Particle exchange failed: $e"
    end
end

"""
    add_received_particles!(parallel_tracker)

Add received particles to the local particle collection.
"""
function add_received_particles!(parallel_tracker::ParallelParticleTracker{T}) where T
    particles = parallel_tracker.base_tracker.particles
    
    for rank_data in parallel_tracker.recv_particles
        if !isempty(rank_data)
            n_new = length(rank_data) ÷ 6
            
            for i in 1:n_new
                idx = (i-1) * 6
                push!(particles.x, rank_data[idx + 1])
                push!(particles.y, rank_data[idx + 2])
                push!(particles.z, rank_data[idx + 3])
                push!(particles.u, rank_data[idx + 4])
                push!(particles.v, rank_data[idx + 5])
                push!(particles.w, rank_data[idx + 6])
            end
            
            particles.np += n_new
        end
    end
end

"""
    interpolate_velocity_parallel(x, y, z, tracker, state, grid)

Interpolate velocity accounting for domain decomposition and halo exchanges.
"""
function interpolate_velocity_parallel(x::T, y::T, z::T,
                                     parallel_tracker::ParallelParticleTracker{T},
                                     state::State, grid::Grid) where T
    
    # Check if particle is in local domain
    local_domain = parallel_tracker.local_domain
    
    if x >= local_domain.x_start && x <= local_domain.x_end
        # Particle is local - use standard interpolation
        return interpolate_velocity(x, y, z, 
                                  parallel_tracker.base_tracker.u_field,
                                  parallel_tracker.base_tracker.v_field,
                                  parallel_tracker.base_tracker.w_field,
                                  parallel_tracker.base_tracker)
    else
        # Particle needs data from another domain
        # For now, return zero velocity (would need halo exchange in full implementation)
        @warn "Particle outside local domain - using zero velocity"
        return (0.0, 0.0, 0.0)
    end
end

"""
    gather_particles_for_io(parallel_tracker)

Gather all particles to I/O rank for writing.
"""
function gather_particles_for_io(parallel_tracker::ParallelParticleTracker{T}) where T
    
    if parallel_tracker.comm === nothing || parallel_tracker.nprocs == 1
        return parallel_tracker.base_tracker  # Serial case
    end
    
    try
        import MPI
        
        comm = parallel_tracker.comm
        rank = parallel_tracker.rank
        
        # Gather particle counts
        local_count = parallel_tracker.base_tracker.particles.np
        all_counts = MPI.Gather(local_count, 0, comm)
        
        if rank == 0
            # Gather all particle data to rank 0
            total_particles = sum(all_counts)
            
            # Create combined tracker
            combined_tracker = deepcopy(parallel_tracker.base_tracker)
            
            # Resize arrays
            resize!(combined_tracker.particles.x, total_particles)
            resize!(combined_tracker.particles.y, total_particles)
            resize!(combined_tracker.particles.z, total_particles)
            resize!(combined_tracker.particles.u, total_particles)
            resize!(combined_tracker.particles.v, total_particles)
            resize!(combined_tracker.particles.w, total_particles)
            combined_tracker.particles.np = total_particles
            
            # Copy local particles first
            idx = 1
            local_particles = parallel_tracker.base_tracker.particles
            for i in 1:local_particles.np
                combined_tracker.particles.x[idx] = local_particles.x[i]
                combined_tracker.particles.y[idx] = local_particles.y[i]
                combined_tracker.particles.z[idx] = local_particles.z[i]
                combined_tracker.particles.u[idx] = local_particles.u[i]
                combined_tracker.particles.v[idx] = local_particles.v[i]
                combined_tracker.particles.w[idx] = local_particles.w[i]
                idx += 1
            end
            
            # Receive from other ranks
            for src_rank in 1:parallel_tracker.nprocs-1
                if all_counts[src_rank + 1] > 0
                    recv_data = Vector{T}(undef, all_counts[src_rank + 1] * 6)
                    MPI.Recv!(recv_data, src_rank, 1, comm)
                    
                    for i in 1:all_counts[src_rank + 1]
                        data_idx = (i-1) * 6
                        combined_tracker.particles.x[idx] = recv_data[data_idx + 1]
                        combined_tracker.particles.y[idx] = recv_data[data_idx + 2]
                        combined_tracker.particles.z[idx] = recv_data[data_idx + 3]
                        combined_tracker.particles.u[idx] = recv_data[data_idx + 4]
                        combined_tracker.particles.v[idx] = recv_data[data_idx + 5]
                        combined_tracker.particles.w[idx] = recv_data[data_idx + 6]
                        idx += 1
                    end
                end
            end
            
            return combined_tracker
        else
            # Send local particles to rank 0
            if local_count > 0
                send_data = T[]
                local_particles = parallel_tracker.base_tracker.particles
                for i in 1:local_particles.np
                    append!(send_data, [local_particles.x[i], local_particles.y[i], local_particles.z[i],
                                       local_particles.u[i], local_particles.v[i], local_particles.w[i]])
                end
                MPI.Send(send_data, 0, 1, comm)
            end
            
            return nothing  # Non-I/O ranks return nothing
        end
        
    catch e
        @warn "Particle gathering failed: $e"
        return parallel_tracker.base_tracker
    end
end

"""
    write_parallel_particle_trajectories(filename, parallel_tracker; metadata=Dict())

Write particle trajectories from parallel simulation.
"""
function write_parallel_particle_trajectories(filename::String, 
                                             parallel_tracker::ParallelParticleTracker{T};
                                             metadata::Dict=Dict()) where T
    
    if parallel_tracker.is_io_rank
        if parallel_tracker.gather_for_io
            # Gather all particles and write from I/O rank
            combined_tracker = gather_particles_for_io(parallel_tracker)
            if combined_tracker !== nothing
                write_particle_trajectories(filename, combined_tracker, metadata=metadata)
            end
        else
            # Write only local particles (separate files per rank)
            rank_filename = replace(filename, ".nc" => "_rank$(parallel_tracker.rank).nc")
            write_particle_trajectories(rank_filename, parallel_tracker.base_tracker, metadata=metadata)
        end
    end
    
    # Synchronize all ranks
    if parallel_tracker.comm !== nothing
        try
            import MPI
            MPI.Barrier(parallel_tracker.comm)
        catch
        end
    end
    
    return filename
end

end # module ParallelParticleAdvection

using .ParallelParticleAdvection