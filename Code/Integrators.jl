
# =============================================================================
# File: src/integrators.jl
# =============================================================================

"""
    velocity_verlet!(system::System, dt::Float64, force_func, potential_func)

Perform one velocity Verlet integration step.

# Algorithm
1. Update positions: r(t+dt) = r(t) + v(t)dt + 0.5*F(t)/m*dt²
2. Calculate new forces F(t+dt)
3. Update velocities: v(t+dt) = v(t) + 0.5*[F(t)+F(t+dt)]/m*dt
"""
function velocity_verlet!(state::SyState, system::System, dt::Float64)
    n = system.N_atoms
    
    # Save old forces
    old_forces = copy(state.force)
    
    # Update positions
    ## Vectorized code
    Mass_matrix = state.M
    state.Coords .+= state.Vel * dt .+ 0.5 * old_forces ./ Mass_matrix * dt^2
    # Calculate new forces
    calculate_forces!(state, system)
    # Update velocities using average force
    avg_f = 0.5 * (old_forces+ state.force)
    state.Vel .+= avg_f ./ Mass_matrix * dt
    # Total energy update after velocity verlet
    state.TotE = 0.5*sum(state.M .* state.Vel.^2) + state.Ene
    # Update time
    state.t += dt
    return nothing
end

"""
    apply_boundary_conditions!(system::System, boundary_type=:reflective)

Apply boundary conditions to system.

# Boundary types
- `:reflective`: Particles bounce off walls (velocity reversal)
- `:periodic`: Periodic boundary conditions
- `:none`: No boundaries
"""
function apply_boundary_conditions!(State::SyState,system::System, boundary_type=:reflective)
    if boundary_type == :none
        return nothing
    end
    
    bl = system.box_size
    
    if boundary_type == :reflective
        # Check each particle
        for i in 1:size(State.Coords, 1)
            for dim in 1:3
                # Check if particle is outside boundary
                if abs(State.Coords[i, dim]) >= bl
                    # Reverse velocity component
                    State.Vel[i, dim] *= -1.0
                    # Clamp position to boundary
                    State.Coords[i, dim] = sign(State.Coords[i, dim]) * (bl - 1e-6)
                end
            end
        end
    elseif boundary_type == :periodic
        # Wrap positions
        for i in 1:size(State.Coords, 1)
            for dim in 1:3
                if State.Coords[i, dim] > bl
                    State.Coords[i, dim] -= 2*bl
                elseif State.Coords[i, dim] < -bl
                   State.Coords[i, dim] += 2*bl
                end
            end
        end
    end
    
    return nothing
end

"""
    run_simulation!(system::System, params::SimulationParams, 
                   force_func, potential_func; 
                   boundary_type=:reflective, verbose=true)

Run a full MD simulation.

Returns: (trajectory, energies, temperatures, potential_energies)
"""
function run_simulation!(Initial_state::SyState,system::System, params::SimulationParams;
                        boundary_type=:reflective, verbose=true)
    
    n_particles = system.N_atoms#size(Initial_state.Coords, 1)
    n_save = params.n_steps ÷ params.output_freq
    
    # Preallocate storage
    trajectory_xyz = zeros(n_save, n_particles, 3)
    energies = zeros(n_save)
    temperatures = zeros(n_save)
    potential_energies = zeros(n_save)
    Simulation_state = deepcopy(Initial_state)

    if verbose
        println("Running MD simulation")
        println("  Particles: $n_particles")
        println("  Steps: $(params.n_steps)")
        println("  dt: $(params.dt)")
        println("  Output frequency: $(params.output_freq)")
    end
    
    # Initial force calculation
    calculate_forces!(Simulation_state, system)
    
    save_idx = 1
    
    # Main MD loop
    for step in 1:params.n_steps
        # Integration step
        velocity_verlet!(Simulation_state, system, params.dt)
        
        # Apply boundary conditions
        apply_boundary_conditions!(Simulation_state, system, boundary_type)
        
        # Save data at specified intervals
        if step % params.output_freq == 0
            trajectory_xyz[save_idx, :, :] = Simulation_state.Coords
            energies[save_idx] = Simulation_state.TotE
            temperatures[save_idx] = temperature(Simulation_state)
            potential_energies[save_idx] = Simulation_state.Ene
            
            if verbose && (save_idx % 10 == 0 || save_idx == 1)
                @printf("Step %6d: E = %12.6f, T = %12.6f, PE = %12.6f\n",
                       step, energies[save_idx], temperatures[save_idx], 
                       potential_energies[save_idx])
            end
            
            save_idx += 1
        end
    end
    trajectory = Trajectory(trajectory_xyz,energies,temperatures,potential_energies)

    return trajectory, Simulation_state 
end

"""
    run_damped_optimization!(system::System, params::SimulationParams,
                            force_func, potential_func; 
                            verbose=true, tol=1e-12)

Run damped dynamics for geometry optimization.
Velocity is set to zero when power becomes negative.
"""
function run_damped_optimization!(Initial_state::SyState,system::System, params::SimulationParams;
                                 verbose=true, tol=1e-12)
    
    n_particles = system.N_atoms#size(Initial_state.Coords, 1)
    n_save = params.n_steps ÷ params.output_freq
    
    # Preallocate storage
    trajectory_xyz = zeros(n_save, n_particles, 3)
    energies = zeros(n_save)
    temperatures = zeros(n_save)
    potential_energies = zeros(n_save)
    Simulation_state = deepcopy(Initial_state)

    if verbose
        println("Running MD simulation")
        println("  Particles: $n_particles")
        println("  Steps: $(params.n_steps)")
        println("  dt: $(params.dt)")
        println("  Output frequency: $(params.output_freq)")
    end
    
    # Initial force calculation
    calculate_forces!(Simulation_state, system)
    
    save_idx = 1
    
    # Main MD loop
    for step in 1:params.n_steps
        # Integration step
        velocity_verlet!(Simulation_state, system, params.dt)
        ## Check power (F·v)
        power = sum(Simulation_state.force .* Simulation_state.Vel)
        
        # Damping: zero velocity if power is negative
        if power <= 0.0
            fill!(Simulation_state.Vel, 0.0)
        end
        
        # Save data at specified intervals
        if step % params.output_freq == 0
            trajectory_xyz[save_idx, :, :] = Simulation_state.Coords
            energies[save_idx] = Simulation_state.TotE
            temperatures[save_idx] = temperature(Simulation_state)
            potential_energies[save_idx] = Simulation_state.Ene
            
            if verbose && (save_idx % 10 == 0 || save_idx == 1)
                @printf("Step %6d: E = %12.6f, T = %12.6f, PE = %12.6f\n",
                       step, energies[save_idx], temperatures[save_idx], 
                       potential_energies[save_idx])
            end
            
            save_idx += 1
        end
    end
    trajectory = Trajectory(trajectory_xyz,energies,temperatures,potential_energies)

    return trajectory, Simulation_state 
end

function run_damped_optimization_FIRE!(system::System, params::SimulationParams, Fire_params::FIRE_params; 
                                 verbose=true, tol=1e-12)
    
    n_particles = size(system.positions, 1)
    max_saves = params.n_steps ÷ params.output_freq
    
    # Storage (will trim at end)
    trajectory = zeros(max_saves, n_particles, 3)
    energies = zeros(max_saves)
    temperatures = zeros(max_saves)
    potential_energies = zeros(max_saves)
    
    if verbose
        println("Running damped optimization")
        println("  Particles: $n_particles")
        println("  Max steps: $(params.n_steps)")
    end
    
    if params.dt != Fire_params.dt
        print("dt from FIRE params different from that in simulation params\n using dt dt from FIRE params")
    end
    # Initial force calculation
    calculate_forces!(system)
    
    save_idx = 1
    prev_pe = system.potential_energy
    
    for step in 1:params.n_steps
        # Save old forces
        old_forces = copy(system.forces)
        old_velocity = copy(system.velocities)
        ## Update positions
        velocity_verlet!(system,Fire_params.dt)
        ## Check power (F·v)
        power = sum(system.forces .* system.velocities)
        ## FIRE steps
        # Damping: zero velocity if power is negative
        if power <= 0.0
            # system updates
            fill!(system.velocities, 0.0)
            system.positions -= Fire_params.dt*old_velocity + 0.25*Fire_params.dt^2*(old_forces+system.forces)
            # parameter updates
            if Fire_params.ddt >= Fire_params.dt
                Fire_params.ddt = Fire_params.dt/10.0
            end
            Fire_params.a = Fire_params.a0
            Fire_params.dt -= Fire_params.ddt
        else
            system.velocities = (1.0 - Fire_params.a)*system.velocities + Fire_params.a*norm(system.velocities)*system.forces/norm(system.forces)
            if Fire_params.da >= Fire_params.a
                Fire_params.da = Fire_params.a/10.0
            end
            Fire_params.a -= Fire_params.da
            Fire_params.dt += Fire_params.ddt
            if Fire_params.a < Fire_params.a_min
                Fire_params.a = Fire_params.a_min
            end
            if Fire_params.dt > Fire_params.dt_max
                Fire_params.dt = Fire_params.dt_max
            end
        end
        
        # Save data
        if step % params.output_freq == 0
            trajectory[save_idx, :, :] = system.positions
            energies[save_idx] = total_energy(system)
            temperatures[save_idx] = temperature(system)
            potential_energies[save_idx] = system.potential_energy
            
            if verbose && (save_idx % 10 == 0 || save_idx == 1)
                @printf("Step %6d: PE = %12.6f, |F| = %12.6f, |v| = %12.6f\n",
                       step, potential_energies[save_idx], 
                       norm(system.forces), norm(system.velocities))
            end
            
            # Check convergence
            pe_change = abs(system.potential_energy - prev_pe)
            force_norm = norm(system.forces)
            vel_norm = norm(system.velocities)
            
            if pe_change < tol && force_norm < tol && vel_norm < tol
                if verbose
                    println("Convergence reached after $step steps")
                end
                save_idx += 1
                break
            end
            
            prev_pe = system.potential_energy
            save_idx += 1
        end
    end
    
    # Trim arrays
    actual_saves = save_idx - 1
    trajectory = Trajectory(trajectory[1:actual_saves, :, :],energies[1:actual_saves],temperatures[1:actual_saves],potential_energies[1:actual_saves])
    return trajectory
            
end


