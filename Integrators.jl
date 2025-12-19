
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
function velocity_verlet!(system::System, dt::Float64, force_func, potential_func)
    n = size(system.positions, 1)
    
    # Save old forces
    old_forces = copy(system.forces)
    
    # Update positions
    for i in 1:n
        system.positions[i, :] .+= system.velocities[i, :] * dt .+ 
                                    0.5 * old_forces[i, :] / system.masses[i] * dt^2
    end
    
    # Calculate new forces
    calculate_forces!(system, force_func, potential_func)
    
    # Update velocities using average force
    for i in 1:n
        avg_force = 0.5 * (old_forces[i, :] + system.forces[i, :])
        system.velocities[i, :] .+= avg_force / system.masses[i] * dt
    end
    
    # Update time
    system.time += dt
    
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
function apply_boundary_conditions!(system::System, boundary_type=:reflective)
    if boundary_type == :none
        return nothing
    end
    
    bl = system.box_size
    
    if boundary_type == :reflective
        # Check each particle
        for i in 1:size(system.positions, 1)
            for dim in 1:3
                # Check if particle is outside boundary
                if abs(system.positions[i, dim]) >= bl
                    # Reverse velocity component
                    system.velocities[i, dim] *= -1.0
                    # Clamp position to boundary
                    system.positions[i, dim] = sign(system.positions[i, dim]) * (bl - 1e-6)
                end
            end
        end
    elseif boundary_type == :periodic
        # Wrap positions
        for i in 1:size(system.positions, 1)
            for dim in 1:3
                if system.positions[i, dim] > bl
                    system.positions[i, dim] -= 2*bl
                elseif system.positions[i, dim] < -bl
                    system.positions[i, dim] += 2*bl
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
function run_simulation!(system::System, params::SimulationParams,
                        force_func, potential_func;
                        boundary_type=:reflective, verbose=true)
    
    n_particles = size(system.positions, 1)
    n_save = params.n_steps ÷ params.output_freq
    
    # Preallocate storage
    trajectory = zeros(n_save, n_particles, 3)
    energies = zeros(n_save)
    temperatures = zeros(n_save)
    potential_energies = zeros(n_save)
    
    if verbose
        println("Running MD simulation")
        println("  Particles: $n_particles")
        println("  Steps: $(params.n_steps)")
        println("  dt: $(params.dt)")
        println("  Output frequency: $(params.output_freq)")
    end
    
    # Initial force calculation
    calculate_forces!(system, force_func, potential_func)
    
    save_idx = 1
    
    # Main MD loop
    for step in 1:params.n_steps
        # Integration step
        velocity_verlet!(system, params.dt, force_func, potential_func)
        
        # Apply boundary conditions
        apply_boundary_conditions!(system, boundary_type)
        
        # Save data at specified intervals
        if step % params.output_freq == 0
            trajectory[save_idx, :, :] = system.positions
            energies[save_idx] = total_energy(system)
            temperatures[save_idx] = temperature(system)
            potential_energies[save_idx] = system.potential_energy
            
            if verbose && (save_idx % 10 == 0 || save_idx == 1)
                @printf("Step %6d: E = %12.6f, T = %12.6f, PE = %12.6f\n",
                       step, energies[save_idx], temperatures[save_idx], 
                       potential_energies[save_idx])
            end
            
            save_idx += 1
        end
    end
    
    return trajectory, energies, temperatures, potential_energies
end

"""
    run_damped_optimization!(system::System, params::SimulationParams,
                            force_func, potential_func; 
                            verbose=true, tol=1e-12)

Run damped dynamics for geometry optimization.
Velocity is set to zero when power becomes negative.
"""
function run_damped_optimization!(system::System, params::SimulationParams,
                                 force_func, potential_func;
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
    
    # Initial force calculation
    calculate_forces!(system, force_func, potential_func)
    
    save_idx = 1
    prev_pe = system.potential_energy
    
    for step in 1:params.n_steps
        # Save old forces
        old_forces = copy(system.forces)
        
        # Update positions
        for i in 1:n_particles
            system.positions[i, :] .+= system.velocities[i, :] * params.dt .+ 
                                        0.5 * old_forces[i, :] / system.masses[i] * params.dt^2
        end
        
        # Calculate new forces
        calculate_forces!(system, force_func, potential_func)
        
        # Update velocities
        for i in 1:n_particles
            avg_force = 0.5 * (old_forces[i, :] + system.forces[i, :])
            system.velocities[i, :] .+= avg_force / system.masses[i] * params.dt
        end
        
        # Check power (F·v)
        power = sum(system.forces .* system.velocities)
        
        # Damping: zero velocity if power is negative
        if power <= 0.0
            fill!(system.velocities, 0.0)
        end
        
        system.time += params.dt
        
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
    return (trajectory[1:actual_saves, :, :], 
            energies[1:actual_saves],
            temperatures[1:actual_saves],
            potential_energies[1:actual_saves])
end


