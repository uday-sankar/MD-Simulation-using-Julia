# =============================================================================
# File: src/forces.jl
# =============================================================================

"""
calculate_forces!(system::System, force_func, potential_func)

Calculate all pairwise forces in the system and update system.forces.
Also updates system.potential_energy.

# Arguments
- `system`: System object to update
- `force_func`: Function f(r_vec) returning force vector
- `potential_func`: Function V(r_vec) returning potential energy
"""
function calculate_forces!(system::System, force_func, potential_func)
    n = size(system.positions, 1)
    
    # Zero out forces
    fill!(system.forces, 0.0)
    system.potential_energy = 0.0
    
    # Pairwise interactions
    for i in 1:n-1
        for j in i+1:n
            r_vec = system.positions[i, :] - system.positions[j, :]
            
            # Calculate force and potential
            force = force_func(r_vec)
            potential = potential_func(r_vec)
            
            # Newton's third law
            system.forces[i, :] .+= force
            system.forces[j, :] .-= force
            
            # Accumulate potential energy
            system.potential_energy += potential
        end
    end
    
    return nothing
end

"""
    harmonic_potential(r_vec; k=1.0, r_eq=2.0)

Harmonic potential: V = 0.5 * k * |r - r_eq|²
"""
function harmonic_potential(r_vec; k=1.0, r_eq=2.0)
    r = norm(r_vec)
    return 0.5 * k * (r - r_eq)^2
end

"""
    harmonic_force(r_vec; k=1.0, r_eq=2.0)

Force from harmonic potential: F = -k * (r - r_eq) * r_hat
"""
function harmonic_force(r_vec; k=1.0, r_eq=2.0)
    r = norm(r_vec)
    if r < 1e-10
        return zeros(3)
    end
    r_hat = r_vec / r
    return -k * (r - r_eq) * r_hat
end

"""
    lennard_jones_potential(r_vec; epsilon=1.0, sigma=1.0, cutoff=Inf)

Lennard-Jones 12-6 potential: V = 4ε[(σ/r)¹² - (σ/r)⁶]
"""
function lennard_jones_potential(r_vec; epsilon=1.0, sigma=1.0, cutoff=Inf)
    r = norm(r_vec)
    
    if r > cutoff
        return 0.0
    end
    
    sr6 = (sigma / r)^6
    sr12 = sr6^2
    
    return 4.0 * epsilon * (sr12 - sr6)
end

"""
    lennard_jones_force(r_vec; epsilon=1.0, sigma=1.0, cutoff=Inf)

Force from Lennard-Jones potential.
"""
function lennard_jones_force(r_vec; epsilon=1.0, sigma=1.0, cutoff=Inf)
    r = norm(r_vec)
    
    if r > cutoff || r < 1e-10
        return zeros(3)
    end
    
    sr6 = (sigma / r)^6
    sr12 = sr6^2
    
    # F = 24ε/r * [2(σ/r)¹² - (σ/r)⁶] * r_hat
    force_mag = 24.0 * epsilon / r^2 * (2.0 * sr12 - sr6)
    
    return force_mag * r_vec
end


