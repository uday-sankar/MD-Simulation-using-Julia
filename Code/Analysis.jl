
# =============================================================================
# File: src/analysis.jl
# =============================================================================

"""
    kinetic_energy(system::System)

Calculate total kinetic energy: KE = Σ 0.5 * m * v²
"""
function kinetic_energy(system::System)
    ke = 0.0
    for i in 1:size(system.positions, 1)
        v_sq = sum(system.velocities[i, :].^2)
        ke += 0.5 * system.masses[i] * v_sq
    end
    return ke
end

"""
    potential_energy(system::System)

Return stored potential energy.
"""
function potential_energy(system::System)
    return system.potential_energy
end

"""
    total_energy(system::System)

Calculate total energy: E = KE + PE
"""
function total_energy(system::System)
    return kinetic_energy(system) + system.potential_energy
end

"""
    temperature(system::System; kB=1.0)

Calculate instantaneous temperature from kinetic energy.
T = 2*KE / (N_dof * kB), where N_dof = 3N - 3 (removing COM motion)
"""
function temperature(system::System; kB=1.0)
    n_particles = size(system.positions, 1)
    ke = kinetic_energy(system)
    n_dof = 3 * n_particles - 3  # Remove COM translation
    return 2.0 * ke / (n_dof * kB)
end

"""
    center_of_mass(system::System)

Calculate center of mass: R_COM = Σ m_i * r_i / Σ m_i
"""
function center_of_mass(system::System)
    total_mass = sum(system.masses)
    com = zeros(3)
    
    for i in 1:size(system.positions, 1)
        com .+= system.masses[i] * system.positions[i, :]
    end
    
    return com ./ total_mass
end

"""
    radius_of_gyration(system::System)

Calculate radius of gyration: R_g² = Σ m_i |r_i - R_COM|² / Σ m_i
"""
function radius_of_gyration(system::System)
    com = center_of_mass(system)
    total_mass = sum(system.masses)
    rg_sq = 0.0
    
    for i in 1:size(system.positions, 1)
        r_rel = system.positions[i, :] - com
        rg_sq += system.masses[i] * sum(r_rel.^2)
    end
    
    return sqrt(rg_sq / total_mass)
end

"""
    mean_interparticle_distance(system::System)

Calculate mean distance between all pairs of particles.
"""
function mean_interparticle_distance(system::System)
    n = size(system.positions, 1)
    distances = Float64[]
    
    for i in 1:n-1
        for j in i+1:n
            r = norm(system.positions[i, :] - system.positions[j, :])
            push!(distances, r)
        end
    end
    
    return sum(distances) / length(distances)
end