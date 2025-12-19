# =============================================================================
# File: src/initialization.jl
# =============================================================================

"""
    grid_initialize(n_grid, spacing)

Initialize particles on a cubic grid.

# Arguments
- `n_grid`: Vector [Nx, Ny, Nz] of grid points in each dimension
- `spacing`: Vector [dx, dy, dz] of spacing between grid points

# Returns
- `positions`: N×3 matrix of positions
- `velocities`: N×3 matrix of velocities (initialized to zero)
"""
function grid_initialize(n_grid::Vector{Int}, spacing::Vector{Float64})
    Nx, Ny, Nz = n_grid
    dx, dy, dz = spacing
    
    # Create grid ranges centered at origin
    x_range = -(Nx-1)*dx/2 : dx : (Nx-1)*dx/2
    y_range = -(Ny-1)*dy/2 : dy : (Ny-1)*dy/2
    z_range = -(Nz-1)*dz/2 : dz : (Nz-1)*dz/2
    
    n_particles = Nx * Ny * Nz
    positions = zeros(n_particles, 3)
    
    idx = 1
    for i in 1:Nx, j in 1:Ny, k in 1:Nz
        positions[idx, :] = [x_range[i], y_range[j], z_range[k]]
        idx += 1
    end
    
    velocities = zeros(n_particles, 3)
    
    return positions, velocities
end

"""
    random_initialize(n_particles, box_size)

Initialize particles randomly in a box.

# Arguments
- `n_particles`: Number of particles
- `box_size`: Size of cubic box (particles placed in [-L, L]³)

# Returns
- `positions`: N×3 matrix of random positions
- `velocities`: N×3 matrix of velocities (initialized to zero)
"""
function random_initialize(n_particles::Int, box_size::Float64)
    positions = (2 * box_size) .* rand(n_particles, 3) .- box_size
    velocities = zeros(n_particles, 3)
    return positions, velocities
end

"""
    stabilize_grid!(n_grid, spacing, params, force_func, potential_func;
                   masses=nothing, output_dir="./", animate=false)

Initialize grid and run damped optimization to find stable configuration.
"""
function stabilize_grid!(n_grid, spacing, params, force_func, potential_func;
                        masses=nothing, output_dir="./", animate=false)
    
    println("Grid Stabilization Started")
    
    # Initialize
    positions, velocities = grid_initialize(n_grid, spacing)
    n_particles = size(positions, 1)
    
    if masses === nothing
        masses = ones(n_particles)
    end
    
    # Create system
    box_size = max(n_grid[1]*spacing[1], n_grid[2]*spacing[2], n_grid[3]*spacing[3])
    system = System(positions, velocities, masses, box_size, force_func, potential_func)
    
    println("  Number of particles: $n_particles")
    
    # Run optimization
    trajectory_xyz, E, T, PE = run_damped_optimization!(system, params, 
                                                   force_func, potential_func)
    
    println("  PE change: $(PE[end] - PE[1])")
    
    # Final geometry (centered at origin)
    com = center_of_mass(system)
    #final_positions = system.positions .- com'
    trajectory = Trajectory(trajectory_xyz,E,T,PE)
    #
    return trajectory, system
end

"""
    stabilize_random!(n_particles, box_size, params, force_func, potential_func;
                     masses=nothing, output_dir="./")

Initialize random configuration and optimize.
"""
function stabilize_random!(n_particles, box_size, params, force_func, potential_func;
                          masses=nothing, output_dir="./")
    
    println("Random Stabilization Started")
    
    # Initialize
    positions, velocities = random_initialize(n_particles, box_size)
    
    if masses === nothing
        masses = ones(n_particles)
    end
    
    # Create system
    system = System(positions, velocities, masses, box_size,force_func,potential_func)
    
    println("  Number of particles: $n_particles")
    
    # Run optimization
    trajectory_xyz, E, T, PE = run_damped_optimization!(system, params,
                                                   force_func, potential_func)
    
    println("  PE change: $(PE[end] - PE[1])")
    
    # Final geometry (centered at origin)
    com = center_of_mass(system)
    final_positions = system.positions .- com'
    
    trajectory = Trajectory(trajectory_xyz,E,T,PE)

    return final_positions, trajectory
end
