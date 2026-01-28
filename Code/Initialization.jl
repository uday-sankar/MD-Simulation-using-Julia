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
                        masses=nothing,Atoms=nothing, output_dir="./", animate=false)
    
    println("Grid Stabilization Started")
    
    # Initialize
    positions, velocities = grid_initialize(n_grid, spacing)
    n_particles = size(positions, 1)
    
    if masses === nothing
        masses = ones(n_particles)
    end
    if Atoms === nothing
        Atoms = ["Ar" for i in 1:size(positions,1)]
    end
    # Create system
    box_size = max(n_grid[1]*spacing[1], n_grid[2]*spacing[2], n_grid[3]*spacing[3])
    # System defnetion  
    system = MD_System(Atoms, force_func,potential_func,box_size)
    # State defenitions
    mass_mat = hcat(masses,masses,masses)
    State_init = SyState(Atoms,positions,velocities,velocities*0.0,mass_mat,0.0,0.0,0.0)
    calculate_forces!(State_init)
    dt = 0.05
    println("  Number of particles: $n_particles")
    
    # Run optimization
    trajectory, Final_State = run_damped_optimization!(State_init,system, params)
    PE = trajectory.PE
    println("  PE change: $(PE[end] - PE[1])")
    
    # Final geometry (centered at origin)
    #com = center_of_mass(system)
    #final_positions = system.positions .- com'
    #trajectory = Trajectory(trajectory_xyz,E,T,PE)
    #
    return trajectory, system, Final_State
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
    system = MD_System(positions, velocities, masses, box_size,force_func,potential_func)
    
    println("  Number of particles: $n_particles")
    
    # Run optimization
    trajectory= run_damped_optimization!(system, params)
    
    println("  PE change: $(PE[end] - PE[1])")
    
    # Final geometry (centered at origin)
    com = center_of_mass(system)
    final_positions = system.positions .- com'
    
    return final_positions, trajectory
end

"""
    Initialation with zero net momenta
"""
function Recenter_and_vel(MD_State, ET)
    # MD_State: The current state to be modified
    # ET: Total kinetic energy to be given
    # This function recenters the given final MD state.
    # The cenetr of mass (CoM) is kept at the origin.
    # Velocities are given such that the CoM remain fixed.
    New_State = deepcopy(MD_State)
    masses = New_State.M[:,1]
    N = length(New_State.Atoms)# number of particles
    com = center_of_mass(New_State)
    CoM = vcat(com,com,com,com)
    # Receter geometry
    New_State.Coords = New_State.Coords - CoM
    ## Getting velocities for the N-1 particles
    vel_3 = (rand(N-1,3) .- 0.5)*2.0
    m_vel_3 = masses[1:N-1]' * vel_3
    vel_4 = -m_vel_3/masses[N]
    full_vel = vcat(vel_3,vel_4)
    KE = 0.5*sum([masses[i]*(full_vel[i,:]'*full_vel[i,:]) for i in 1:N])
    # scaling velocity to math total energy
    full_vel = full_vel*(ET/KE)^0.5# scale each velocity by ET/KE
    New_State.Vel = deepcopy(full_vel)
    #Moment of inertial test
    Ang_M = 0.0
    for i in 1:N
        Ang_M += masses[i]*norm(cross(New_State.Coords[i,:],full_vel[i,:]))
    end
    KE = 0.5*sum([masses[i]*(full_vel[i,:]'*full_vel[i,:]) for i in 1:N])
    print("Recentering Done:\n\t Total KE:$KE \t Angular Momentum: $Ang_M")
    return New_State
end