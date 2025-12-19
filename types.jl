# =============================================================================
# File: src/types.jl
# =============================================================================

"""
    System

Mutable struct containing all system state information.

# Fields
- `positions::Matrix{Float64}`: Particle positions (N×3)
- `velocities::Matrix{Float64}`: Particle velocities (N×3)
- `forces::Matrix{Float64}`: Forces on particles (N×3)
- `masses::Vector{Float64}`: Particle masses (N,)
- `box_size::Float64`: Boundary box size
- `time::Float64`: Current simulation time
- `potential_energy::Float64`: Current potential energy
"""
mutable struct System
    positions::Matrix{Float64}      # N × 3
    velocities::Matrix{Float64}     # N × 3
    forces::Matrix{Float64}         # N × 3
    masses::Vector{Float64}         # N
    box_size::Float64               # Scalar
    time::Float64                   
    potential_energy::Float64       # Scalar
    Force_func::Any            # Function: vetor (Nx3) -> vector (Nx3)
    Potential_func::Any       # Function: vector (Nx3) -> Float
end

"""
    System(positions, velocities, masses, box_size)

Construct a System with initial positions, velocities, masses, and box size.
Forces are initialized to zero.
"""
function System(positions::Matrix{Float64}, 
                velocities::Matrix{Float64},
                masses::Vector{Float64},
                box_size::Float64,F::Any,V::Any)
    n_particles = size(positions, 1)
    forces = zeros(n_particles, 3)
    system = System(positions, velocities, forces, masses, box_size, 0.0, 0.0,F,V)
    calculate_forces!(system)
    return system
end

"""
    SimulationParams

Immutable struct for simulation parameters.
"""
struct SimulationParams
    n_steps::Int 
    dt::Float64
    output_freq::Int 
    boundary_size::Float64
    function SimulationParams(;n_steps = 100, dt=0.1, output_freq = 5, boundary_size = 10)
        new(n_steps,dt,output_freq,boundary_size)
    end
end

mutable struct Trajectory
    Trajectory_coords::Array{Float64,3}   # T x N x 3
    Tot_Energy::Array{Float64,1}   # T x 1
    Temperature::Array{Float64,1}  # T x 1
    PE::Array{Float64,1}  # T x 1
end

