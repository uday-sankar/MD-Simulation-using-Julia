# =============================================================================
# File: src/types.jl
# =============================================================================

"""
    SyState

    Mutable struct which contain all the information about the current state of the system
    Fields:
    - Coords : Positions
    - Vel : Velocity
    - force : Forces
    - M : mass
    - Ene : Potential Energy
    - TotE : Total Energy of the system if needed 
    - t : time 
"""
mutable struct SyState 
    Atoms::Vector{String}
    Coords::Matrix{Float64}     # N x 3
    Vel::Matrix{Float64}        # N × 3
    force::Matrix{Float64}      # N × 3
    M::Matrix{Float64}          # N x 1
    Ene::Float64                # Scalar
    TotE::Float64               # Scalar
    t::Float64                  # Scalar 
end

"""
SysState 
"""
#function SyState(Atoms::Vector{String},Coords::Matrix{Float64},Vel::Matrix{Float64},force::Matrix{Float64},M::Vector{Float64}, Ene::Float64, TotE::Float64, t::Float64)
#    return SyState(Atoms)
#end

"""
    System

Mutable struct containing all system state information.

# Fields
- Atoms : A vector containing the symbol of atoms involved.
- `box_size::Float64`: Boundary box size
- Force_func : The function describing the force
- Potenial_func : The potential energy surface
- Inter_atomic_tag : A tag that determines whether the force and potential functions are inter atomic or global 
"""
mutable struct System
#   positions::Matrix{Float64}      # N × 3
#   velocities::Matrix{Float64}     # N × 3
#   forces::Matrix{Float64}         # N × 3
#   masses::Vector{Float64}         # N
#   time::Float64                   # Scalar
#   potential_energy::Float64       # Scalar
#   State_init::SyState             # Intial state of the system
    Atoms::Vector{String}           # N x 1
    N_atoms::Int64                  # Scalar (N)
    Force_func::Function                 # Function: vetor (Nx3) -> vector (Nx3)
    Potential_func::Function            # Function: vector (Nx3) -> Float
    box_size::Float64               # Scalar
    Inter_atomic_tag::Bool          # true/false

end

"""
    System( Atoms, box_size, force_func, P_func, Inter_atomic_tag)

Construct a System with initial positions, velocities, masses, and box size.
Forces are initialized to zero.
"""
function System(Atoms::Vector{String}, F::Function, V::Function, box_size::Float64, inter_atomic::Bool=true)
    N_atoms = size(Atoms, 1)
    system = System(Atoms, N_atoms, F, V, box_size, inter_atomic)
    #calculate_forces!(system)
    return system
end

# """
#     System(positions, velocities, masses, box_size)
# 
# Construct a System with initial positions, velocities, masses, and box size.
# Forces are initialized to zero.
# """
# function System(positions::Matrix{Float64}, 
#                 velocities::Matrix{Float64},
#                 masses::Vector{Float64},
#                 box_size::Float64,F::Any,V::Any,inter_atomic::Bool=true)
#     n_particles = size(positions, 1)
#     forces = zeros(n_particles, 3)
#     system = System(positions, velocities, forces, masses, box_size, 0.0, 0.0, F, V, inter_atomic)
#     calculate_forces!(system)
#     return system
# end

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
    Atoms::Vector{String}
end

mutable struct FIRE_params
    a0::Float64
    a::Float64
    da::Float64
    a_min::Float64
    dt::Float64
    ddt::Float64
    dt_max::Float64
    function FIRE_params(;a0=1.0,a=1.0,da=0.01,a_min=0.1,dt=0.1,ddt=0.01,dt_max=0.1)
        new(a0,a,da,a_min,dt,ddt,dt_max)
    end
end
