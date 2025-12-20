using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra

include("Code/MD_Base.jl")
using .MD_Base
##
# Your data vectors (replace with yours)
# time = 0:0.001:100.0  # ps
# T = randn(length(time)) .+ 300.0  # K, noisy example
# E_total = cumsum(randn(length(time))) .- 1000.0  # kJ/mol
# PE = cumsum(randn(length(time))) .- 5000.0       # kJ/mol
using Plots
##
# =============================================================================
fold = "Version 4"
# Define force and potential functions
function harmonic_force_example(r_vec)
    return MD_Base.harmonic_force(r_vec; k=1.0, r_eq=4.0)
end

function harmonic_potential_example(r_vec)
    return MD_Base.harmonic_potential(r_vec; k=1.0, r_eq=4.0)
end

# Initialize system on a grid
positions, velocities = grid_initialize([3, 4, 3], [2.0, 2.0, 2.0])
n_particles = size(positions, 1)
masses = ones(n_particles)

# Create system
system = System(positions, velocities, masses, 15.0,harmonic_force_example,harmonic_potential_example)

# Set simulation parameters
dt = 0.05
params = SimulationParams(n_steps=1000,dt=dt,output_freq=1, boundary_size=15)
##

# Run simulation
trajectory= run_simulation!(
    system, params,
    boundary_type = :reflective,
    verbose = true
)
# Save results
write_trajectory(trajectory, "$fold/harmonic_trajectory.xyz")
write_enhanced_trajectory(trajectory, "$fold/harmonic_detailed.xyz",
                         dt=params.dt)

# Visualize (if Plots is available)
try
    animate_trajectory(trajectory, filename="$fold/harmonic.gif", 
                      box_size=params.boundary_size, frames=50)
    plot_snapshots(trajectory, [1, :mid, :end], 
                  filename="$fold/harmonic_snapshots.png")
catch e
    println("Visualization skipped (Plots.jl may not be available)")
end

energies = trajectory.Tot_Energy
T = trajectory.Temperature
PE = trajectory.PE
##
t_end = system.time
Plot_data(trajectory,0:dt:t_end; smooth_flag=false,fold=fold,window_smooth=40)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
# =======================================================================================
##
# Grid stabilization
using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra

include("Code/MD_Base.jl")
using .MD_Base
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
fold = "Version 4/Vel-Ver"
println("\n" * "="^70)
println("Example: Lennard-Jones Cluster Optimization")
println("="^70)

# Define LJ force and potential
function lj_force_example(r_vec)
    return MD_Base.lennard_jones_force(r_vec; epsilon=1.0, sigma=3.0, cutoff=8.0)
end

function lj_potential_example(r_vec)
    return MD_Base.lennard_jones_potential(r_vec; epsilon=1.0, sigma=3.0, cutoff=8.0)
end

# Grid stabilization
opt_params = SimulationParams(
    n_steps = 100,
    dt = 0.01,
    output_freq = 1,
    boundary_size = 20.0
)

pos, vel = grid_initialize([3,3,3], [4.0,4.0,4.0])
mass = ones(size(pos,1))
system = System(pos, vel, mass, 10.0, lj_force_example, lj_potential_example)

Fire_params = FIRE_params(dt=0.01)

trajectory_opt = run_damped_optimization_FIRE!(system, opt_params, Fire_params)

trajectory_opt = run_damped_optimization!(system, opt_params)

final_geom = system.positions
# Save optimized geometry
write_xyz(final_geom, "$fold/lj_optimized_vel-ver.xyz")

t_end = system.time
Time = 1:1:size(trajectory_opt.PE,1)
Plot_data(trajectory_opt,Time; smooth_flag=false,fold=fold,window_smooth=40,plt_tot_E=false)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(trajectory_opt.Tot_Energy[end] - trajectory_opt.Tot_Energy[1]))")
println("  Mean temperature: $(sum(trajectory_opt.Temperature)/length(trajectory_opt.Temperature))")

write_trajectory(trajectory_opt,"$fold/trajectory_opt.xyz")
##
# ==============================================================================================================================================
##
using Molly

# Setup Molly system for force calculations only
function setup_molly_system(positions::Matrix{Float64}, atom_types::Vector{String})
    n_atoms = size(positions, 1)
    
    # Create Molly atoms with LJ parameters
    # You can load these from forcefield files
    atoms = [Atom(
        index=i, 
        charge=0.0,
        mass=get_mass(atom_types[i]),
        σ=get_sigma(atom_types[i]),  # LJ sigma
        ϵ=get_epsilon(atom_types[i])  # LJ epsilon
    ) for i in 1:n_atoms]
    
    # Convert positions to Molly format (with units)
    coords = [SVector(positions[i,:]...)u"Å" for i in 1:n_atoms]
    
    # Define interactions
    lj_inter = LennardJones(
        use_neighbors=true,
        cutoff=DistanceCutoff(10.0u"Å")
    )
    
    # Create system
    sys = System(
        atoms=atoms,
        coords=coords,
        boundary=CubicBoundary(20.0u"Å"),
        pairwise_inters=(lj_inter,)
    )
    
    return sys
end

# Wrapper for your MD code
function molly_force_wrapper(positions::Matrix{Float64}, molly_sys)
    n_atoms = size(positions, 1)
    
    # Update positions in Molly system
    new_coords = [SVector(positions[i,:]...)u"Å" for i in 1:n_atoms]
    updated_sys = Molly.System(
        molly_sys;
        coords=new_coords
    )
    
    # Calculate forces and energy
    forces_molly = Molly.forces(updated_sys)
    energy = Molly.potential_energy(updated_sys)
    
    # Convert back to your format (N×3 matrix, no units)
    forces = zeros(n_atoms, 3)
    for i in 1:n_atoms
        forces[i, :] = ustrip.(u"kJ * mol^-1 * Å^-1", forces_molly[i])
    end
    energy_val = ustrip(u"kJ/mol", energy)
    
    return energy_val, forces
end

# Use in your existing code
molly_sys = setup_molly_system(initial_positions, atom_types)

# Create closure for your integrator
force_func = (pos) -> molly_force_wrapper(pos, molly_sys)[2]
potential_func = (pos) -> molly_force_wrapper(pos, molly_sys)[1]

# Use with your existing MDSimulation code
trajectory, E, T, PE = run_simulation!(
    system, params,
    force_func,
    potential_func
)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
fold = "Version 4"
println("\n" * "="^70)
println("Example: Lennard-Jones Cluster Optimization")
println("="^70)

# Define LJ force and potential
function lj_force_example(r_vec)
    return MD_Base.lennard_jones_force(r_vec; epsilon=1.0, sigma=3.0, cutoff=8.0)
end

function lj_potential_example(r_vec)
    return MD_Base.lennard_jones_potential(r_vec; epsilon=1.0, sigma=3.0, cutoff=8.0)
end

# Grid stabilization
opt_params = SimulationParams(
    n_steps = 5000,
    dt = 0.01,
    output_freq = 1,
    boundary_size = 20.0
)

opt_trajectory, system = stabilize_grid!(
    [4, 4, 4],           # Grid dimensions
    [5.5, 5.5, 5.5],     # Spacing
    opt_params,
    lj_force_example,
    lj_potential_example;
    output_dir = "./"
)
final_geom = system.positions
# Save optimized geometry
write_xyz(final_geom, "$fold/lj_optimized.xyz")

t_end = system.time
Time = 0:0.01:t_end#size(opt_trajectory.PE,1)
Plot_data(opt_trajectory,Time; smooth_flag=false,fold=fold,window_smooth=40,plt_tot_E=false)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
