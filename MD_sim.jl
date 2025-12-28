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
fold = "Testing/Version 5:  Molecules"
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
masses = 40.0*ones(n_particles)
boxsize = 10.0
##
Atoms = ["Ar" for i in 1:size(positions,1)]
# Create system
system = System(Atoms, harmonic_force_example,harmonic_potential_example,boxsize)
##
mass_mat = hcat(masses,masses,masses)
State_init = MD_Base.SyState(Atoms,positions,velocities,velocities*0.0,mass_mat,0.0,0.0,0.0)
# Set simulation parameters
dt = 0.05
params = SimulationParams(n_steps=1000,dt=dt,output_freq=1, boundary_size=15)
##
calculate_forces!(State_init,system)

# Run simulation
trajectory, Final_State = run_simulation!(State_init,
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
    animate_trajectory(trajectory, filename="$fold/harmonic.mp4", 
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
t_end = Final_State.t
Plot_data(trajectory,0:dt:t_end; smooth_flag=false,fold=fold,window_smooth=40)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
# =======================================================================================
##
# Grid stabilization
dt = 0.001
opt_params = SimulationParams(
    n_steps = 500,
    dt = dt,
    output_freq = 1,
    boundary_size = 20.0
)

opt_trajectory, system, Final_State = stabilize_grid!(
    [4, 4, 4],           # Grid dimensions
    [2.5, 3.5, 8.5],     # Spacing
    opt_params,
    harmonic_force_example,
    harmonic_potential_example;
    output_dir = "./"
)
final_geom = Final_State.Coords
# Save optimized geometry
write_xyz(final_geom, "$fold/lj_optimized.xyz")
##
p = plot( 1:size(opt_trajectory.PE,1), opt_trajectory.PE)#,xlabel="Steps",ylabel="PE"
savefig(p,"$fold/LJ_opt_PE.png")
##
println("\nOptimized geometry saved to $fold/lj_optimized.xyz")

write_trajectory(opt_trajectory,"$fold/lj_optimization_traj.xyz")

t_end = Final_State.t
Plot_data(opt_trajectory,1:1:500; smooth_flag=false,fold=fold,window_smooth=40)
