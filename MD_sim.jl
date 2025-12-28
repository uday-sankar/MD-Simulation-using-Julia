using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra

include("Code/MD_Base.jl")
using .MD_Base
##
include("/Users/upm5017/Documents/Work/FORTRAN_Julia/Working files/H2CO_triplet.jl")
##
# =============================================================================
fold = "/Users/upm5017/Documents/Work/MD_sim/Testing/Version 5: Molecules Implem"
# Define force and potential functions
function V(xyz::Matrix{Float64})
    return ch2o_pes(xyz)
end
##
dx = 0.01
function Force(xyz::Matrix{Float64})
    shape_xyz = size(xyz)
    f = zeros(shape_xyz)
    dummy_xyz = deepcopy(xyz)
    for i in 1:shape_xyz[1]
        for j in 1:shape_xyz[2]
            Xn = copy(dummy_xyz)
            Xp = copy(dummy_xyz)
            Xn[i,j] += dx
            Xp[i,j] -= dx
            f[i,j] = ( V(Xp) - V(Xn) )/(2*dx)
        end
    end
    return f
end
##
xyz_guess = [-1.0 0.2 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0; 0.10 0.0 1.3]
##
# Initialize system on a grid
positions = xyz_guess
velocities = xyz_guess*0.0
n_particles = size(positions, 1)
masses = [1.0, 1.0, 12.0, 16.0]
boxsize = 5.0
##
Atoms = ["H", "H", "C", "O"]#["Ar" for i in 1:size(positions,1)]
# Create system
system = System(Atoms, Force,V,boxsize,false)
##
mass_mat = hcat(masses,masses,masses)
State_init = MD_Base.SyState(Atoms,positions,velocities,velocities*0.0,mass_mat,0.0,0.0,0.0)
##
calculate_forces!(State_init,system)
# Set simulation parameters
dt = 0.01
params = SimulationParams(n_steps=1000,dt=dt,output_freq=1, boundary_size=15)
##
fire_params = FIRE_params(dt=0.01,ddt=0.002,da=0.01)
trajectory, Final_State = run_damped_optimization_FIRE!(State_init,
    system, params, fire_params;
    verbose = true, tol=1e-4)

# Save results
write_trajectory(trajectory, "$fold/H2CO_FIRE_stab.xyz")
write_enhanced_trajectory(trajectory, "$fold/H2CO_FIRE_stab_detailed.xyz",
                         dt=params.dt)

# Visualize (if Plots is available)
try
    animate_trajectory(trajectory, filename="$fold/H2CO_FIRE_stab.mp4", 
                      box_size=params.boundary_size, frames=50)
    plot_snapshots(trajectory, [1, :mid, :end], 
                  filename="$fold/H2CO_stab_snapshot.png")
catch e
    println("Visualization skipped (Plots.jl may not be available)")
end

energies = trajectory.Tot_Energy
T = trajectory.Temperature
PE = trajectory.PE
##
t_end = Final_State.t
Plot_data(trajectory,0:dt:t_end+dt; smooth_flag=false,fold=fold,window_smooth=40,plt_tot_E=false)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
# =======================================================================================
##
##
params = SimulationParams(n_steps=1000,dt=0.01,output_freq=1, boundary_size=5)
##
com = center_of_mass(Final_State)
CoM = vcat(com,com,com,com)
Final_State.Coords = Final_State.Coords - CoM
Final_State.Vel = [ -1.0 0.0 -0.0; +1.0 -0.0 0.0; -1.0 0.0 0.0; -1.0 -0.0 0.0]
# Run simulation
trajectory_MD_run, Final_State_MD = run_simulation!(Final_State,
    system, params,
    boundary_type = :reflective,
    verbose = true
)
write_trajectory(trajectory_MD_run, "$fold/H2CO_MD_run.xyz")
