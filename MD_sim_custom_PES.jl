using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra
#using Unitful

include("./Code/MD_Base.jl") # Base code
using .MD_Base
##
include("./FORTRAN_Test_PES/H2CO_triplet.jl")# custom force for H2CO system
##
# =============================================================================
fold = "./" # Path to Current folder (or the folder to keep all the information)
# Define force and potential functions
function V(xyz::Matrix{Float64})# custom potential
    return ch2o_pes(xyz)#u"eV"
end
##
dx = 0.01#u"Å" Arbitary dx to evaluate numerical forces
dt = 0.01#u"fs" Arbitary time step
# Both can be varied to determine the optimal values which allow energy convergance
##
function Force(xyz::Matrix{Float64})# custom force 
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
xyz_guess = [-1.0 0.2 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0; 0.10 0.0 1.3] # Guess Geometry
##
# Initialize system on a grid
positions = xyz_guess
velocities = xyz_guess*0.0 ./dt
n_particles = size(positions, 1)
##
Atoms = ["H", "H", "C", "O"] # Atoms involved in the exact order as forces
masses = [1.0, 1.0, 12.0, 16.0] # Masses for each atom
boxsize = 5.0
##
system = MD_System(Atoms,Force,V,boxsize,false)
##
mass_mat = copy(masses)#hcat(masses,masses,masses)
State_init = SyState(Atoms,positions,Vel = velocities,force = velocities*0.0,M = mass_mat)
##
calculate_forces!(State_init)
# Set simulation parameters
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

# Visualize (if Plots is available). Current vizualization works with julia. Python Vizualization coming soon
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
Plot_data(trajectory,0:dt:t_end+dt; smooth_flag=false,filename="$fold/H2CO_stab_data.png",window_smooth=40,plt_tot_E=false)
##
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
# ======================MD simulation with an arbitary energy ================================================
# The following code is not excat. A code to do properly smapling of all vibrational modes is yet to be integrated into my MD suite.
# So the energy is dritubuted by keeping the center of mass fixed with finate rotational energy. (Check Codes/Initilizations.jl)
##
# Assiging velocity
ET = 8.0#Total kinetic energy to be distributed as K.E.
New_MD_state = Recenter_and_vel(Final_State,ET)
##
params = SimulationParams(n_steps=1000,dt=0.02,output_freq=1, boundary_size=3)
# Run simulation
trajectory_MD_run, Final_State_MD = run_simulation!(New_MD_state,
    system, params,
    boundary_type = :reflective,
    verbose = true
)
write_trajectory(trajectory_MD_run, "$fold/H2CO_MD_run.xyz")

t_end = Final_State_MD.t

Plot_data(trajectory_MD_run,1:size(trajectory_MD_run.PE,1); smooth_flag=false,filename="$fold/MD_run_data.png",window_smooth=40,plt_tot_E=true)
