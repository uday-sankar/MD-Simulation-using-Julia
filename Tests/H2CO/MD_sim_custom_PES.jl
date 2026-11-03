#using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra
#using Unitful
import USM_MD as MD_Base

println("Current Directoty",pwd())
##
include("../../FORTRAN_Test_PES/H2CO_triplet.jl")# custom force for H2CO system
##
# =============================================================================
Outdir = "FIRE_stab" # Path to Current folder (or the folder to keep all the information)
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
system = MD_Base.MD_System(Atoms,Force,V,boxsize,false)
##
mass_mat = copy(masses)#hcat(masses,masses,masses)
State_init = MD_Base.SyState(Atoms,positions,Vel = velocities,force = velocities*0.0,M = mass_mat)
##
MD_Base.calculate_forces!(State_init)
# Set simulation parameters
params = MD_Base.SimulationParams(n_steps=1000,dt=dt,output_freq=1,de_norm=1e-3,f_norm=1e-2,f_max=5e-2,Out_Dir=Outdir)
##
fire_params = MD_Base.FIRE_params(dt=0.1,ddt_inc=1.2,ddt_dec=0.5,dt_max=0.8,da=0.9,a=0.1,a0=0.1,dt_min=0.001)
trajectory, Final_State = MD_Base.run_damped_optimization_FIRE!(State_init,
    system, params, fire_params;
    verbose = true)

# Save results
MD_Base.write_trajectory(trajectory, "$Outdir/H2CO_FIRE_stab.xyz")

energies = trajectory.Tot_Energy
T = trajectory.Temperature
PE = trajectory.PE
##
t_end = Final_State.t
MD_Base.Plot_data(trajectory,0:dt:t_end+dt; smooth_flag=false,fold=Outdir,window_smooth=40,plt_tot_E=false)
##
println("\nMinimization Final Statistics:")
println("  Energy Change: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
# ======================MD simulation with an arbitary energy ================================================
# The following code is not excat. A code to do properly smapling of all vibrational modes is yet to be integrated into my MD suite.
# So the energy is dritubuted by keeping the center of mass fixed with finate rotational energy. (Check Codes/Initilizations.jl)
##
# Assiging velocity
ET = 1.0#Total kinetic energy, in eV, to be distributed as K.E.
New_MD_state = MD_Base.Recenter_and_vel(Final_State,ET)
##
Outdir = "MD_run" # Path to Current folder (or the folder to keep all the information)
##
params = MD_Base.SimulationParams(n_steps=1000,dt=0.02,output_freq=1,Out_Dir=Outdir)
# Run simulation
trajectory_MD_run, Final_State_MD = MD_Base.run_simulation!(New_MD_state,
    system, params,
    boundary_type = :reflective
)
##
energies = trajectory_MD_run.Tot_Energy
T = trajectory_MD_run.Temperature
PE = trajectory_MD_run.PE
println("\nMD run Final Statistics:")
println("  Energy Change: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
##
MD_Base.write_trajectory(trajectory_MD_run, "$Outdir/H2CO_MD_run.xyz")


MD_Base.Plot_data(trajectory_MD_run,1:size(trajectory_MD_run.PE,1); smooth_flag=false,window_smooth=40,plt_tot_E=true,fold=Outdir)
