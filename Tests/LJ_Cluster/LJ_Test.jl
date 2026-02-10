using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra
##
include("../../Code/MD_Base.jl")
using .MD_Base
##
Out_put_fold = "./"
##
# getting positions and veclocities
pos, vel = grid_initialize([6,6,6],[5.0,5.0,5.0]) # 5 particles along x,y,z with seperation of 4.0
num_part = size(pos,1)
Atoms = ["Au" for i in 1:num_part]
mass = ones(num_part)
## 
F(xyz) = MD_Base.lennard_jones_force(xyz; epsilon=2.0, sigma = 3.5, cutoff = Inf)
V(xyz) = MD_Base.lennard_jones_potential(xyz; epsilon=2.0, sigma = 3.5, cutoff = Inf)
box_size = 25.0
system = MD_System(Atoms,F,V,box_size,true)
##
State_init = SyState(Atoms, pos; Vel=vel, force = vel*0.0,M=mass)
##
dt = 0.01
out_freq = 1
N_steps = 2000
params = SimulationParams(n_steps=N_steps,dt=dt,output_freq=out_freq, boundary_size=box_size)
##
a, da = 1.0, 0.05
fire_params = FIRE_params(dt=dt,ddt=dt/10.0,da=da,a=a)
##
opt_trajectory, opt_State = run_damped_optimization_FIRE!(State_init,
    system, params, fire_params;
    verbose = true, tol=1e-4)
## Save results
write_trajectory(opt_trajectory, "$Out_put_fold/LJ_stab.xyz")
##
try
    animate_trajectory(opt_trajectory, filename="$Out_put_fold/LJ_stab.mp4", 
                      box_size=params.boundary_size, frames=100)
catch e
    println("Visualization skipped (Plots.jl may not be available)")
end
##
energies = opt_trajectory.Tot_Energy
T = opt_trajectory.Temperature
PE = opt_trajectory.PE
t_end = opt_State.t
Plot_data(opt_trajectory,0:dt:t_end; smooth_flag=false,fold=Out_put_fold,window_smooth=40,plt_tot_E=false)
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
## MD Run
# Assiging velocity
ET = 10.0#Total kinetic energy to be distributed as K.E.
New_MD_state = Recenter_and_vel(opt_State,ET)
##
z_vel = repeat([0.0,0.0,-4.0]',num_part,1)
New_MD_state.Vel = New_MD_state.Vel + z_vel
params = SimulationParams(n_steps=N_steps*2.0,dt=dt,output_freq=out_freq, boundary_size=box_size)
trajectory_MD_run, Final_State_MD = run_simulation!(New_MD_state,
    system, params,
    boundary_type = :reflective,
    verbose = true
)
##
write_trajectory(trajectory_MD_run, "$Out_put_fold/LJ_MD_run.xyz")
##
t_end = Final_State_MD.t
Plot_data(trajectory_MD_run,1:size(trajectory_MD_run.PE,1); smooth_flag=false,fold=Out_put_fold,window_smooth=40,plt_tot_E=true)
##
animate_trajectory(trajectory_MD_run, filename="$Out_put_fold/LJ_MD_run.mp4", 
                      box_size=params.boundary_size, frames=200)
