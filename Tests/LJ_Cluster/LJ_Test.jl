using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
import USM_MD as MD_Base
#using Statistics
#using LinearAlgebra
##
Out_put_fold = "Minimize"
if !isdir(Out_put_fold)
    mkdir(Out_put_fold)
end
##
# getting positions and veclocities 
# A grid based initlization of particles in a box with given dimensions and spacing
pos, vel = MD_Base.grid_initialize([4,4,4],[5.0,5.0,5.0]) # 4 particles along x,y,z with seperation of 2.0
num_part = size(pos,1)
Atoms = ["Au" for i in 1:num_part]
mass = ones(num_part)
## 
F(xyz) = MD_Base.lennard_jones_force(xyz; epsilon=2.0, sigma = 3.5, cutoff = Inf)
V(xyz) = MD_Base.lennard_jones_potential(xyz; epsilon=2.0, sigma = 3.5, cutoff = Inf)
box_size = 15.0
system = MD_Base.MD_System(Atoms,F,V,box_size,true)
##
State_init = MD_Base.SyState(Atoms, pos; Vel=vel, force = vel*0.0,M=mass)
##
dt = 0.1
out_freq = 1
N_steps = 2000
params = MD_Base.SimulationParams(n_steps=N_steps,dt=dt,output_freq=out_freq)
##
a, da = 0.5, 0.9
fire_params = MD_Base.FIRE_params(dt=dt,ddt_inc=1.1,ddt_dec=0.5,da=da,a=a,a0=0.5,dt_min=0.001)
##
opt_trajectory, opt_State = MD_Base.run_damped_optimization_FIRE!(State_init,
    system, params, fire_params;
    verbose = true, tol=1e-4)
## Save results
MD_Base.write_trajectory(opt_trajectory, "$Out_put_fold/LJ_stab.xyz")
##
#try
MD_Base.animate_trajectory(opt_trajectory, filename="$Out_put_fold/LJ_stab.mp4", box_size=system.box_size, frames=100)
#catch e
#    println("Visualization skipped (Plots.jl may not be available)")
#end
##
energies = opt_trajectory.Tot_Energy
T = opt_trajectory.Temperature
PE = opt_trajectory.PE
t_end = opt_State.t
##
#print("Temp Len: ",length(T))
#print("Time Len: ",length(0:dt:t_end+dt))
#
MD_Base.Plot_data(opt_trajectory,opt_trajectory.Time; smooth_flag=false,fold=Out_put_fold,window_smooth=40,plt_tot_E=false)
println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(T)/length(T))")
## MD Run
# Assiging velocity
ET = 100.0#Total kinetic energy to be distributed as K.E.
New_MD_state = MD_Base.Recenter_and_vel(opt_State,ET)
##
Out_put_fold = "MD"
if !isdir(Out_put_fold)
    mkdir(Out_put_fold)
end
##
N_steps = 40000
dt = 0.0005
##
z_vel = repeat([0.0,0.0,-2.0]',num_part,1)
New_MD_state.Vel = New_MD_state.Vel +  z_vel
out_freq = 10
#print(New_MD_state.Vel)
params = MD_Base.SimulationParams(n_steps=N_steps,dt=dt,output_freq=out_freq)
trajectory_MD_run, Final_State_MD = MD_Base.run_simulation!(New_MD_state,
    system, params,
    boundary_type = :reflective,
    verbose = true
)
##
MD_Base.write_trajectory(trajectory_MD_run, "$Out_put_fold/LJ_MD_run.xyz")
##
t_end = Final_State_MD.t
MD_Base.Plot_data(trajectory_MD_run,1:size(trajectory_MD_run.PE,1); smooth_flag=false,fold=Out_put_fold,window_smooth=40,plt_tot_E=true)
##
MD_Base.animate_trajectory(trajectory_MD_run, filename="$Out_put_fold/LJ_MD_run.mp4", 
                      box_size=system.box_size, frames=400)
