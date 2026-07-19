using PythonCall

np = pyimport("numpy")
tb = pyimport("tblite.interface")
##
using Revise
using Plots  # ]add Plots if needed
gr()         # Use GR backend for high-quality output # Optional, for visualization
using Statistics
using LinearAlgebra
#using Unitful

println("Current Directoty",pwd())

##
include("../../Code/MD_Base.jl") # Base code
using .MD_Base
##
const SYMBOL_TO_Z = Dict(
    "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7,
    "O" => 8, "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13,
    "Si" => 14, "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18,
    # extend as needed for your systems
)
##
atoms, coords = read_xyz("./Init_simple.xyz")
##
atomic_numbers = [SYMBOL_TO_Z[a] for a in atoms]
coords_bohr = coords# .* 1.8897259886
##
##
calc = tb.Calculator(
    "GFN2-xTB",
    np.array(atomic_numbers),
    np.array(coords_bohr),
    charge = 0,      # total molecular charge
    uhf    = 0,      # number of unpaired electrons (mult - 1)
)
calc.set("verbosity", 0)   # suppress the SCF cycle printout every call
##
# --- called every MD step ---
function energy_forces!(positions_bohr::Matrix{Float64})
    calc.update(positions = np.array(positions_bohr).* 1.8897259886)   # reuses prior density as SCF guess (faster + more stable than rebuilding Calculator)
    res = calc.singlepoint()
    #energy_hartree = res["energy"]
    energy_hartree = pyconvert(Float64, res["energy"].item())
    gradient = pyconvert(Matrix{Float64}, res["gradient"])   # dE/dR, Hartree/Bohr, shape (natoms,3)
    forces = -gradient                                       # F = -∇E
    return energy_hartree, forces
end
##
function Get_Energy(coords)
    calc.update(positions = np.array(coords).* 1.8897259886)   # reuses prior density as SCF guess (faster + more stable than rebuilding Calculator)
    res = calc.singlepoint()
    #energy_hartree = res["energy"]
    energy_hartree = pyconvert(Float64, res["energy"].item())
    #gradient = pyconvert(Matrix{Float64}, res["gradient"])   # dE/dR, Hartree/Bohr, shape (natoms,3)
    #forces = -gradient                                       # F = -∇E
    return energy_hartree
end
##
function Get_force(coords)
    calc.update(positions = np.array(coords.* 1.8897259886))   # reuses prior density as SCF guess (faster + more stable than rebuilding Calculator)
    res = calc.singlepoint()
    gradient = pyconvert(Matrix{Float64}, res["gradient"])   # dE/dR, Hartree/Bohr, shape (natoms,3)
    forces = -gradient                                       # F = -∇E
    return forces#./0.529
end
##
Get_Energy(coords_bohr)
Get_force(coords_bohr)
##
using PeriodicTable

function get_atomic_masses(atoms::Vector{String})
    return [elements[Symbol(a)].atomic_mass.val for a in atoms]
end
##
boxsize = 100.0
##
system = MD_System(atoms,Get_force,Get_Energy,boxsize,false)
##
mass_mat = get_atomic_masses(atoms)#copy(masses)#hcat(masses,masses,masses)
##
State_init = SyState(atoms,coords_bohr,Vel = coords_bohr*0.0,force = coords_bohr*0.0,M = mass_mat)
# Assiging velocity
ET = 0.005#Total kinetic energy to be distributed as K.E.
New_MD_state = Recenter_and_vel(State_init,ET)
##
params = SimulationParams(n_steps=1000,dt=0.1,output_freq=1, boundary_size=20)
# Run simulation
trajectory_MD_run, Final_State_MD = run_simulation!(New_MD_state,
    system, params,
    boundary_type = :reflective,
    verbose = true
)
fold = "Tests/Methane_O2"
write_trajectory(trajectory_MD_run, "$fold/MD_run.xyz")
##
Plot_data(trajectory_MD_run,1:size(trajectory_MD_run.PE,1); smooth_flag=false,window_smooth=40,plt_tot_E=true)