# =============================================================================
# File: examples/run_simulation.jl
# Example usage of the MDSimulation module
# =============================================================================
#using Pkg
#using Revise
using Plots 
using Plots  # Optional, for visualization
using Statistics
using LinearAlgebra

include("Code/MD_Base.jl")
using .MD_Base
# =============================================================================
# Example 1: Simple Harmonic Oscillator System
# =============================================================================
fold = "Testing" # cuurent folder of interest

println("="^70)
println("Example 1: Harmonic Oscillator")
println("="^70)

# Define force and potential functions
function harmonic_force_example(r_vec)
    return MD_Base.harmonic_force(r_vec; k=1.0, r_eq=2.0)
end

function harmonic_potential_example(r_vec)
    return MD_Base.harmonic_potential(r_vec; k=1.0, r_eq=2.0)
end

# Initialize system on a grid
positions, velocities = grid_initialize([3, 2, 3], [2.0, 2.0, 2.0])
n_particles = size(positions, 1)
masses = ones(n_particles)

# Create system
system = System(positions, velocities, masses, 15.0,harmonic_force_example,harmonic_potential_example)

# Set simulation parameters
params = SimulationParams(n_steps=1000,dt=0.05,output_freq=1, boundary_size=15)
##

# Run simulation
trajectory, energies, temps, PE = run_simulation!(
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

println("\nFinal Statistics:")
println("  Energy drift: $(abs(energies[end] - energies[1]))")
println("  Mean temperature: $(sum(temps)/length(temps))")


# =============================================================================
# Example 2: Lennard-Jones Cluster Optimization
# =============================================================================

println("\n" * "="^70)
println("Example 2: Lennard-Jones Cluster Optimization")
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
##
p = plot( 1:size(opt_trajectory.PE,1), opt_trajectory.PE)#,xlabel="Steps",ylabel="PE"
savefig(p,"$fold/LJ_opt_PE.png")
##
println("\nOptimized geometry saved to $fold/lj_optimized.xyz")

write_trajectory(opt_trajectory,"$fold/lj_optimization_traj.xyz")

# =============================================================================
# Example 3: Production MD Run with LJ Potential
# =============================================================================

println("\n" * "="^70)
println("Example 3: Production MD with LJ Potential")
println("="^70)

# Use the optimized geometry as starting point
# Add some initial velocities
n_atoms = size(final_geom, 1)
initial_velocities = 0.1 * ones(n_atoms, 3)#randn(n_atoms, 3)  # Random thermal velocities

# Create new system
lj_system = System(final_geom, initial_velocities, ones(n_atoms), 10.0,lj_force_example,
    lj_potential_example)

# Production run parameters
prod_params = SimulationParams(
    n_steps = 10000,
    dt = 0.05,
    output_freq = 5,
    boundary_size = 10.0
)

# Run production simulation
lj_traj, lj_E, lj_T, lj_PE = run_simulation!(
    lj_system, prod_params;
    boundary_type = :reflective,
    verbose = true
)

# Save trajectory
write_enhanced_trajectory(lj_traj, "$fold/lj_production.xyz",
                         dt=prod_params.dt)

# Analysis
println("\nProduction Run Statistics:")
println("  Total steps: $(prod_params.n_steps)")
println("  Energy drift: $(abs(lj_E[end] - lj_E[1]))")
println("  Mean temperature: $(sum(lj_T)/length(lj_T)) ± $(std(lj_T))")
println("  Mean energy: $(sum(lj_E)/length(lj_E)) ± $(std(lj_E))")

# Calculate additional properties
println("\nFinal Configuration:")
println("  Radius of gyration: $(radius_of_gyration(lj_system))")
println("  Mean interparticle distance: $(mean_interparticle_distance(lj_system))")

# Create animation and analysis plots
#try
    animate_trajectory(lj_traj, filename="$fold/lj_production.gif",
                      box_size=20.0, frames=100,
                      trail_length=30, show_bonds=true, bond_cutoff=1.5)
    
    analyze_cluster_dynamics(lj_traj, lj_E, lj_T, lj_PE,
                            filename_prefix="lj_analysis")
#catch e
#    println("Visualization skipped: $e")
#end


# =============================================================================
# Example 4: Random Initialization and Optimization
# =============================================================================

println("\n" * "="^70)
println("Example 4: Random Initialization")
println("="^70)
opt_params = SimulationParams(
    n_steps = 5000,
    dt = 0.001,
    output_freq = 50,
    boundary_size = 20.0
)
# Random stabilization
random_geom, random_traj = stabilize_random!(
    20,                  # Number of particles
    10.0,                # Box size
    opt_params,
    lj_force_example,
    lj_potential_example
)

write_xyz(random_geom, "$fold/random_optimized.xyz")
println("Random optimized geometry saved")


# =============================================================================
# Example 5: Custom Potential
# =============================================================================

println("\n" * "="^70)
println("Example 5: Custom Potential")
println("="^70)

# Define custom Morse potential
function morse_potential(r_vec; D=1.0, a=1.0, r0=1.5)
    r = norm(r_vec)
    if r < 1e-10
        return 0.0
    end
    exp_term = exp(-a * (r - r0))
    return D * (1.0 - exp_term)^2
end

function morse_force(r_vec; D=1.0, a=1.0, r0=1.5)
    r = norm(r_vec)
    if r < 1e-10
        return zeros(3)
    end
    r_hat = r_vec / r
    exp_term = exp(-a * (r - r0))
    force_mag = 2.0 * D * a * (1.0 - exp_term) * exp_term / r
    return force_mag * r_vec
end

# Run with Morse potential
morse_positions, morse_velocities = grid_initialize([2, 2, 2], [2.0, 2.0, 2.0])
morse_system = System(morse_positions, morse_velocities, 
                      ones(size(morse_positions, 1)), 10.0, morse_force,morse_potential)

morse_params = SimulationParams(n_steps=1000, dt=0.01, output_freq=10, boundary_size=10.0)

morse_traj, morse_E, morse_T, morse_PE = run_simulation!(
    morse_system, morse_params,;
    verbose = true
)

write_trajectory(morse_traj, "$fold/morse_trajectory.xyz")
println("Morse potential simulation complete")


# =============================================================================
# Example 6: Analysis Only (Load from File)
# =============================================================================

println("\n" * "="^70)
println("Example 6: Post-Processing Analysis")
println("="^70)

# This example shows how to analyze an existing trajectory
# Assuming you have saved a trajectory earlier

# Read geometry
# existing_geom = read_xyz("lj_optimized.xyz")
# println("Loaded $(size(existing_geom, 1)) atoms from file")

# You could then create a system and analyze it
# analysis_system = System(existing_geom, zeros(size(existing_geom)), 
#                          ones(size(existing_geom, 1)), 20.0)

# println("Center of mass: $(center_of_mass(analysis_system))")
# println("Radius of gyration: $(radius_of_gyration(analysis_system))")


println("\n" * "="^70)
println("All examples completed successfully!")
println("="^70)
println("\nGenerated files:")
println("  - harmonic_trajectory.xyz")
println("  - harmonic_detailed.xyz")
println("  - lj_optimized.xyz")
println("  - lj_production.xyz")
println("  - random_optimized.xyz")
println("  - morse_trajectory.xyz")
println("  - Various .gif and .png files (if Plots.jl is available)")


# =============================================================================
# Quick Reference Guide
# =============================================================================

"""
QUICK REFERENCE FOR MDSimulation MODULE

1. INITIALIZATION
   - grid_initialize([Nx,Ny,Nz], [dx,dy,dz])    # Grid setup
   - random_initialize(n_particles, box_size)    # Random setup

2. SYSTEM CREATION
   - System(positions, velocities, masses, box_size)

3. SIMULATION PARAMETERS
   - SimulationParams(n_steps, dt, output_freq, boundary_size)

4. FORCE FUNCTIONS (examples included)
   - harmonic_force, harmonic_potential
   - lennard_jones_force, lennard_jones_potential
   - Custom: define your own!

5. RUNNING SIMULATIONS
   - run_simulation!(system, params, force_func, potential_func)
   - run_damped_optimization!(system, params, force_func, potential_func)
   - stabilize_grid!(n_grid, spacing, params, force_func, potential_func)
   - stabilize_random!(n_particles, box_size, params, force_func, potential_func)

6. ANALYSIS
   - temperature(system)
   - kinetic_energy(system), potential_energy(system), total_energy(system)
   - center_of_mass(system)
   - radius_of_gyration(system)
   - mean_interparticle_distance(system)

7. I/O
   - read_xyz(filename)
   - write_xyz(positions, filename)
   - write_trajectory(trajectory, filename)
   - write_enhanced_trajectory(trajectory, filename; energies, temps, dt)

8. VISUALIZATION (requires Plots.jl)
   - plot_frame(positions)
   - animate_trajectory(trajectory)
   - plot_snapshots(trajectory, frames)
   - analyze_cluster_dynamics(trajectory, E, T, PE)
"""