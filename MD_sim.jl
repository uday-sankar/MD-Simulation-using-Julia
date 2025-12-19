# =============================================================================
# File: examples/run_simulation.jl
# Example usage of the MDSimulation module
# =============================================================================

include("path/to/your/ModuleName.jl")

using .MD_Base
using Plots  # Optional, for visualization

# =============================================================================
# Example 1: Simple Harmonic Oscillator System
# =============================================================================

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
positions, velocities = grid_initialize([3, 3, 3], [2.0, 2.0, 2.0])
n_particles = size(positions, 1)
masses = ones(n_particles)

# Create system
system = System(positions, velocities, masses, 15.0)

# Set simulation parameters
params = SimulationParams(1000,0.01,10,15.0,)
##

# Run simulation
trajectory, energies, temps, PE = run_simulation!(
    system, params,
    harmonic_force_example,
    harmonic_potential_example;
    boundary_type = :reflective,
    verbose = true
)

# Save results
write_trajectory(trajectory, "harmonic_trajectory.xyz")
write_enhanced_trajectory(trajectory, "harmonic_detailed.xyz",
                         energies=energies,
                         temperatures=temps,
                         dt=params.dt)

# Visualize (if Plots is available)
try
    animate_trajectory(trajectory, filename="harmonic.gif", 
                      box_size=params.boundary_size, frames=50)
    plot_snapshots(trajectory, [1, :mid, :end], 
                  filename="harmonic_snapshots.png")
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
    return lennard_jones_force(r_vec; epsilon=1.0, sigma=1.0, cutoff=3.0)
end

function lj_potential_example(r_vec)
    return lennard_jones_potential(r_vec; epsilon=1.0, sigma=1.0, cutoff=3.0)
end

# Grid stabilization
opt_params = SimulationParams(
    n_steps = 5000,
    dt = 0.001,
    output_freq = 50,
    boundary_size = 20.0
)

final_geom, opt_trajectory = stabilize_grid!(
    [3, 3, 3],           # Grid dimensions
    [2.5, 2.5, 2.5],     # Spacing
    opt_params,
    lj_force_example,
    lj_potential_example;
    output_dir = "./"
)

# Save optimized geometry
write_xyz(final_geom, "lj_optimized.xyz")

println("\nOptimized geometry saved to lj_optimized.xyz")


# =============================================================================
# Example 3: Production MD Run with LJ Potential
# =============================================================================

println("\n" * "="^70)
println("Example 3: Production MD with LJ Potential")
println("="^70)

# Use the optimized geometry as starting point
# Add some initial velocities
n_atoms = size(final_geom, 1)
initial_velocities = 0.1 * randn(n_atoms, 3)  # Random thermal velocities

# Create new system
lj_system = System(final_geom, initial_velocities, ones(n_atoms), 20.0)

# Production run parameters
prod_params = SimulationParams(
    n_steps = 10000,
    dt = 0.005,
    output_freq = 50,
    boundary_size = 20.0
)

# Run production simulation
lj_traj, lj_E, lj_T, lj_PE = run_simulation!(
    lj_system, prod_params,
    lj_force_example,
    lj_potential_example;
    boundary_type = :reflective,
    verbose = true
)

# Save trajectory
write_enhanced_trajectory(lj_traj, "lj_production.xyz",
                         energies=lj_E,
                         temperatures=lj_T,
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
try
    animate_trajectory(lj_traj, filename="lj_production.gif",
                      box_size=20.0, frames=100,
                      trail_length=30, show_bonds=true, bond_cutoff=1.5)
    
    analyze_cluster_dynamics(lj_traj, lj_E, lj_T, lj_PE,
                            filename_prefix="lj_analysis")
catch e
    println("Visualization skipped: $e")
end


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

write_xyz(random_geom, "random_optimized.xyz")
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
                      ones(size(morse_positions, 1)), 10.0)

morse_params = SimulationParams(1000, 0.01, 10, 10.0)

morse_traj, morse_E, morse_T, morse_PE = run_simulation!(
    morse_system, morse_params,
    morse_force, morse_potential;
    verbose = true
)

write_trajectory(morse_traj, "morse_trajectory.xyz")
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