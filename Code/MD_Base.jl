# Recommended file structure:
# 
# MDSimulation/
# ├── src/
# │   ├── MDSimulation.jl      # Main module file
# │   ├── types.jl             # Data structures
# │   ├── integrators.jl       # Time integration (velocity-verlet, langevin, etc.)
# │   ├── forces.jl            # Force calculations
# │   ├── initialization.jl    # System initialization functions
# │   ├── analysis.jl          # Analysis tools
# │   ├── io.jl                # Input/output functions
# │   └── visualization.jl     # Plotting and animation
# └── examples/
#     └── run_simulation.jl    # Example usage

# =============================================================================
# File: src/MDSimulation.jl (Main module file)
# =============================================================================

module MD_Base
using Plots  
using LinearAlgebra
using Printf

# Export types
export System, SimulationParams, FIRE_params

# Export integrators
export velocity_verlet!, run_simulation!, run_damped_optimization!, run_damped_optimization_FIRE!

# Export force calculation
export calculate_forces!

# Export initialization
export grid_initialize, random_initialize
export stabilize_grid!, stabilize_random!

# Export analysis
export kinetic_energy, potential_energy, total_energy
export temperature, center_of_mass
export radius_of_gyration, mean_interparticle_distance

# Export I/O
export read_xyz, write_xyz, write_trajectory
export write_enhanced_trajectory

# Export visualization (conditional on Plots)
export animate_trajectory, plot_frame, plot_snapshots
export analyze_cluster_dynamics
export Plot_data, smooth_boxcar

# Include all submodules
include("types.jl")
include("forces.jl")
include("Integrators.jl")
include("Initialization.jl")
include("Analysis.jl")
include("io.jl")
include("Visualization.jl")

end # module MDSimulation

