A molecular dynamics simulation code written in Julia.

MD_sim.jl : A sample MD code with LJ Potential
MD_sim_custom_PES.jl: A sample code which works with a custom user defined PES (Here it uses a FORTRAN PES for formaldehyde)

The code defines three main object types (implemented in types.jl):
- System: Contains positions, velocities, and all system-related information
- SimulationParams: Includes number of steps, timestep (dt), output frequency, boundaries, and other simulation parameters
- Trajectory: Stores positions during the simulation, total energy, potential energy, temperature, and all trajectory-related data

Capabilities:
- Initialize geometry as a grid in 3D
- Optimize geometries using FIRE or dampened velocity Verlet 
- Perform MD simulations for any user-defined interatomic potential or custom PES
- Outputs VMD-compatible files for visualization

Limitations:
- No explicit way to grab a force field. (The user can define a way to grab a force field from another library and use it with this code)
- No Thermostats, No PBC
- No explicit parallelization (The MD steps are vectorized natively in Julia (but not parallelized using MPI))

References:
The sample H2CO PES was obtained from the chempotpy library which used a PIP-NN from another publication:
PIP NN: https://pubs.rsc.org/en/content/articlelanding/2017/cp/c7cp04578f, 
Li, Jun, Changjian Xie, and Hua Guo. "Kinetics and dynamics of the C (3 P)+ H 2 O reaction on a full-dimensional accurate triplet state potential energy surface." Physical Chemistry Chemical Physics 19, no. 34 (2017): 23280-23288.

Chempotpy: https://pubs.acs.org/doi/full/10.1021/acs.jpca.3c05Shu, 
Yinan, Zoltan Varga, Dayou Zhang, and Donald G. Truhlar. "ChemPotPy: A Python Library for Analytic Representations of Potential Energy Surfaces and Diabatic Potential Energy Matrices." The Journal of Physical Chemistry A 127, no. 45 (2023): 9635-9640.899, 