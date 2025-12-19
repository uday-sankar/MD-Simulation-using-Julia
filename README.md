# Molecular Dynamics using Julia
An Molecular dynamics simulation code written julia.

The code has three main type of objects (implmented in types.jl):
- System: Positions, velocity etc. all infomration about the system
- SimulationParams: No. of Steps, dt, output frequecy, boundary etc. Parameters for the simulation
- Trajectory: Positions during trajectory, Tot. E, Pote. E, Temperature etc. All trajectory data

System

Capabilities:
- Initiaze geometry as a grid in 3D
- Optimize geometries using dampened velocity verlet 
- Perform MD simulations for any user defined interatomic potential
- Outputs VMD comaptible files for vizualization

Limitations:
- Current implementations assumes all particles are identical 
- No explicit way to grab a forcce field, as of now
- No Thermostats, No PBC

