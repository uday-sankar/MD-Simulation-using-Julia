# =============================================================================
# File: src/io.jl
# =============================================================================

"""
    read_xyz(filename)

Read XYZ coordinates from file.
Expects format: atom_id  x  y  z

Returns: positions (N×3 matrix)
"""
function read_xyz(filename::String)
    positions = []
    
    open(filename, "r") do f
        for line in eachline(f)
            columns = split(line)
            if length(columns) >= 4
                x = parse(Float64, columns[2])
                y = parse(Float64, columns[3])
                z = parse(Float64, columns[4])
                push!(positions, [x, y, z])
            end
        end
    end
    
    return hcat(positions...)' |> Matrix{Float64}
end

"""
    write_xyz(positions, filename; atom_types="Ar")

Write positions to XYZ file.

# Arguments
- `positions`: N×3 matrix of positions
- `filename`: Output filename
- `atom_types`: String or Vector of atom type labels
"""
function write_xyz(positions::Matrix{Float64}, filename::String; 
                   atom_types="Ar")
    
    n_atoms = size(positions, 1)
    
    # Handle atom types
    if typeof(atom_types) == String
        types = fill(atom_types, n_atoms)
    else
        types = atom_types
        if length(types) != n_atoms
            error("atom_types length must match number of atoms")
        end
    end
    
    open(filename, "w") do f
        write(f, "$n_atoms\n \n")
        for i in 1:n_atoms
            @printf(f, "%s %.6f %.6f %.6f\n", 
                   types[i], positions[i, 1], positions[i, 2], positions[i, 3])
        end
    end
    
    println("Geometry written to $filename")
end

"""
    write_trajectory(trajectory, filename; atom_types="Ar")

Write trajectory in XYZ format for visualization (e.g., VMD).

XYZ format:
```
<number of atoms>
comment line (frame info)
<atom> <x> <y> <z>
...
```

# Arguments
- `trajectory`: Array (n_frames, n_atoms, 3)
- `filename`: Output filename
- `atom_types`: String or Vector of atom type labels
"""
function write_trajectory(trajectory::Array{Float64, 3}, filename::String;
                         atom_types="Ar")
    
    n_frames = size(trajectory, 1)
    n_atoms = size(trajectory, 2)
    
    # Handle atom types
    if typeof(atom_types) == String
        types = fill(atom_types, n_atoms)
    else
        types = atom_types
        if length(types) != n_atoms
            error("atom_types length must match number of atoms")
        end
    end
    
    println("Writing trajectory to $filename")
    println("  Frames: $n_frames")
    println("  Atoms: $n_atoms")
    
    open(filename, "w") do f
        for frame in 1:n_frames
            # Write number of atoms
            write(f, "$n_atoms\n")
            
            # Write comment line
            write(f, "Frame $frame\n")
            
            # Write coordinates
            for atom in 1:n_atoms
                @printf(f, "%-4s  %12.6f  %12.6f  %12.6f\n",
                       types[atom],
                       trajectory[frame, atom, 1],
                       trajectory[frame, atom, 2],
                       trajectory[frame, atom, 3])
            end
        end
    end
    
    println("Trajectory written successfully!")
end

"""
    write_enhanced_trajectory(trajectory, filename;
                             atom_types="Ar",
                             energies=nothing,
                             temperatures=nothing,
                             dt=nothing)

Write trajectory with additional metadata in comment lines.

# Additional Arguments
- `energies`: Vector of energies for each frame
- `temperatures`: Vector of temperatures for each frame
- `dt`: Timestep size
"""
function write_enhanced_trajectory(trajectory::Array{Float64, 3}, filename::String;
                                  atom_types="Ar",
                                  energies=nothing,
                                  temperatures=nothing,
                                  dt=nothing)
    
    n_frames = size(trajectory, 1)
    n_atoms = size(trajectory, 2)
    
    # Handle atom types
    if typeof(atom_types) == String
        types = fill(atom_types, n_atoms)
    else
        types = atom_types
    end
    
    println("Writing enhanced trajectory to $filename")
    
    open(filename, "w") do f
        for frame in 1:n_frames
            write(f, "$n_atoms\n")
            
            # Enhanced comment line
            comment = "Frame $frame"
            if dt !== nothing
                time = (frame - 1) * dt
                comment *= @sprintf(" | Time: %.4f", time)
            end
            if energies !== nothing && frame <= length(energies)
                comment *= @sprintf(" | E: %.6f", energies[frame])
            end
            if temperatures !== nothing && frame <= length(temperatures)
                comment *= @sprintf(" | T: %.6f", temperatures[frame])
            end
            write(f, "$comment\n")
            
            # Write coordinates
            for atom in 1:n_atoms
                @printf(f, "%-4s  %12.6f  %12.6f  %12.6f\n",
                       types[atom],
                       trajectory[frame, atom, 1],
                       trajectory[frame, atom, 2],
                       trajectory[frame, atom, 3])
            end
        end
    end
    
    println("Enhanced trajectory written successfully!")
end

