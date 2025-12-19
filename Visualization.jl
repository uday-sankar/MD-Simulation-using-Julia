
# =============================================================================
# File: src/visualization.jl
# =============================================================================

# Only load visualization if Plots is available
if isdefined(Main, :Plots) || (isdefined(Main, :__PLOTTING_LOADED__) && Main.__PLOTTING_LOADED__)
    using Plots
    
    """
        plot_frame(positions; filename="frame.png", masses=nothing, box_size=20.0)
    
    Plot a single frame of particle positions.
    """
    function plot_frame(positions::Matrix{Float64}; 
                       filename="frame.png",
                       masses=nothing,
                       box_size=20.0)
        
        n_atoms = size(positions, 1)
        
        if masses === nothing
            masses = ones(n_atoms)
        end
        
        # Create box
        bl = box_size
        plot([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], -bl.*ones(5),
             color=:black, linewidth=3, legend=false,
             xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl),
             xlabel="X", ylabel="Y", zlabel="Z")
        plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], bl.*ones(5),
             color=:black, linewidth=3)
        plot!(bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl],
             color=:black, linewidth=3)
        plot!(-bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl],
             color=:black, linewidth=3)
        
        # Plot particles
        for i in 1:n_atoms
            scatter!([positions[i, 1]], [positions[i, 2]], [positions[i, 3]],
                    color=:blue, markersize=2*masses[i])
        end
        
        savefig(filename)
        println("Frame saved to $filename")
    end
    
    """
        animate_trajectory(trajectory; filename="animation.gif", 
                          masses=nothing, frames=100, box_size=20.0,
                          trail_length=20, show_bonds=false, bond_cutoff=1.5)
    
    Create animated GIF of trajectory.
    """
    function animate_trajectory(trajectory::Array{Float64, 3};
                               filename="animation.gif",
                               masses=nothing,
                               frames=100,
                               box_size=20.0,
                               trail_length=20,
                               show_bonds=false,
                               bond_cutoff=1.5)
        
        n_timesteps = size(trajectory, 1)
        n_atoms = size(trajectory, 2)
        
        if masses === nothing
            masses = ones(n_atoms)
        end
        
        println("Creating animation: $filename")
        println("  Frames: $frames")
        println("  Timesteps: $n_timesteps")
        println("  Atoms: $n_atoms")
        
        bl = box_size
        step_size = max(1, n_timesteps ÷ frames)
        
        anim = @animate for j in 1:frames
            i = min(j * step_size, n_timesteps)
            
            # Draw box
            plot([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], -bl.*ones(5),
                 color=:black, linewidth=3, legend=false,
                 xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl),
                 xlabel="X", ylabel="Y", zlabel="Z",
                 title="Frame $i/$n_timesteps")
            plot!([-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl], bl.*ones(5),
                 color=:black, linewidth=3)
            plot!(bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl],
                 color=:black, linewidth=3)
            plot!(-bl.*ones(5), [-bl, bl, bl, -bl, -bl], [-bl, -bl, bl, bl, -bl],
                 color=:black, linewidth=3)
            
            # Draw bonds if requested
            if show_bonds
                for pi in 1:n_atoms
                    for pj in pi+1:n_atoms
                        r = norm(trajectory[i, pi, :] - trajectory[i, pj, :])
                        if r < bond_cutoff
                            plot!([trajectory[i, pi, 1], trajectory[i, pj, 1]],
                                  [trajectory[i, pi, 2], trajectory[i, pj, 2]],
                                  [trajectory[i, pi, 3], trajectory[i, pj, 3]],
                                  color=:gray, alpha=0.3, linewidth=1)
                        end
                    end
                end
            end
            
            # Draw particles
            for p in 1:n_atoms
                scatter!([trajectory[i, p, 1]], [trajectory[i, p, 2]], [trajectory[i, p, 3]],
                        color=:blue, markersize=2*masses[p])
            end
            
            # Draw trails
            trail_start = max(1, i - trail_length)
            if i > trail_start
                for p in 1:n_atoms
                    plot!(trajectory[trail_start:i, p, 1],
                          trajectory[trail_start:i, p, 2],
                          trajectory[trail_start:i, p, 3],
                          color=:lightblue, alpha=0.4, linewidth=1)
                end
            end
        end
        
        gif(anim, filename, fps=30)
        println("Animation saved to $filename")
    end
    
    """
        plot_snapshots(trajectory, frames; filename="snapshots.png",
                      masses=nothing, box_size=20.0, show_bonds=false,
                      bond_cutoff=1.5)
    
    Create snapshot comparison at specified frames.
    
    # Arguments
    - `frames`: Vector like [1, :mid, :end] or [1, 50, 100]
    """
    function plot_snapshots(trajectory::Array{Float64, 3}, frames;
                           filename="snapshots.png",
                           masses=nothing,
                           box_size=20.0,
                           show_bonds=false,
                           bond_cutoff=1.5)
        
        n_timesteps = size(trajectory, 1)
        n_atoms = size(trajectory, 2)
        
        if masses === nothing
            masses = ones(n_atoms)
        end
        
        # Convert frame specifications
        frame_indices = Int[]
        frame_labels = String[]
        
        for f in frames
            if f == :mid
                push!(frame_indices, n_timesteps ÷ 2)
                push!(frame_labels, "Middle")
            elseif f == :end
                push!(frame_indices, n_timesteps)
                push!(frame_labels, "Final")
            else
                push!(frame_indices, f)
                push!(frame_labels, "t=$f")
            end
        end
        
        n_frames = length(frame_indices)
        plots_array = []
        
        bl = box_size
        
        for (idx, i) in enumerate(frame_indices)
            p = plot(legend=false,
                    xlabel="X", ylabel="Y", zlabel="Z",
                    xlims=(-bl, bl), ylims=(-bl, bl), zlims=(-bl, bl),
                    title=frame_labels[idx],
                    camera=(30, 30))
            
            # Draw bonds
            if show_bonds
                for pi in 1:n_atoms
                    for pj in pi+1:n_atoms
                        r = norm(trajectory[i, pi, :] - trajectory[i, pj, :])
                        if r < bond_cutoff
                            plot!([trajectory[i, pi, 1], trajectory[i, pj, 1]],
                                  [trajectory[i, pi, 2], trajectory[i, pj, 2]],
                                  [trajectory[i, pi, 3], trajectory[i, pj, 3]],
                                  color=:gray, alpha=0.3, linewidth=1)
                        end
                    end
                end
            end
            
            # Draw particles
            scatter!(trajectory[i, :, 1], trajectory[i, :, 2], trajectory[i, :, 3],
                    color=:blue, markersize=10 .* masses)
            
            push!(plots_array, p)
        end
        
        combined = plot(plots_array..., layout=(1, n_frames), 
                       size=(400*n_frames, 400))
        savefig(combined, filename)
        println("Snapshots saved to $filename")
    end
    
    """
        analyze_cluster_dynamics(trajectory, energies, temperatures, potential_energies;
                                filename_prefix="analysis", masses=nothing)
    
    Create comprehensive analysis plots.
    """
    function analyze_cluster_dynamics(trajectory::Array{Float64, 3},
                                     energies::Vector{Float64},
                                     temperatures::Vector{Float64},
                                     potential_energies::Vector{Float64};
                                     filename_prefix="analysis",
                                     masses=nothing)
        
        n_frames = size(trajectory, 1)
        n_atoms = size(trajectory, 2)
        
        if masses === nothing
            masses = ones(n_atoms)
        end
        
        # Calculate radius of gyration
        Rg = zeros(n_frames)
        for i in 1:n_frames
            # Calculate CoM for this frame
            total_mass = sum(masses)
            com = zeros(3)
            for p in 1:n_atoms
                com .+= masses[p] .* trajectory[i, p, :]
            end
            com ./= total_mass
            
            # Calculate Rg
            rg_sq = 0.0
            for p in 1:n_atoms
                r_rel = trajectory[i, p, :] - com
                rg_sq += masses[p] * sum(r_rel.^2)
            end
            Rg[i] = sqrt(rg_sq / total_mass)
        end
        
        # Calculate mean interparticle distance
        mean_dist = zeros(n_frames)
        for i in 1:n_frames
            dists = Float64[]
            for pi in 1:n_atoms
                for pj in pi+1:n_atoms
                    r = norm(trajectory[i, pi, :] - trajectory[i, pj, :])
                    push!(dists, r)
                end
            end
            mean_dist[i] = sum(dists) / length(dists)
        end
        
        # Create plots
        steps = 1:n_frames
        
        p1 = plot(steps, energies, xlabel="Step", ylabel="Total Energy",
                 title="Energy Conservation", linewidth=2, legend=false)
        
        p2 = plot(steps, temperatures, xlabel="Step", ylabel="Temperature",
                 title="Temperature", linewidth=2, color=:red, legend=false)
        
        p3 = plot(steps, potential_energies, xlabel="Step", ylabel="PE",
                 title="Potential Energy", linewidth=2, color=:green, legend=false)
        
        p4 = plot(steps, Rg, xlabel="Step", ylabel="Rg",
                 title="Radius of Gyration", linewidth=2, color=:purple, legend=false)
        
        p5 = plot(steps, mean_dist, xlabel="Step", ylabel="⟨r⟩",
                 title="Mean Distance", linewidth=2, color=:orange, legend=false)
        
        p6 = histogram(energies, xlabel="Energy", ylabel="Count",
                      title="Energy Distribution", bins=50, legend=false)
        
        combined = plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(1200, 1200))
        savefig(combined, "$(filename_prefix)_dynamics.png")
        
        println("\nCluster Analysis:")
        println("  Initial Rg: $(Rg[1])")
        println("  Final Rg: $(Rg[end])")
        println("  Energy drift: $(abs(energies[end] - energies[1]))")
        println("  Mean temperature: $(sum(temperatures)/length(temperatures))")
        println("  Plots saved to $(filename_prefix)_dynamics.png")
        
        return Rg, mean_dist
    end
    
else
    # Provide stub functions if Plots is not available
    function plot_frame(args...; kwargs...)
        @warn "Plots.jl not loaded. Visualization functions unavailable."
    end
    
    animate_trajectory = plot_frame
    plot_snapshots = plot_frame
    analyze_cluster_dynamics = plot_frame
end