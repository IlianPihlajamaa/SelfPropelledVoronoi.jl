
"""
    create_visualization_callback(arrays, plot_frequency, dt, box)

Creates a visualization callback for the simulation.
The callback will visualize the simulation with the frequency set by the second argument, plotting the positions of the particles, the energy, and the mean squared displacement (MSD).

The visualization will be displayed in a figure with two axes: one for the particle configuration and one for the energy and MSD.
The first axis shows the particle positions and Voronoi edges, while the second axis shows the energy and MSD over time.

# Arguments
- `arrays`: The `ArrayStruct` containing the particle positions and orientations.
- `plot_frequency`: The frequency at which to plot the visualization.
- `dt`: The time step of the simulation.
- `box`: The `SimulationBox` containing the box sizes for periodic boundary conditions.

# Returns
 - A callback function that takes `parameters`, `arrays`, and `output` as arguments and updates the visualization.

# Notes
- The printing is necessary to ensure that the visualization updates correctly during the simulation.
"""
function create_visualization_callback(arrays,  plot_frequency, dt, box)
    function compute_msd(displacement_array)
        msd = 0.0
        for i in eachindex(displacement_array)
            msd += displacement_array[i][1]^2 + displacement_array[i][2]^2
        end
        return msd / length(displacement_array)
    end

    function apply_periodic_boundary_conditions!(displacement_array, box_sizes)
        for i in eachindex(displacement_array)
            dr = displacement_array[i]
            dr -= round.(dr ./ box_sizes) .* box_sizes
            displacement_array[i] = dr
        end
    end

    # This function creates a callback that will be called during the simulation
    # It will visualize the simulation at each step
    energy_list = Observable(Float64[])
    mean_squared_displacement_list = Observable(Float64[])
    displacement_array = [SVector{2, Float64}(0.0, 0.0) for _ in 1:N]
    previous_positions = copy(arrays.positions)
    t_array = @lift dt*eachindex($energy_list)

    plot_positions = Observable(SVector{2, Float64}[])
    for i in 1:N
        push!(plot_positions.val, previous_positions[i])
    end

    plot_voronoi_edges = Observable(Tuple{SVector{2, Float64}, SVector{2, Float64}}[])

    linesegs2 = Tuple{SVector{2, Float64}, SVector{2, Float64}}[]
    for particle_i in 1:length(arrays.positions)
        indices = arrays.neighborlist.voronoi_vertex_indices[particle_i]
        vor_positions = arrays.neighborlist.voronoi_vertices[indices]
        for i in 1:length(indices)
            j = i % length(indices) + 1
            posi = vor_positions[i]
            posj = vor_positions[j]
            # apply periodic boundary conditions

            push!(linesegs2, (posi, posj))
        end
    end
    plot_voronoi_edges.val = linesegs2


    Lx, Ly = box.box_sizes

    fig = Figure(size=(1000,500))
    ax1 = Axis(fig[1:2, 1], title="configuration", limits=(0, Lx, 0, Ly), xlabel="x", ylabel="y")
    scatter!(ax1, plot_positions, markersize=12, color =:blue)
    linesegments!(ax1, plot_voronoi_edges, color=:black)

    linesegments!(ax1, [
        (0, 0), (Lx, 0),
        (Lx, 0), (Lx, Ly),
        (Lx, Ly), (0, Ly),
        (0, Ly), (0, 0)
    ], color=:black, linewidth=4)
    
    max_energy = @lift(
        if isempty($energy_list)
            1.0
        else
            maximum($energy_list)*1.1
        end
    )

    min_energy = @lift(
        if isempty($energy_list)
            0.0
        else
            minimum($energy_list)*0.9
        end
        )


    xmax_energy = @lift(
        if isempty($t_array)
            1.0
        else
            maximum($t_array)
        end
    )
    ax2 = Axis(fig[1, 2], title="energy", ylabel="energy", xlabel="t", limits=@lift((0, $xmax_energy, $min_energy, $max_energy)))
    lines!(ax2, t_array, energy_list, color=:blue)

    ax4 = Axis(fig[2, 2], title="msd", ylabel="msd", xlabel="t", yscale=log10, xscale=log10, limits=(dt, 1e5, 1e-6, 1e2))
    lines!(ax4, @lift($t_array[2:end]), mean_squared_displacement_list, color=:blue, label="msd")
    display(fig)

    function visualize(parameters, arrays, output)

        # apply periodic boundary conditions
        push!(energy_list.val, SelfPropelledVoronoi.compute_energy(parameters, arrays, output))
        if !(length(energy_list[]) == 1) # no movement yet
            dr = arrays.positions .- previous_positions
            apply_periodic_boundary_conditions!(dr, parameters.box.box_sizes) 
            displacement_array .+= dr
            push!(mean_squared_displacement_list.val, compute_msd(displacement_array))
        end
        previous_positions .= arrays.positions

        if !(output.steps_done % plot_frequency == 0) 
            return
        end

        if SelfPropelledVoronoi.verify_tessellation(parameters, arrays, output) == false
            SelfPropelledVoronoi.voronoi_tesselation!(parameters, arrays, output)
        end

        # visualize the initial configuration
        println("Visualizing step $(output.steps_done)")

        energy_list[] = copy(energy_list[])
        mean_squared_displacement_list[] = copy(mean_squared_displacement_list[])

        plot_positions[] = SVector.(arrays.positions)

        # draw voronoi edges
        linesegs2 = Tuple{SVector{2, Float64}, SVector{2, Float64}}[]
        for particle_i in 1:length(arrays.positions)
            indices = arrays.neighborlist.voronoi_vertex_indices[particle_i]
            vor_positions = arrays.neighborlist.voronoi_vertices[indices]
            for i in 1:length(indices)
                j = i % length(indices) + 1
                posi = vor_positions[i]
                posj = vor_positions[j]
                # apply periodic boundary conditions

                push!(linesegs2, (posi, posj))
            end
        end
        plot_voronoi_edges[] = linesegs2
        
    end

    return (parameters, arrays, output) -> visualize(parameters, arrays, output)
end