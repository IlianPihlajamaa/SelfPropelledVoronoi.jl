# This file is for my personal testing and development of the SelfPropelledVoronoi.jl package.

import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, StaticArrays, Random
using Statistics, BenchmarkTools, GLMakie


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

    function compute_mqd(displacement_array)
        msd = 0.0
        for i in eachindex(displacement_array)
            msd += displacement_array[i][1]^4 + displacement_array[i][2]^4
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
    N = length(arrays.positions)
    energy_list = Observable(Float64[])
    mean_squared_displacement_list = Observable(Float64[])
    mean_quartic_displacement_list = Observable(Float64[])
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
    ax1 = Axis(fig[1:3, 1], title="configuration", limits=(0, Lx, 0, Ly), xlabel="x", ylabel="y")
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

    tmin = dt
    tmax = @lift(
        if isempty($t_array)
            2dt
        else
            maximum($t_array)
        end
    )
    ymin = @lift(
        if isempty($mean_squared_displacement_list)
            1e-6
        else
            minimum($mean_squared_displacement_list)*0.9
        end
    )
    ymax = @lift(
        if isempty($mean_squared_displacement_list)
            1e2
        else
            maximum($mean_squared_displacement_list)*1.1
        end
    )

    ax4 = Axis(fig[2, 2], title="msd", ylabel="msd", xlabel="t", yscale=log10, xscale=log10, limits=@lift(($tmin, $tmax, $ymin, $ymax)))
    lines!(ax4, @lift($t_array[2:end]), mean_squared_displacement_list, color=:blue, label="msd")
    display(fig)
    α2 = @lift 0.5 * ($mean_quartic_displacement_list ./ $mean_squared_displacement_list .^ 2) .- 1
    ax3 = Axis(fig[3, 2], title="α2", ylabel="α2", xlabel="t", xscale=log10, limits=(0.01, 1000.0, -1.0, 1.0))
    lines!(ax3, @lift($t_array[2:end]), α2, color=:blue, label="α2")

    display(fig)

    function visualize(parameters, arrays, output)
        push!(energy_list[], SelfPropelledVoronoi.compute_energy(parameters, arrays, output))

        # apply periodic boundary conditions
        if !(length(energy_list[]) == 1) # no movement yet

            dr = arrays.positions .- previous_positions
            apply_periodic_boundary_conditions!(dr, parameters.box.box_sizes) 
            displacement_array .+= dr
            push!(mean_squared_displacement_list[], compute_msd(displacement_array))
            push!(mean_quartic_displacement_list[], compute_mqd(displacement_array))
        end

        α2new = 0.5 * (mean_quartic_displacement_list[] ./ mean_squared_displacement_list[] .^ 2) .- 1
        α2.val = α2new
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

logN = 1.6:0.2:4
run_times = Float64[]
for  Nappr = 10 .^ logN
    N = round(Int, Nappr)
    println("N = ", N)
    rho = 1.0
    L = sqrt(N/rho)
    Lx = L
    Ly = L
    dt = 0.05
    pbc_layer_depth = 2.5

    # Create a box
    box = SimulationBox(Lx, Ly)
    # Create a VoronoiCells object
    target_perimeters = 3.8*ones(N)
    target_areas = ones(N)
    K_P = ones(N)
    K_A = ones(N)
    active_force_strengths = ones(N)*0.01
    D_r = ones(N)
    voronoi_cells = VoronoiCells(
        target_perimeters,
        target_areas,
        K_P,
        K_A,
        active_force_strengths,
        D_r
    )

    # Create a ParameterStruct object
    kBT = 1.0
    frictionconstant = 1.0
    random_seed = 5435
    Random.seed!(random_seed)
    dump_info = DumpInfo(
        save=false,
        filename="Data/dump_$(random_seed).h5",
        when_to_save_array=0:1000:1000000,
    )

    rng = Random.MersenneTwister(random_seed)
    verbose=false
    cb(x...) = nothing 
    parameter_struct = ParameterStruct(N = N, dt = dt, 
        kBT = kBT, frictionconstant = frictionconstant, 
        periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
        dump_info = dump_info, callback = cb, RNG = rng)


    # Create an ArrayStruct object
    arrays = ArrayStruct(N)

    arrays.positions .= [rand(SVector{2, Float64}) .* box.box_sizes for _ in 1:N]
    #put particles on cubic lattice

    arrays.orientations .= 2π*rand(Float64, N) 

    # Create an Output object
    output = Output()
    arrays.neighborlist.check_tesselation = true
    # Run the simulation

    Nsteps = 100 ÷ dt
    println("Running simulation for N = ", N)
    trun = @elapsed run_simulation!(parameter_struct, arrays, output, Nsteps)
    trun_per_step = trun / Nsteps
    push!(run_times, trun_per_step)
    println("Run time for N = $N: ", trun_per_step, " seconds per step")
    @show output.N_voronoi_tesselations
end

fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1], title="Run time per step vs N", xlabel="N", ylabel="Run time (seconds per step)", xscale=log10, yscale=log10, limits=(60, 20000, 0.0001, 10))
scatterlines!(ax, 10 .^ logN, run_times, color=:blue, label="My code")

# extracted from paper
CellGPUtimes = [70.82071603229362	0.0003319015273144204
228.67384926584032	0.0010138931166180027
1134.3036212658983	0.005115239553097779
6317.65386596967	0.027982303949282764
33093.72671689607	0.14117482957733088
101878.95042348774	0.44544867545171346
] 

lines!(ax, CellGPUtimes[:, 1], CellGPUtimes[:, 2], color=:black, label="CellGPU times (CPU)", linestyle=:solid)

SPV2dtimes = [60.54647256691812	0.0018156274641954692
236.5999493550163	0.002765501715090403
763.960380327748	0.006111982405245609
3468.2346710721613	0.029853717700365313
13369.495472949164	0.15557137861690418
44665.20714587935	0.8107015041703827
]

lines!(ax, SPV2dtimes[:, 1], SPV2dtimes[:, 2], color=:gray, label="2dSPV times", linestyle=:dash)
axislegend(ax, position=:lt, fontsize=10)
display(fig)
save("example/run_time_vs_N.png", fig)

# visualize = create_visualization_callback(arrays, 1000, dt, box)

# # with a visualization callback
# parameter_struct2 = ParameterStruct(N = N, dt = dt, 
#     kBT = kBT, frictionconstant = frictionconstant, 
#     periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
#     dump_info = dump_info, callback = visualize, RNG = rng)

# Nsteps = 100000 ÷ dt
# @profview run_simulation!(parameter_struct2, arrays, output, Nsteps)
1+1