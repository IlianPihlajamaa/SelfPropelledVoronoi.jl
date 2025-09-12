# This file is for my personal testing and development of the SelfPropelledVoronoi.jl package.

import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, StaticArrays, Random
using Statistics, BenchmarkTools, GLMakie
import Quickhull, SIMD


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

    plot_positions = Observable(SVector{2, Float64}[])
    for i in 1:arrays.neighborlist.N_positions_with_pbc
        push!(plot_positions.val, arrays.neighborlist.positions_with_pbc[i])
    end

    vertices = Observable(SVector{2, Float64}[])
    for i in 1:length(arrays.neighborlist.voronoi_vertices)
        push!(vertices.val, arrays.neighborlist.voronoi_vertices[i])
    end

    plot_voronoi_edges = Observable(Tuple{SVector{2, Float64}, SVector{2, Float64}}[])

    linesegs2 = Tuple{SVector{2, Float64}, SVector{2, Float64}}[]
    for particle_i in 1:length(arrays.positions)
        indices = arrays.neighborlist.voronoi_vertex_indices[particle_i]
        vor_positions = arrays.neighborlist.voronoi_vertices[indices]
        N_indices = arrays.neighborlist.N_voronoi_vertices_pp[particle_i]
        for i in 1:N_indices
            j = i % length(indices) + 1
            posi = vor_positions[i]
            posj = vor_positions[j]
            # apply periodic boundary conditions

            push!(linesegs2, (posi, posj))
        end
    end
    plot_voronoi_edges.val = linesegs2


    Lx, Ly = box.box_sizes

    fig = Figure(size=(1000,1000))
    ax1 = Axis(fig[1:3, 1], title="configuration", limits=(-5, Lx+5, -5, Ly+5), xlabel="x", ylabel="y")
    scatter!(ax1, plot_positions, markersize=12, color =:blue)
    text!(ax1, plot_positions, text = string.(1:length(arrays.positions)))

    linesegments!(ax1, plot_voronoi_edges, color=:black, linewidth=0.5)
    text!(ax1, vertices, text = string.(1:length(vertices[])), color=:red, align = (:left, :bottom), fontsize=8)
    linesegments!(ax1, [
        (0, 0), (Lx, 0),
        (Lx, 0), (Lx, Ly),
        (Lx, Ly), (0, Ly),
        (0, Ly), (0, 0)
    ], color=:black, linewidth=4)
    
    display(fig)

    function visualize(parameters, arrays, output; do_anyway=false)
        if !do_anyway
            if !(output.steps_done % plot_frequency == 0) 
                return
            end
        end

        if SelfPropelledVoronoi.verify_tessellation(parameters, arrays, output) == false
            SelfPropelledVoronoi.voronoi_tessellation!(parameters, arrays, output)
        end

        # visualize the initial configuration
        println("Visualizing step $(output.steps_done)")

        plot_positions[] = SVector.(arrays.neighborlist.positions_with_pbc[1:arrays.neighborlist.N_positions_with_pbc])
        vertices[] = SVector.(arrays.neighborlist.voronoi_vertices)
        # draw voronoi edges
        linesegs2 = Tuple{SVector{2, Float64}, SVector{2, Float64}}[]
        for particle_i in 1:arrays.neighborlist.N_positions_with_pbc
            n_indices = arrays.neighborlist.N_voronoi_vertices_pp[particle_i]
            indices = arrays.neighborlist.voronoi_vertex_indices[particle_i][1:n_indices]
            vor_positions = arrays.neighborlist.voronoi_vertices[indices]
            for i in 1:n_indices
                j = i % n_indices + 1
                posi = vor_positions[i]
                posj = vor_positions[j]

                push!(linesegs2, (posi, posj))
            end
        end
        plot_voronoi_edges[] = linesegs2

    end

    return (parameters, arrays, output; do_anyway=false) -> visualize(parameters, arrays, output; do_anyway=do_anyway)
end





N = 50
rho = 1.0
L = sqrt(N/rho)
Lx = L
Ly = L
dt = 0.01
pbc_layer_depth = 3.5

# Create a box
box = SimulationBox(Lx, Ly)
# Create a VoronoiCells object
target_perimeters = 4.0ones(N)
target_areas = ones(N)
K_P = ones(N)
K_A = ones(N)
active_force_strengths = ones(N)*0.1
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
random_seed = 564574564
Random.seed!(random_seed)
dump_info = DumpInfo(
    save=true,
    filename="Data/dump_$(random_seed).h5",
    when_to_save_array=0:1000:10000000,
)

rng = Random.MersenneTwister(random_seed)
verbose=true
cb = function (x...)
    return nothing 
end
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
arrays.neighborlist.check_tessellation = true
# Run the simulation

Nsteps = 10000
dt = 0.0001
run_simulation!(parameter_struct, arrays, output, Nsteps)
visualize = create_visualization_callback(arrays, round(Int, 0.1/dt), dt, box)

# with a visualization callback
parameter_struct2 = ParameterStruct(N = N, dt = dt, 
    kBT = kBT, frictionconstant = frictionconstant, 
    periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
    dump_info = dump_info, callback = visualize, RNG = rng)

Nsteps = 10000
@time @profview run_simulation!(parameter_struct2, arrays, output, Nsteps) # 4.54 - 4.73 seconds (35.37 M allocations: 2.616 GiB, 8.61% gc time)

