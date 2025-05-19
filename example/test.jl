import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, CairoMakie, StaticArrays, Random, Quickhull

N = 50
rho = 0.7
L = sqrt(N/rho)
Lx = L
Ly = L
dt = 0.01
Nsteps = 1000

# Create a box
box = SimulationBox(Lx, Ly)
# Create a VoronoiCells object
target_perimeters = ones(N)
target_areas = ones(N)
K_P = ones(N)
K_A = ones(N)
active_force_strengths = ones(N)
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
random_seed = rand(UInt32)
Random.seed!(random_seed)
dump_info = DumpInfo(
    save=true,
    filename="dump_$(random_seed).h5",
    when_to_save_array=0:1000:1000000,
    save_r=true,
    save_F=false,
    save_Epot=false
)

function visualize(parameters, arrays, output)
    if !(output.steps_done % 100 == 0)
        return 
    end
    Lx, Ly = parameters.box.box_sizes
    # visualize the initial configuration
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1, 1], title="configuration at step $(output.steps_done)", limits=(0.0, Lx, 0.0, Ly))
    scatter!(ax, arrays.positions, markersize=25)

    D_pbc = 2.5
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc
    positions_with_pbc = copy(arrays.positions)
    for i in 1:N
        x, y = arrays.positions[i]
        if x < D_pbc
            push!(positions_with_pbc, SVector(x + Lx, y))
        end
        if x > Lx - D_pbc
            push!(positions_with_pbc, SVector(x - Lx, y))
        end
        if y < D_pbc
            push!(positions_with_pbc, SVector(x, y + Ly))
        end
        if y > Ly - D_pbc
            push!(positions_with_pbc, SVector(x, y - Ly))
        end
    end



    @time tri = Quickhull.delaunay(positions_with_pbc)
    @time vor_vertices = Quickhull.voronoi_centers(tri)
    @time vor_edges = Quickhull.voronoi_edges(tri)
    @time vor_edge_points = [(vor_vertices[i], vor_vertices[j]) for (i, j) in vor_edges]

    scatter!(ax, vor_vertices, markersize=15, color=:green)
    linesegments!(vor_edge_points, color=:black)

    for facet in facets(tri)
        i = facet[1]
        j = facet[2]
        k = facet[3]

        if i < j
            linesegments!(ax, [
                (positions_with_pbc[i], positions_with_pbc[j]),
            ], color=:red, linestyle=:dash)

        else
            linesegments!(ax, [
                (positions_with_pbc[j], positions_with_pbc[i]),
            ], color=:red, linestyle=:dash)
        end

        if j < k
            linesegments!(ax, [
                (positions_with_pbc[j], positions_with_pbc[k]),
            ], color=:red, linestyle=:dash)

        else
            linesegments!(ax, [
                (positions_with_pbc[k], positions_with_pbc[j]),
            ], color=:red, linestyle=:dash)
        end

        if k < i
            linesegments!(ax, [
                (positions_with_pbc[k], positions_with_pbc[i]),
            ], color=:red, linestyle=:dash)

        else
            linesegments!(ax, [
                (positions_with_pbc[i], positions_with_pbc[k]),
            ], color=:red, linestyle=:dash)
        end

        # color every facet with a different color
  
        poly!(ax, [positions_with_pbc[i], positions_with_pbc[j], positions_with_pbc[k]], alpha=0.2)

    end

    display(fig)
end

verbose=true
callback = visualize 
parameter_struct = ParameterStruct(
    N,
    dt,
    Nsteps,
    kBT,
    frictionconstant,
    verbose,
    box,
    voronoi_cells,
    dump_info,
    callback,
    MersenneTwister(random_seed)
)

# Create an ArrayStruct object
arrays = ArrayStruct(N, 16)

# arrays.positions .= [rand(SVector{2, Float64}) .* box.box_sizes for _ in 1:N]
#put particles on cubic lattice
x = LinRange(0, box.box_sizes[1], ceil(Int, sqrt(N)))
y = LinRange(0, box.box_sizes[2], ceil(Int, sqrt(N)))
x = repeat(x, inner=(ceil(Int, sqrt(N)), 1))
y = repeat(y, inner=(1, ceil(Int, sqrt(N))))
x = x[1:N]
y = y[1:N]
arrays.positions .= SVector.(x, y)


arrays.orientations .= 2π*rand(Float64, N) 


# Create an Output object
output = Output()

# Run the simulation
run_simulation!(parameter_struct, arrays, output)



