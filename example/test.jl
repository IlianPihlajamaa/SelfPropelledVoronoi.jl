import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, CairoMakie, StaticArrays, Random, Quickhull

N = 1000
Lx = 10.0
Ly = 10.0
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
voronoi_cells = VoronoiCells(
    target_perimeters,
    target_areas,
    K_P,
    K_A,
    active_force_strengths
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
    if !(step % 100 == 0)
        return 
    end
    # visualize the initial configuration
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1, 1], title="configuration at step $(output.steps_done)", limits=(0.0, Lx, 0.0, Ly))
    scatter!(ax, arrays.positions, markersize=15)
    tri = delaunay(arrays.positions)
    vor_edges = voronoi_edge_points(tri)
    linesegments!(ax, vor_edges, color=:red)
    display(fig)
end


callback = visualize 
parameter_struct = ParameterStruct(
    N,
    dt,
    Nsteps,
    kBT,
    frictionconstant,
    box,
    voronoi_cells,
    dump_info,
    callback,
    MersenneTwister(random_seed)
)

# Create an ArrayStruct object
arrays = ArrayStruct(N, 16)
arrays.positions .= [rand(SVector{2, Float64}) .* box.box_sizes for _ in 1:N]
arrays.orientations .= 2π*rand(Float64, N) 


# Create an Output object
output = Output()





