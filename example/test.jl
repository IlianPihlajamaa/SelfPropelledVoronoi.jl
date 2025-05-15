import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, CairoMakie, StaticArrays, Random

N = 100
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
random_seed = 1234
Random.seed!(random_seed)
dump_info = DumpInfo(
    save=true,
    filename="dump_$(random_seed).h5",
    when_to_save_array=0:1000:1000000,
    save_r=true,
    save_F=false,
    save_Epot=false
)
callback = (x) -> nothing # Placeholder for a callback function
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
arrays = ArrayStruct(N)
arrays.positions .= [rand(SVector{2, Float64}) .* box.box_sizes for _ in 1:N]
arrays.orientations .= 2π*rand(Float64, N) 


# Create an Output object
output = Output()


# visualize the initial configuration
fig = Figure()
ax = Axis(fig[1, 1], title="Initial Configuration")
scatter!(ax, arrays.positions, markersize=5)
display(fig)


