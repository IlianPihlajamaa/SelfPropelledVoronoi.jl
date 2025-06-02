# This file is for my personal testing and development of the SelfPropelledVoronoi.jl package.

import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, StaticArrays, Random
using Statistics, BenchmarkTools

N = 400
rho = 1.0
L = sqrt(N/rho)
Lx = L
Ly = L
dt = 0.1
pbc_layer_depth = 2.5

# Create a box
box = SimulationBox(Lx, Ly)
# Create a VoronoiCells object
target_perimeters = 3.85*ones(N)
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
    save=false,
    filename="Data/dump_$(random_seed).h5",
    when_to_save_array=0:1000:1000000,
)

rng = Random.MersenneTwister(random_seed)
verbose=true
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

Nsteps = 500 ÷ dt
run_simulation!(parameter_struct, arrays, output, Nsteps)
visualize = SelfPropelledVoronoi.create_visualization_callback(arrays, 25, dt, box)

# with a visualization callback
parameter_struct2 = ParameterStruct(N = N, dt = dt, 
    kBT = kBT, frictionconstant = frictionconstant, 
    periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
    dump_info = dump_info, callback = visualize, RNG = rng)

Nsteps = 1000 ÷ dt
run_simulation!(parameter_struct2, arrays, output, Nsteps)