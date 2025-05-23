import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, CairoMakie, StaticArrays, Random, Quickhull, ColorSchemes
using Statistics
N = 500
rho = 0.9
L = sqrt(N/rho)
Lx = L
Ly = L
dt = 0.01
Nsteps = 1000

pbc_layer_depth = 2.5


# Create a box
box = SimulationBox(Lx, Ly)
# Create a VoronoiCells object
target_perimeters = 4.1*ones(N)
target_areas = ones(N)
K_P = 10*ones(N)
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
random_seed = 564574564
Random.seed!(random_seed)
dump_info = DumpInfo(
    save=true,
    filename="dump_$(random_seed).h5",
    when_to_save_array=0:1000:1000000,
    save_r=true,
    save_F=false,
    save_Epot=false
)


energy_list = Float64[]
area_std_list = Float64[]
mean_perimeter_list = Float64[]

function visualize(parameters, arrays, output)
    push!(energy_list, SelfPropelledVoronoi.compute_energy(parameters, arrays, output))
    push!(mean_perimeter_list, mean(arrays.perimeters))
    push!(area_std_list, std(arrays.areas))
    if !(output.steps_done % 100 == 0)
        return 
    end
    SelfPropelledVoronoi.voronoi_tesselation!(parameters, arrays, output)
    Lx, Ly = parameters.box.box_sizes
    # visualize the initial configuration
    fig = Figure(size=(2000,1000))
    

    ax1 = Axis(fig[1:2, 1], title="from neighborlist", limits=(-pbc_layer_depth, Lx+pbc_layer_depth, -pbc_layer_depth, Ly+pbc_layer_depth))
    #draw original box 
   
    # draw the voronoi cells
    n_particles_pbc = length(arrays.neighborlist.positions_with_pbc)
    scatter!(ax1, arrays.neighborlist.positions_with_pbc, markersize=25, color =:blue)

    #draw voronoi vertices 
    scatter!(ax1, arrays.neighborlist.voronoi_vertices, markersize=15, color=:green)
    # draw delaunay edges (voronoi neighbors)
    for particle_i in 1:length(arrays.positions)
        for particle_j in arrays.neighborlist.voronoi_neighbors[particle_i]
            if particle_i > particle_j
                continue
            end
            posi = arrays.positions[particle_i]
            posj = arrays.positions[particle_j]

            dr = posj - posi
            # apply periodic boundary conditions
            dr = dr - round.(dr ./ parameters.box.box_sizes) .* parameters.box.box_sizes
            posj2 = posi + dr
            linesegments!(ax1, [
                (posi, posj2),
            ], color=:red, linestyle=:dash)

            posi2 = posj - dr
            linesegments!(ax1, [
                (posi2, posj),
            ], color=:red, linestyle=:dash)
        end
    end
    # draw voronoi edges
    for particle_i in 1:length(arrays.positions)
        indices = arrays.neighborlist.voronoi_vertex_indices[particle_i]
        vor_positions = arrays.neighborlist.voronoi_vertex_positions_per_particle[particle_i]

        for i in 1:length(indices)
            j = i % length(indices) + 1

            
            posi = vor_positions[i]
            posj = vor_positions[j]

            # only if both are in the original box
            # if !( 0 < posi[1] < Lx && 0 < posi[2] < Ly && 0 < posj[1] < Lx && 0 < posj[2] < Ly)
            #     continue
            # end


            # apply periodic boundary conditions
            linesegments!(ax1, [
                (posi, posj),
            ], color=:black)

        end
    end

    linesegments!(ax1, [
        (0, 0), (Lx, 0),
        (Lx, 0), (Lx, Ly),
        (Lx, Ly), (0, Ly),
        (0, Ly), (0, 0)
    ], color=:black, linewidth=4)

    ax2 = Axis(fig[1, 2], title="energy", ylabel="energy")
    lines!(ax2, 1:length(energy_list), energy_list, color=:blue)

    ax3 = Axis(fig[2, 2], title="perimeter")
    band!(ax3, 1:length(area_std_list), mean_perimeter_list - area_std_list, mean_perimeter_list + area_std_list, color=(:blue,0.3), label="std perimeter")
    lines!(ax3, 1:length(mean_perimeter_list), mean_perimeter_list, color=:blue, label="mean perimeter")
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
    pbc_layer_depth,
    verbose,
    box,
    voronoi_cells,
    dump_info,
    callback,
    MersenneTwister(random_seed)
)

# Create an ArrayStruct object
arrays = ArrayStruct(N)

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



