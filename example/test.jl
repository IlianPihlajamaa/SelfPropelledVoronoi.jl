import Pkg; Pkg.activate("example")

using Revise
using SelfPropelledVoronoi, CairoMakie, StaticArrays, Random, ColorSchemes
using Statistics


N = 400
rho = 1.0
L = sqrt(N/rho)
Lx = L
Ly = L
dt = 0.01

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
    save=true,
    filename="Data/dump_$(random_seed).h5",
    when_to_save_array=0:1000:1000000,
)


energy_list = Float64[]
area_std_list = Float64[]
mean_perimeter_list = Float64[]
mean_squared_displacement_list = Float64[]
t_array = Float64[]
x0 = rand(Float64, N) .* box.box_sizes[1]
y0 = rand(Float64, N) .* box.box_sizes[2]
displacement_array = [SVector{2, Float64}(0.0, 0.0) for _ in 1:N]
previous_positions = SVector{2, Float64}.(x0, y0)


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

function visualize(parameters, arrays, output)

    # apply periodic boundary conditions

    push!(energy_list, SelfPropelledVoronoi.compute_energy(parameters, arrays, output))
    push!(mean_perimeter_list, mean(arrays.perimeters))
    push!(area_std_list, std(arrays.areas))
    push!(t_array, output.steps_done * parameters.dt)
    if !(length(t_array) == 1) # no movement yet
        dr = arrays.positions .- previous_positions
        apply_periodic_boundary_conditions!(dr, parameters.box.box_sizes) 
        displacement_array .+= dr
        push!(mean_squared_displacement_list, compute_msd(displacement_array))
    end
    previous_positions .= arrays.positions

    if !(output.steps_done % 1000 == 0) 
        return
    end


    Lx, Ly = parameters.box.box_sizes
    # visualize the initial configuration
    fig = Figure(size=(1000,500))
    

    ax1 = Axis(fig[1:2, 1], title="configuration", limits=(-pbc_layer_depth, Lx+pbc_layer_depth, -pbc_layer_depth, Ly+pbc_layer_depth))

    # # draw the voronoi cells
    scatter!(ax1, arrays.neighborlist.positions_with_pbc, markersize=12, color =:blue)


    # draw voronoi edges
    linesegs2 = Tuple{SVector{2, Float64}, SVector{2, Float64}}[]
    for particle_i in 1:length(arrays.neighborlist.positions_with_pbc)
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
    linesegments!(ax1, linesegs2, color=:black)
    
    # draw original box
    linesegments!(ax1, [
        (0, 0), (Lx, 0),
        (Lx, 0), (Lx, Ly),
        (Lx, Ly), (0, Ly),
        (0, Ly), (0, 0)
    ], color=:black, linewidth=4)

    ax2 = Axis(fig[1, 2], title="energy", ylabel="energy")
    lines!(ax2, dt*(100:length(energy_list)), energy_list[100:end], color=:blue)


    ax4 = Axis(fig[2, 2], title="msd", ylabel="msd", xlabel="t", yscale=log10, xscale=log10)
    lines!(ax4, t_array[2:end] .- t_array[1], (mean_squared_displacement_list), color=:blue, label="msd")
    display(fig)
end
rng = Random.MersenneTwister(random_seed)
verbose=true
cb(x...) = nothing 
parameter_struct = ParameterStruct(N = N, dt = dt, 
    kBT = kBT, frictionconstant = frictionconstant, 
    periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
    dump_info = dump_info, callback = cb, RNG = rng)


# Create an ArrayStruct object
arrays = ArrayStruct(N)

# arrays.positions .= [rand(SVector{2, Float64}) .* box.box_sizes for _ in 1:N]
#put particles on cubic lattice

arrays.positions .= SVector.(x0, y0)
arrays.orientations .= 2π*rand(Float64, N) 


# Create an Output object
output = Output()

# Run the simulation
Nsteps = 5000
run_simulation!(parameter_struct, arrays, output, Nsteps)
previous_positions .= arrays.positions
parameter_struct2 = ParameterStruct(N = N, dt = dt, 
    kBT = kBT, frictionconstant = frictionconstant, 
    periodic_boundary_layer_depth = pbc_layer_depth, verbose = verbose, box = box, particles= voronoi_cells,
    dump_info = dump_info, callback = visualize, RNG = rng)

Nsteps = 100000
run_simulation!(parameter_struct2, arrays, output, Nsteps)