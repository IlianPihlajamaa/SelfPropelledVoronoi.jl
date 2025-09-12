
abstract type Particles end

"""
    VoronoiCells(target_perimeters, target_areas, K_P, K_A, active_force_strengths, rotational_diffusion_constants)

Represents the properties of particles modeled as Voronoi cells.
This struct stores parameters related to the mechanical properties (target area and perimeter, and associated spring constants)
and active behavior (active force strength and rotational diffusion) of each particle/cell.

# Fields
- `target_perimeters::Vector{Float64}`: Target perimeters for each Voronoi cell. Deviations from these values incur an energy penalty.
- `target_areas::Vector{Float64}`: Target areas for each Voronoi cell. Deviations from these values incur an energy penalty.
- `K_P::Vector{Float64}`: Spring constants for perimeter deviations. Determines how strongly a cell resists changes to its perimeter.
- `K_A::Vector{Float64}`: Spring constants for area deviations. Determines how strongly a cell resists changes to its area.
- `active_force_strengths::Vector{Float64}`: Magnitudes of the active force for each particle. This force propels the particle in the direction of its orientation.
- `rotational_diffusion_constants::Vector{Float64}`: Rotational diffusion rates for each particle. This determines how quickly the particle's orientation changes randomly.
"""
struct VoronoiCells <: Particles
    target_perimeters::Vector{Float64}
    target_areas::Vector{Float64}
    K_P::Vector{Float64}
    K_A::Vector{Float64}
    active_force_strengths::Vector{Float64}
    rotational_diffusion_constants::Vector{Float64}
    VoronoiCells(p0s, A0s, KPs, KAs, f0s, Drs) = new(
        p0s,
        A0s,
        KPs,
        KAs,
        f0s,
        Drs
    )
end

"""
    SimulationBox(Lx, Ly)

Defines the simulation domain, typically a rectangular box with periodic boundary conditions.

# Fields
- `box_sizes::SVector{2, Float64}`: A 2D static vector representing the dimensions (Lx, Ly) of the simulation box. `Lx` is the width and `Ly` is the height.
"""
struct SimulationBox
    box_sizes::SVector{2, Float64}
    SimulationBox(Lx::Float64, Ly::Float64) = new(SVector{2, Float64}(Lx, Ly))
end



"""
    NeighborList
A structure to hold the neighbor list and Voronoi edges for each particle.
The `neighbors` field is a vector of vectors, where each inner vector contains the indices of the neighboring particles for each particle.
The `voronoi_vertices` field is a vector of vectors, where each inner vector contains the index of the Voronoi vertices for each particle. The coordinates are accessed by `delaunay_circumcenters` and are stored in counterclockwise order
The `delaunay_circumcenters` field is a vector of vectors, where each inner vector contains the coordinates of the Delaunay circumcenters for each particle. The Delaunay circumcenters are represented as `SVector{2, Float64}`.
"""
mutable struct VoronoiNeighborList
    voronoi_neighbors::Vector{Vector{Int64}}            # List of Voronoi neighbors for each particle
    voronoi_vertex_indices::Vector{Vector{Int64}}     # List of Voronoi vertices for each particle sorted counterclockwise
    voronoi_vertices::Vector{SVector{2,Float64}}                   # List of voronoi_vertices
    positions_with_pbc::Vector{SVector{2, Float64}}               # List of positions with periodic boundary conditions
    position_indices::Vector{Int64}                               # List of indices for the positions with periodic boundary conditions
    delaunay_facet_triplets::Vector{NTuple{3, Int}}               # Stores original particle indices forming Delaunay facets
    N_voronoi_neighbors_pp::Vector{Int64}                   # Number of Voronoi neighbors for each particle
    N_voronoi_vertices_pp::Vector{Int64}                      # Number of Voronoi vertices for each particle
    N_voronoi_vertices::Int64                   # Number of Voronoi vertices
    N_triplets::Int64                     # Number of Delaunay triplets
    N_positions_with_pbc::Int64                     # Number of positions with periodic boundary conditions
    temp_buffer1::Vector{Int64} # Temporary buffer for storing indices during computations
    temp_buffer2::Vector{Int64} # Temporary buffer for storing indices during computations
    temp_buffer3::Vector{Float64} # Temporary buffer for storing distances during computations
    check_tessellation::Bool # Flag which controls whether tessellation is checked or if tessellation 
                # is performed every time step regardless of its validity. For debugging. Set to true for performance.
    VoronoiNeighborList() = new(
        Vector{Vector{Int64}}(),
        Vector{Vector{Int64}}(),
        Vector{SVector{2, Float64}}(),
        Vector{SVector{2, Float64}}(),
        Vector{Int64}(),
        Vector{NTuple{3, Int}}(),
        Vector{Int64}(),
        Vector{Int64}(),
        0,
        0,
        0,
        Vector{Int64}(),
        Vector{Int64}(),
        Vector{Float64}(),
        true
    )
end

"""
    ArrayStruct(N)

A structure that holds various arrays for storing particle properties and the simulation state.
The type parameter `NB` refers to the type of the neighbor list structure used (e.g., `VoronoiNeighborList`).

# Fields
- `positions::Vector{SVector{2, Float64}}`: Current 2D positions of all particles.
- `old_positions::Vector{SVector{2, Float64}}`: Positions of all particles at the previous time step. Used for some integration schemes or state tracking.
- `forces::Vector{SVector{2, Float64}}`: Current 2D forces acting on all particles.
- `old_forces::Vector{SVector{2, Float64}}`: Forces acting on all particles at the previous time step. Used for some integration schemes.
- `orientations::Vector{Float64}`: Current orientations of all particles, represented as an angle.
- `old_orientations::Vector{Float64}`: Orientations of all particles at the previous time step.
- `areas::Vector{Float64}`: Current areas of the Voronoi cell corresponding to each particle.
- `perimeters::Vector{Float64}`: Current perimeters of the Voronoi cell corresponding to each particle.
- `random_forces::Vector{Float64}`: Stores random numbers or forces, often used for implementing stochastic elements like rotational diffusion.
- `neighborlist::NB`: A structure (e.g., `VoronoiNeighborList`) that holds neighbor information for each particle, essential for calculating interactions and Voronoi properties.
"""
struct ArrayStruct{NB}
    positions::Vector{SVector{2, Float64}}
    old_positions::Vector{SVector{2, Float64}}
    forces::Vector{SVector{2, Float64}}
    old_forces::Vector{SVector{2, Float64}}
    orientations::Vector{Float64}
    old_orientations::Vector{Float64}
    areas::Vector{Float64}
    perimeters::Vector{Float64}
    random_forces::Vector{Float64}
    neighborlist::NB
    ArrayStruct(N) = new{VoronoiNeighborList}(
        zeros(SVector{2, Float64}, N),
        zeros(SVector{2, Float64}, N),
        zeros(SVector{2, Float64}, N),
        zeros(SVector{2, Float64}, N),
        zeros(Float64, N),
        zeros(Float64, N),
        zeros(Float64, N),
        zeros(Float64, N),
        zeros(Float64, N),
        VoronoiNeighborList()
    )
end

"""
    Output()

A mutable struct that stores various simulation output quantities. These quantities are updated throughout the simulation and can be used for monitoring progress, logging, or analysis.

# Fields
- `potential_energy::Float64`: The total potential energy of the system at the current state.
- `steps_done::Int64`: The number of simulation time steps that have been completed.
- `N_voronoi_tessellations::Int64`: The total number of Voronoi tessellations that have been performed during the simulation.
"""
mutable struct Output
    potential_energy::Float64
    steps_done::Int64
    N_voronoi_tessellations::Int64
    Output() = new(
        0.0,
        0,
        0
    )
end

"""
    DumpInfo(; save=true, filename="dump_...", when_to_save_array=0:1000:1000000, save_r=true, save_F=false, save_Epot=false)

A mutable struct that holds parameters for controlling the dumping (saving) of simulation data to a file.
The type parameter `A` corresponds to the type of `when_to_save_array`, which is typically an array or range of integers.

# Fields
- `save::Bool`: A boolean flag to enable (`true`) or disable (`false`) the saving of simulation data. Defaults to `true`.
- `filename::String`: The name of the file where simulation data will be saved. Defaults to a randomly generated HDF5 filename (e.g., "dump_12345.h5").
- `when_to_save_array::A`: An array or range specifying the simulation steps at which data should be saved. `A` is the type of this field, typically `AbstractArray`. Defaults to saving every 1000 steps up to 1,000,000.
- `save_r::Bool`: A boolean flag indicating whether to save particle positions (`true`) or not (`false`). Defaults to `true`.
- `save_u`::Bool`: A boolean flag indicating whether to save particle orientations (`true`) or not (`false`). Defaults to `true`.
- `save_F::Bool`: A boolean flag indicating whether to save particle forces (`true`) or not (`false`). Defaults to `false`.
- `save_Epot::Bool`: A boolean flag indicating whether to save the total potential energy of the system (`true`) or not (`false`). Defaults to `false`.
- `save_areas::Bool`: A boolean flag indicating whether to save the areas of Voronoi cells (`true`) or not (`false`). Defaults to `true`.
- `save_perimeters::Bool`: A boolean flag indicating whether to save the perimeters of Voronoi cells (`true`) or not (`false`). Defaults to `true`.
"""
mutable struct DumpInfo{A<:AbstractArray}
    save::Bool
    filename::String
    when_to_save_array::A
    save_r::Bool
    save_u::Bool
    save_F::Bool
    save_Epot::Bool
    save_areas::Bool
    save_perimeters::Bool
    DumpInfo(;
        save::Bool=true,
        filename::String="dump_$(rand(Int)).h5",
        when_to_save_array=0:1000:1000000,
        save_r::Bool=true,
        save_u::Bool=true,
        save_F::Bool=false,
        save_Epot::Bool=false,
        save_areas::Bool=true,
        save_perimeters::Bool=true
    ) = new{typeof(when_to_save_array)}(
        save,
        filename,
        when_to_save_array,
        save_r,
        save_u,
        save_F,
        save_Epot,
        save_areas,
        save_perimeters
    )
end

mutable struct TrajectoryData
    positions_trajectory::Vector{Vector{SVector{2, Float64}}}
    orientations_trajectory::Vector{Vector{Float64}}
    forces_trajectory::Vector{Vector{SVector{2, Float64}}}
    potential_energy_trajectory::Vector{Float64}
    areas_trajectory::Vector{Vector{Float64}}
    perimeters_trajectory::Vector{Vector{Float64}}
    steps_saved::Vector{Int64}
    TrajectoryData() = new(
        Vector{Vector{SVector{2, Float64}}}(),
        Vector{Vector{Float64}}(),
        Vector{Vector{SVector{2, Float64}}}(),
        Vector{Float64}(),
        Vector{Vector{Float64}}(),
        Vector{Vector{Float64}}(),
        Vector{Int64}()
    )
end


"""
    ParameterStruct(N =100, dt = 0.01, N_steps = 10000, kBT = 1.0, frictionconstant = 1.0,
                    periodic_boundary_layer_depth = 1.0, verbose = false, box = SimulationBox(10.0, 10.0),
                    particles = VoronoiCells(zeros(Float64, 100), zeros(Float64, 100), zeros(Float64, 100),
                                              zeros(Float64, 100), ones(Float64, 100), ones(Float64, 100)),
                    dump_info = DumpInfo(), callback = x -> nothing,
                    RNG = Random.MersenneTwister(1234))

A struct that holds all essential parameters for running a simulation.
This includes physical parameters, numerical parameters, and settings for I/O and control.

Type parameters:
- `P`: The type of the particle properties structure (e.g., `VoronoiCells`), which must be a subtype of `Particles`.
- `CB`: The type of the callback function.

# Fields
- `N::Int`: The total number of particles in the simulation.
- `dt::Float64`: The time step size for the integration algorithm.
- `N_steps::Int`: The total number of simulation steps to perform.
- `kBT::Float64`: The thermal energy, given by the product of Boltzmann's constant (kB) and temperature (T).
- `frictionconstant::Float64`: The friction coefficient, determining the damping of particle motion.
- `periodic_boundary_layer_depth::Float64`: The depth of the layer around the simulation box used for periodic boundary condition checks and neighbor finding.
- `verbose::Bool`: A boolean flag to enable (`true`) or disable (`false`) verbose output during the simulation.
- `box::SimulationBox`: A `SimulationBox` struct defining the dimensions and properties of the simulation domain.
- `particles::P`: A struct (of type `P`) holding the specific properties of the particles (e.g., `VoronoiCells` which includes target areas, perimeters, etc.).
- `dump_info::DumpInfo`: A `DumpInfo` struct containing parameters for controlling how simulation data is saved to files.
- `callback::CB`: A callback function (of type `CB`) that can be executed at specified intervals during the simulation (e.g., for custom analysis or logging).
- `RNG::Random.MersenneTwister`: An instance of a random number generator (specifically `MersenneTwister`) used for stochastic processes in the simulation.
- `seed::Int64`: The seed value used to initialize the random number generator, stored for reproducibility.
"""
mutable struct ParameterStruct{P<:Particles, CB}
    N::Int
    dt::Float64
    kBT::Float64
    frictionconstant::Float64
    periodic_boundary_layer_depth::Float64
    verbose::Bool
    box::SimulationBox
    particles::P
    dump_info::DumpInfo
    callback::CB
    RNG::Random.MersenneTwister
    seed::Int64
    function ParameterStruct(;
        N::Int=100,
        dt::Float64=0.01,
        kBT::Float64=1.0,
        frictionconstant::Float64=1.0,
        periodic_boundary_layer_depth::Float64=2.5,
        verbose::Bool=false,
        box::SimulationBox=SimulationBox(10.0, 10.0),
        particles=VoronoiCells(
            zeros(Float64, 100),  # target_perimeters
            zeros(Float64, 100),  # target_areas
            zeros(Float64, 100),  # K_P
            zeros(Float64, 100),  # K_A
            ones(Float64, 100),   # active_force_strengths
            ones(Float64, 100)    # rotational_diffusion_constants
        ),
        dump_info::DumpInfo=DumpInfo(),
        callback=x->nothing,  # Default callback does nothing
        RNG::Random.MersenneTwister=Random.MersenneTwister(1234),
    ) 
        return new{typeof(particles), typeof(callback)}(
            N,
            dt,
            kBT,
            frictionconstant,
            periodic_boundary_layer_depth,
            verbose,
            box,
            particles,
            dump_info,
            callback,
            RNG, 
            Int(RNG.seed)
        )
    end
        
end

