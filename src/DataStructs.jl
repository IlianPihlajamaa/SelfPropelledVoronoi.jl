abstract type Particles end

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

struct SimulationBox
    box_sizes::SVector{2, Float64}
    SimulationBox(Lx::Float64, Ly::Float64) = new(SVector{2, Float64}(Lx, Ly))
end



"""
    NeighborList
A structure to hold the neighbor list and Voronoi edges for each particle.
The `neighbors` field is a vector of vectors, where each inner vector contains the indices of the neighboring particles for each particle.
The `voronoi_vertices` field is a vector of vectors, where each inner vector contains the index of the Voronoi vertices for each particle. The coordinates are accessed by `delauney_circumcenters` and are stored in counterclockwise order
The `delauney_circumcenters` field is a vector of vectors, where each inner vector contains the coordinates of the Delauney circumcenters for each particle. The Delauney circumcenters are represented as `SVector{2, Float64}`.
"""
mutable struct VoronoiNeighborList
    voronoi_neighbors::Vector{Vector{Int64}}            # List of Voronoi neighbors for each particle
    voronoi_vertex_indices::Vector{Vector{Int}}     # List of Voronoi vertices for each particle sorted counterclockwise
    voronoi_vertex_positions_per_particle::Vector{Vector{SVector{2, Float64}}}  # List of Voronoi vertex positions for each particle
    voronoi_vertices::Vector{SVector{2,Float64}}                   # List of voronoi_vertices
    positions_with_pbc::Vector{SVector{2, Float64}}               # List of positions with periodic boundary conditions
    position_indices::Vector{Int64}                               # List of indices for the positions with periodic boundary conditions
    VoronoiNeighborList(N) = new(
        [Vector{Int64}() for _ in 1:N],
        [Vector{Int64}() for _ in 1:N],
        [SVector{2, Float64}[] for _ in 1:N],
        [zero(SVector{2, Float64}) for _ in 1:N],
        SVector{2, Float64}[],
        Int64[]
    )
end


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
        VoronoiNeighborList(N)
    )
end





mutable struct Output
    potential_energy::Float64                   # Total potential energy of all particles
    steps_done::Int64                           # Number of MC steps done
    N_voronoi_tesselations::Int64               # Number of Voronoi tesselations done
    Output() = new(
        0.0,
        0,
        0
    )
end


mutable struct DumpInfo{A<:AbstractArray}
    save::Bool
    filename::String
    when_to_save_array::A
    save_r::Bool
    save_F::Bool
    save_Epot::Bool
    DumpInfo(;
        save::Bool=true,
        filename::String="dump_$(rand(Int)).h5",
        when_to_save_array=0:1000:1000000,
        save_r::Bool=true,
        save_F::Bool=false,
        save_Epot::Bool=false
    ) = new{typeof(when_to_save_array)}(
        save,
        filename,
        when_to_save_array,
        save_r,
        save_F,
        save_Epot
    )
end



struct ParameterStruct{P<:Particles, CB}
    N::Int
    dt::Float64
    N_steps::Int
    kBT::Float64
    frictionconstant::Float64
    periodic_boundary_layer_depth::Float64
    verbose::Bool
    box::SimulationBox
    particles::P
    dump_info::DumpInfo
    callback::CB
    RNG::Random.MersenneTwister
end
