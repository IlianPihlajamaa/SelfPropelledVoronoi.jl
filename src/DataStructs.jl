struct ParameterStruct{P<:Particles, CB}
    N::Int
    kBT::Float64
    frictionconstant::Float64
    random_seed::Int
    box::Box
    particles::P
    dump_info::DumpInfo
    callback::CB
end

abstract type Particles end

struct VoronoiCells <: Particles
    target_perimeters::Vector{Float64}
    target_areas::Vector{Float64}
    K_P::Vector{Float64}
    K_A::Vector{Float64}
    active_force_strengths::Vector{Float64}
    VoronoiCells(p0s, A0s, KPs, KAs, f0s) = new(
        p0s,
        A0s,
        KPs,
        KAs,
        f0s
    )
end

struct Box
    box_sizes::SVector{2, Float64}
    Box(Lx::Float64, Ly::Float64) = new(SVector{2, Float64}(Lx, Ly))
end

struct ArrayStruct
    positions::Vector{SVector{2, Float64}}
    forces::Vector{SVector{2, Float64}}
    orientations::Vector{Float64}
    areas::Vector{Float64}
    perimeters::Vector{Float64}
    ArrayStruct(N) = new(
        zeros(SVector{2, Float64}, N),
        zeros(SVector{2, Float64}, N),
        zeros(Float64, N),
        zeros(Float64, N),
        zeros(Float64, N)
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


