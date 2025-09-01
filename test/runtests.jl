using SelfPropelledVoronoi
using Test
using HDF5
using StaticArrays
import Random

# Helper function to set up mock objects for testing



@testset "apply_periodic_boundary_conditions" begin
    @testset "2D" begin
        @test SelfPropelledVoronoi.apply_periodic_boundary_conditions(SVector(0.5, 0.5), SVector(1.0, 1.0)) ≈ SVector(0.5, 0.5)
        @test SelfPropelledVoronoi.apply_periodic_boundary_conditions(SVector(1.5, -0.5), SVector(1.0, 1.0)) ≈ SVector(0.5, 0.5)
        @test SelfPropelledVoronoi.apply_periodic_boundary_conditions(SVector(1.0, 0.0), SVector(1.0, 1.0)) ≈ SVector(0.0, 0.0)
    end
end

@testset "compute_pair_distance_vector" begin

    @testset "2D" begin
        @test SelfPropelledVoronoi.compute_pair_distance_vector(SVector(0.2, 0.2), SVector(0.4, 0.4), SVector(1.0, 1.0)) ≈ SVector(0.2, 0.2)
        @test SelfPropelledVoronoi.compute_pair_distance_vector(SVector(0.8, 0.2), SVector(0.2, 0.4), SVector(1.0, 1.0)) ≈ SVector(0.4, 0.2)
        @test SelfPropelledVoronoi.compute_pair_distance_vector(SVector(0.8, 0.8), SVector(0.2, 0.2), SVector(1.0, 1.0)) ≈ SVector(0.4, 0.4)
    end
end



include("test_forces.jl")
include("test_tessellation.jl")
include("test_loading_and_saving.jl")
              
              
              
              