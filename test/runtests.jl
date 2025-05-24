using SelfPropelledVoronoi
using Test

@testset "SelfPropelledVoronoi.jl" begin

    @test SelfPropelledVoronoi.test_force(1) == 2

    include("test_forces.jl")
end
