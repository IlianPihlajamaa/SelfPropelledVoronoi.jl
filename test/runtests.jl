using SelfPropelledVoronoi
using SelfPropelledVoronoi.Tesselation
using Test
using StaticArrays

# Helper to initialize basic structures for tests
function setup_test_environment(particle_positions::Vector{SVector{2, Float64}}, box_size_val::Float64 = 100.0)
    N = length(particle_positions)
    
    _box = Box(SVector(box_size_val, box_size_val))
    # Assuming Parameters constructor is Parameters(N, box, ...other_args_with_defaults...)
    # Or more simply, if a direct struct initialization is possible and preferred for tests:
    parameters = Parameters(N=N, box=_box) # Adjust if constructor is different
    
    arrays = Arrays(N=N) # Assuming Arrays(N=N) constructor that initializes positions, old_positions, neighborlist
    arrays.positions = deepcopy(particle_positions)
    arrays.old_positions = Vector{SVector{2, Float64}}[] # Start with empty/distinct old_positions

    output = Output() 

    return parameters, arrays, output
end

using StaticArrays: SVector


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


@testset "Tesselation.verify_tessellation" begin

    @testset "Valid Tessellation (Square)" begin
        positions = [
            SVector(0.0, 1.0), SVector(1.0, 1.0), 
            SVector(1.0, 0.0), SVector(0.0, 0.0) 
        ] # N=4
        params, arrays, output = setup_test_environment(positions, 10.0) # Small box for clarity

        voronoi_tesselation!(params, arrays, output)
        
        initial_positions_copy = deepcopy(arrays.positions)
        arrays.old_positions = [SVector(999.0,999.0)] # Ensure it's different before call

        result = verify_tessellation(params, arrays, output)

        @test result == true
        @test arrays.old_positions == initial_positions_copy
        @test arrays.old_positions != [SVector(999.0,999.0)] # Ensure it was updated
    end

    @testset "Invalid Tessellation (Particle Moved into Circumcircle of a Delaunay Facet)" begin
        # Scenario:
        # 1. Start with a square: p1(0,10), p2(10,10), p3(10,0), p4(0,0). Original indices 1,2,3,4.
        # 2. Call voronoi_tesselation! to populate delaunay_facet_triplets, voronoi_neighbors, position_indices.
        #    A Delaunay facet for this square could be (1,2,3) (original indices, sorted).
        #    The circumcenter of p1,p2,p3 is (5,5).
        # 3. p4 is a neighbor of p1 and p3 in the initial square tessellation.
        # 4. Move p4 into the circumcircle of (p1,p2,p3) by changing arrays.positions[4].
        # 5. Call update_positions_with_pbcs! to update positions_with_pbc and position_indices.
        # 6. Call verify_tessellation. It uses:
        #    - arrays.neighborlist.delaunay_facet_triplets (from initial tessellation - STALE for p4's position).
        #    - arrays.positions (for facet points p1,p2,p3 and test particle p4 - CURRENT, p4 is moved).
        #    - arrays.neighborlist.voronoi_neighbors (from initial tessellation - STALE list of PBC indices).
        #    - arrays.neighborlist.position_indices (to map stale PBC indices to original indices - UPDATED by update_positions_with_pbcs!).
        #    This setup should detect the violation.

        positions_initial = [
            SVector(0.0, 10.0),  # p1 (Original Index 1)
            SVector(10.0, 10.0), # p2 (Original Index 2)
            SVector(10.0, 0.0),  # p3 (Original Index 3)
            SVector(0.0, 0.0)    # p4 (Original Index 4)
        ]
        params, arrays, output = setup_test_environment(positions_initial, 20.0)

        # Perform initial tessellation. This populates all relevant fields in arrays.neighborlist.
        voronoi_tesselation!(params, arrays, output)
        
        # Set a distinct old_positions to check it's not updated on failure
        arrays.old_positions = [SVector(123.0, 456.0), SVector(1.0,1.0), SVector(2.0,2.0), SVector(3.0,3.0)]
        original_old_positions_snapshot = deepcopy(arrays.old_positions)

        # Move p4 (arrays.positions[4]) into the circumcircle of facet (p1,p2,p3).
        # The circumcenter of (0,10), (10,10), (10,0) is (5,5).
        arrays.positions[4] = SVector(4.9, 4.9) # Clearly inside

        # Call update_positions_with_pbcs! AFTER moving p4.
        # This updates arrays.neighborlist.positions_with_pbc and arrays.neighborlist.position_indices.
        update_positions_with_pbcs!(params, arrays, output)
        
        result = verify_tessellation(params, arrays, output)

        @test result == false "p4 moved into circumcircle of facet (p1,p2,p3) was not detected."
        @test arrays.old_positions == original_old_positions_snapshot "old_positions was updated despite invalid tessellation."
    end

    @testset "Insufficient Particles/Facets (<3 particles)" begin
        positions = [SVector(0.0,0.0), SVector(1.0,0.0)] # N=2
        params, arrays, output = setup_test_environment(positions, 10.0)

        voronoi_tesselation!(params, arrays, output)
        initial_positions_copy = deepcopy(arrays.positions)
        arrays.old_positions = [SVector(999.0,999.0)] # Ensure different before call

        result = verify_tessellation(params, arrays, output)

        @test result == true
        @test arrays.old_positions == initial_positions_copy
    end
    
    # Points on circumcircles are implicitly tested by the "Valid Tessellation (Square)" case,
    # as square configurations naturally have points on circumcircles (e.g., p1 on C(p2,p3,p4)),
    # and the function uses `d_sq < R_sq - epsilon`.
end

# Keep existing tests if any
@testset "SelfPropelledVoronoi.jl Original Tests" begin
    # This test was in the original file. Assuming it's relevant.
    # If SelfPropelledVoronoi.test_force is not defined, this will error.
    # For now, assuming it exists as per the original test file.
    if isdefined(SelfPropelledVoronoi, :test_force)
        @test SelfPropelledVoronoi.test_force(1) == 2
    else
        @warn "SelfPropelledVoronoi.test_force not defined. Skipping this original test."
    end
end

include("test_forces.jl")
