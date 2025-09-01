# Helper to initialize basic structures for tests
function setup_test_environment(particle_positions::Vector{SVector{2, Float64}}, box_size_val::Float64 = 100.0)
    N = length(particle_positions)
    
    _box = SimulationBox(box_size_val, box_size_val)

    parameters = ParameterStruct(N=N, box=_box) 
    
    arrays = ArrayStruct(N) 
    arrays.positions .= deepcopy(particle_positions)

    output = Output() 

    return parameters, arrays, output
end




@testset "Tessellation.verify_tessellation" begin

    @testset "Valid Tessellation (Square)" begin

        Lx, Ly = 10.0, 10.0
        N = 100
        x = rand(Float64, N) .* Lx
        y = rand(Float64, N) .* Ly
        positions_initial = [SVector(x[i], y[i]) for i in 1:N]
        params, arrays, output = setup_test_environment(positions_initial, 20.0)
        arrays.neighborlist.check_tessellation = true
        # Perform initial tessellation. This populates all relevant fields in arrays.neighborlist.
        SelfPropelledVoronoi.voronoi_tessellation!(params, arrays, output)

        positions_new = deepcopy(arrays.positions)
        arrays.positions .= positions_new
        # Call verify_tessellation with the current positions and the stale tessellation data.
        result = SelfPropelledVoronoi.verify_tessellation(params, arrays, output)
        @test result == true
    end

    @testset "Invalid Tessellation (Particle Moved into Circumcircle of a Delaunay Facet)" begin
        # Scenario:
        # 1. Start with a square: p1(0,10), p2(10,10), p3(10,0), p4(0,0). Original indices 1,2,3,4.
        # 2. Call voronoi_tessellation! to populate delaunay_facet_triplets, voronoi_neighbors, position_indices.
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
        arrays.neighborlist.check_tessellation = true

        # Perform initial tessellation. This populates all relevant fields in arrays.neighborlist.
        SelfPropelledVoronoi.voronoi_tessellation!(params, arrays, output)
        
        # Set a distinct old_positions to check it's not updated on failure
        arrays.positions .= [SVector(123.0, 456.0), SVector(1.0,1.0), SVector(2.0,2.0), SVector(3.0,3.0)]
        original_old_positions_snapshot = deepcopy(arrays.positions)

        # Move p4 (arrays.positions[4]) into the circumcircle of facet (p1,p2,p3).
        # The circumcenter of (0,10), (10,10), (10,0) is (5,5).
        arrays.positions[4] = SVector(4.9, 4.9) # Clearly inside
    
        result = SelfPropelledVoronoi.verify_tessellation(params, arrays, output)

        @test result == false 
    end



end