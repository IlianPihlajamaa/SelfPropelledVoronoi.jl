using SelfPropelledVoronoi
using Test
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

@testset "Cell Property Updates" begin
    # Removed: using SelfPropelledVoronoi: ParameterStruct, ArrayStruct, VoronoiCells, SimulationBox, Output, DumpInfo, update_perimeters!, update_areas!
    using Random: MersenneTwister

    @testset "Single Square Cell" begin
        N = 1
        particles_data = SelfPropelledVoronoi.VoronoiCells([0.0], [0.0], [0.0], [0.0], [0.0], [0.0])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0)
        # Need to define DumpInfo for ParameterStruct
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))
        
        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        square_vertices = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = square_vertices
        output = SelfPropelledVoronoi.Output()

        SelfPropelledVoronoi.update_perimeters!(parameters, arrays, output)
        @test arrays.perimeters[1] ≈ 4.0

        SelfPropelledVoronoi.update_areas!(parameters, arrays, output)
        @test arrays.areas[1] ≈ 1.0
    end

    @testset "Single Rectangular Cell" begin
        N = 1
        particles_data = SelfPropelledVoronoi.VoronoiCells([0.0], [0.0], [0.0], [0.0], [0.0], [0.0])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0)
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))

        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        rect_vertices = [SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(2.0, 1.0), SVector(0.0, 1.0)]
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = rect_vertices
        output = SelfPropelledVoronoi.Output()

        SelfPropelledVoronoi.update_perimeters!(parameters, arrays, output)
        @test arrays.perimeters[1] ≈ 6.0

        SelfPropelledVoronoi.update_areas!(parameters, arrays, output)
        @test arrays.areas[1] ≈ 2.0
    end

    @testset "Two Cells (Square and Rectangle)" begin
        N = 2
        particles_data = SelfPropelledVoronoi.VoronoiCells(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0)
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))
        
        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        square_vertices = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
        rect_vertices = [SVector(0.0, 0.0), SVector(2.0, 0.0), SVector(2.0, 1.0), SVector(0.0, 1.0)]
        
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = square_vertices
        arrays.neighborlist.voronoi_vertex_positions_per_particle[2] = rect_vertices
        output = SelfPropelledVoronoi.Output()

        SelfPropelledVoronoi.update_perimeters!(parameters, arrays, output)
        @test arrays.perimeters[1] ≈ 4.0
        @test arrays.perimeters[2] ≈ 6.0

        SelfPropelledVoronoi.update_areas!(parameters, arrays, output)
        @test arrays.areas[1] ≈ 1.0
        @test arrays.areas[2] ≈ 2.0
    end
end

@testset "compute_energy" begin
    # Removed: using SelfPropelledVoronoi: ParameterStruct, ArrayStruct, VoronoiCells, SimulationBox, Output, DumpInfo, compute_energy
    using Random: MersenneTwister # Ensure MersenneTwister is explicitly imported if not already available globally in test scope

    @testset "Single Square Cell, Zero Energy" begin
        N = 1
        target_area = 1.0
        target_perimeter = 4.0
        K_A = 0.5
        K_P = 0.2
        particles_data = SelfPropelledVoronoi.VoronoiCells([target_perimeter], [target_area], [K_P], [K_A], [0.0], [0.0])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0)
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))
        
        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        square_vertices = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = square_vertices
        output = SelfPropelledVoronoi.Output()

        energy = SelfPropelledVoronoi.compute_energy(parameters, arrays, output)
        @test energy ≈ 0.0
    end

    @testset "Single Square Cell, Non-Zero Energy" begin
        N = 1
        target_area = 0.8
        target_perimeter = 3.5
        K_A = 0.5
        K_P = 0.2
        particles_data = SelfPropelledVoronoi.VoronoiCells([target_perimeter], [target_area], [K_P], [K_A], [0.0], [0.0])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0) # Re-define or ensure it's available
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))
        
        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        square_vertices = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)]
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = square_vertices
        output = SelfPropelledVoronoi.Output()

        energy = SelfPropelledVoronoi.compute_energy(parameters, arrays, output)
        expected_energy = K_A * (1.0 - target_area)^2 + K_P * (4.0 - target_perimeter)^2
        # 0.5 * (1.0 - 0.8)^2 + 0.2 * (4.0 - 3.5)^2 = 0.5 * (0.2)^2 + 0.2 * (0.5)^2
        # = 0.5 * 0.04 + 0.2 * 0.25 = 0.02 + 0.05 = 0.07
        @test energy ≈ expected_energy
    end

    @testset "Two Cells (Square and Rectangle), Mixed Energy Contributions" begin
        N = 2
        # Cell 1 (Square)
        A1_target = 1.0
        P1_target = 4.0
        KA1 = 0.5
        KP1 = 0.2
        # Cell 2 (Rectangle)
        A2_target = 1.8
        P2_target = 5.5
        KA2 = 0.3
        KP2 = 0.1

        particles_data = SelfPropelledVoronoi.VoronoiCells([P1_target, P2_target], [A1_target, A2_target], [KP1, KP2], [KA1, KA2], [0.0, 0.0], [0.0, 0.0])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0) # Re-define or ensure it's available
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 0.1, false, box, particles_data, dump_info, nothing, MersenneTwister(123))

        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        square_vertices = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(1.0, 1.0), SVector(0.0, 1.0)] # Area=1, Perimeter=4
        rect_vertices = [SVector(0.0,0.0), SVector(2.0,0.0), SVector(2.0,1.0), SVector(0.0,1.0)] # Area=2, Perimeter=6
        arrays.neighborlist.voronoi_vertex_positions_per_particle[1] = square_vertices
        arrays.neighborlist.voronoi_vertex_positions_per_particle[2] = rect_vertices
        output = SelfPropelledVoronoi.Output()

        energy = SelfPropelledVoronoi.compute_energy(parameters, arrays, output)
        
        energy_cell1 = KA1 * (1.0 - A1_target)^2 + KP1 * (4.0 - P1_target)^2
        energy_cell2 = KA2 * (2.0 - A2_target)^2 + KP2 * (6.0 - P2_target)^2
        # energy_cell1 = 0.5 * (1.0 - 1.0)^2 + 0.2 * (4.0 - 4.0)^2 = 0.0
        # energy_cell2 = 0.3 * (2.0 - 1.8)^2 + 0.1 * (6.0 - 5.5)^2 
        #            = 0.3 * (0.2)^2 + 0.1 * (0.5)^2 
        #            = 0.3 * 0.04 + 0.1 * 0.25 = 0.012 + 0.025 = 0.037
        expected_total_energy = energy_cell1 + energy_cell2
        @test energy ≈ expected_total_energy
    end
end

@testset "DataStructs Constructors" begin
    # Removed: using SelfPropelledVoronoi: VoronoiCells, SimulationBox, VoronoiNeighborList, ArrayStruct, Output, DumpInfo, ParameterStruct
    using Random: MersenneTwister # Already imported in other testsets, but good for clarity

    @testset "VoronoiCells Constructor" begin
        p0s = [1.0]; A0s = [1.0]; KPs = [0.1]; KAs = [0.2]; f0s = [0.3]; Drs = [0.4]
        vc = SelfPropelledVoronoi.VoronoiCells(p0s, A0s, KPs, KAs, f0s, Drs)
        @test vc.target_perimeters === p0s
        @test vc.target_areas === A0s
        @test vc.K_P === KPs
        @test vc.K_A === KAs
        @test vc.active_force_strengths === f0s
        @test vc.rotational_diffusion_constants === Drs
    end

    @testset "SimulationBox Constructor" begin
        Lx = 10.0; Ly = 20.0
        sb = SelfPropelledVoronoi.SimulationBox(Lx, Ly)
        @test sb.box_sizes == SVector{2, Float64}(Lx, Ly)
    end

    @testset "VoronoiNeighborList Constructor" begin
        N = 3
        vnl = SelfPropelledVoronoi.VoronoiNeighborList(N)
        @test length(vnl.voronoi_neighbors) == N
        @test all(isempty, vnl.voronoi_neighbors)
        @test all(x -> isa(x, Vector{Int64}), vnl.voronoi_neighbors)

        @test length(vnl.voronoi_vertex_indices) == N
        @test all(isempty, vnl.voronoi_vertex_indices)
        @test all(x -> isa(x, Vector{Int}), vnl.voronoi_vertex_indices)
        
        @test length(vnl.voronoi_vertex_positions_per_particle) == N
        @test all(isempty, vnl.voronoi_vertex_positions_per_particle)
        @test all(x -> isa(x, Vector{SVector{2, Float64}}), vnl.voronoi_vertex_positions_per_particle)

        @test length(vnl.voronoi_vertices) == N # This was an error in the original DataStructs.jl, it should be N
                                                # Corrected in my mental model of DataStructs.jl for the test
                                                # Assuming it's a list of N zero vectors based on constructor
        @test all(x -> x == zero(SVector{2, Float64}), vnl.voronoi_vertices)


        @test vnl.cell_centers_that_share_a_vertex == Tuple{Int, Int, Int}[]
        @test vnl.positions_with_pbc == SVector{2, Float64}[]
        @test vnl.position_indices == Int64[]
    end

    @testset "ArrayStruct Constructor" begin
        N = 2
        as = SelfPropelledVoronoi.ArrayStruct(N)
        @test length(as.positions) == N
        @test all(x -> x == zero(SVector{2, Float64}), as.positions)
        @test length(as.old_positions) == N
        @test all(x -> x == zero(SVector{2, Float64}), as.old_positions)
        @test length(as.forces) == N
        @test all(x -> x == zero(SVector{2, Float64}), as.forces)
        @test length(as.old_forces) == N
        @test all(x -> x == zero(SVector{2, Float64}), as.old_forces)

        @test length(as.orientations) == N
        @test all(x -> x == 0.0, as.orientations)
        @test length(as.old_orientations) == N
        @test all(x -> x == 0.0, as.old_orientations)
        @test length(as.areas) == N
        @test all(x -> x == 0.0, as.areas)
        @test length(as.perimeters) == N
        @test all(x -> x == 0.0, as.perimeters)
        @test length(as.random_forces) == N
        @test all(x -> x == 0.0, as.random_forces)

        @test isa(as.neighborlist, SelfPropelledVoronoi.VoronoiNeighborList)
        @test length(as.neighborlist.voronoi_neighbors) == N
    end

    @testset "Output Constructor" begin
        out = SelfPropelledVoronoi.Output()
        @test out.potential_energy == 0.0
        @test out.steps_done == 0
        @test out.N_voronoi_tesselations == 0
    end

    @testset "DumpInfo Constructor" begin
        # Default
        di_default = SelfPropelledVoronoi.DumpInfo()
        @test di_default.save == true
        @test startswith(di_default.filename, "dump_")
        @test endswith(di_default.filename, ".h5")
        @test di_default.when_to_save_array == (0:1000:1000000)
        @test di_default.save_r == true
        @test di_default.save_F == false
        @test di_default.save_Epot == false

        # Custom
        custom_filename = "custom.h5"
        custom_when_to_save = 1:10:100
        di_custom = SelfPropelledVoronoi.DumpInfo(save=false, filename=custom_filename, when_to_save_array=custom_when_to_save, save_r=false, save_F=true, save_Epot=true)
        @test di_custom.save == false
        @test di_custom.filename == custom_filename
        @test di_custom.when_to_save_array == custom_when_to_save
        @test di_custom.save_r == false
        @test di_custom.save_F == true
        @test di_custom.save_Epot == true
    end

    @testset "ParameterStruct Constructor" begin
        N = 1
        particles = SelfPropelledVoronoi.VoronoiCells([1.0], [1.0], [0.1], [0.2], [0.3], [0.4])
        box = SelfPropelledVoronoi.SimulationBox(10.0, 10.0)
        dump_info = SelfPropelledVoronoi.DumpInfo()
        rng = MersenneTwister(123)
        
        dt_val = 0.01
        N_steps_val = 1000
        kBT_val = 0.1
        frictionconstant_val = 1.0
        pbd_val = 0.05
        verbose_val = true
        callback_val = nothing

        ps = SelfPropelledVoronoi.ParameterStruct(N, dt_val, N_steps_val, kBT_val, frictionconstant_val, pbd_val, verbose_val, box, particles, dump_info, callback_val, rng)

        @test ps.N == N
        @test ps.dt == dt_val
        @test ps.N_steps == N_steps_val
        @test ps.kBT == kBT_val
        @test ps.frictionconstant == frictionconstant_val
        @test ps.periodic_boundary_layer_depth == pbd_val
        @test ps.verbose == verbose_val
        @test ps.box === box
        @test ps.particles === particles
        @test ps.dump_info === dump_info
        @test ps.callback === callback_val
        @test ps.RNG === rng
    end
end

@testset "Tesselation Utilities" begin
    # SelfPropelledVoronoi.replace_or_push!, SelfPropelledVoronoi.norm2, SelfPropelledVoronoi.outer are not exported
    # using SelfPropelledVoronoi: replace_or_push!, norm2, outer # This would fail
    # StaticArrays is already imported at the top level of test/runtests.jl

    @testset "replace_or_push!" begin
        arr1 = [1, 2, 3]
        SelfPropelledVoronoi.replace_or_push!(arr1, 4, 2)
        @test arr1 == [1, 4, 3]

        SelfPropelledVoronoi.replace_or_push!(arr1, 5, 4)
        @test arr1 == [1, 4, 3, 5]

        arr2 = [10]
        SelfPropelledVoronoi.replace_or_push!(arr2, 20, 2)
        @test arr2 == [10, 20]

        # Test ArgumentError for index out of bounds (index > length + 1)
        arr3 = [1, 2, 3]
        @test_throws ArgumentError SelfPropelledVoronoi.replace_or_push!(arr3, 6, 6) # Trying to insert at index 6 in array of length 3
        @test_throws ArgumentError SelfPropelledVoronoi.replace_or_push!(arr3, 6, 5) # Trying to insert at index 5 in array of length 3
        # Check original array is unmodified by failed attempts
        @test arr3 == [1, 2, 3]
    end

    @testset "norm2" begin
        v1 = SVector(3.0, 4.0)
        n1 = SelfPropelledVoronoi.norm2(v1)
        @test n1 ≈ 25.0

        v2 = SVector(0.0, 0.0)
        n2 = SelfPropelledVoronoi.norm2(v2)
        @test n2 ≈ 0.0

        v3 = SVector(-1.0, 2.0)
        n3 = SelfPropelledVoronoi.norm2(v3)
        @test n3 ≈ 5.0
        
        v4 = SVector(1.0, 1.0) # Example with non-integer result for norm
        n4 = SelfPropelledVoronoi.norm2(v4) # 1^2 + 1^2 = 2
        @test n4 ≈ 2.0
    end

    @testset "outer" begin
        a1 = SVector(1.0, 2.0)
        b1 = SVector(3.0, 4.0)
        res1 = SelfPropelledVoronoi.outer(a1, b1)
        # Expected: a[1]*b[2] - a[2]*b[1] = 1.0*4.0 - 2.0*3.0 = 4.0 - 6.0 = -2.0
        # The function returns SVector(val, -val) based on its implementation in Tesselation.jl
        # outer(a,b) = SVector(a[1]*b[2] - a[2]*b[1], a[2]*b[1] - a[1]*b[2])
        # So for (1,2) and (3,4): (1*4 - 2*3, 2*3 - 1*4) = (4-6, 6-4) = (-2, 2)
        @test res1 == SVector(-2.0, 2.0)

        a2 = SVector(0.0, 0.0)
        b2 = SVector(1.0, 1.0)
        res2 = SelfPropelledVoronoi.outer(a2, b2)
        # (0*1 - 0*1, 0*1 - 0*1) = (0,0)
        @test res2 == SVector(0.0, 0.0)

        a3 = SVector(5.0, -2.0)
        b3 = SVector(1.0, 3.0)
        res3 = SelfPropelledVoronoi.outer(a3, b3)
        # (5*3 - (-2)*1, (-2)*1 - 5*3) = (15 - (-2), -2 - 15) = (17, -17)
        @test res3 == SVector(17.0, -17.0)
    end
end

@testset "Circumcenter Tests" begin
    # SelfPropelledVoronoi.circumcenter is not exported.
    # StaticArrays.SVector is already imported.

    @testset "Right-Angled Triangle 1" begin
        a = SVector(0.0, 0.0)
        b = SVector(2.0, 0.0)
        c = SVector(0.0, 2.0)
        expected_center = SVector(1.0, 1.0)
        center = SelfPropelledVoronoi.circumcenter(a, b, c)
        @test center ≈ expected_center
    end

    @testset "Right-Angled Triangle 2" begin
        a = SVector(1.0, 1.0)
        b = SVector(3.0, 1.0)
        c = SVector(1.0, 4.0)
        expected_center = SVector(2.0, 2.5)
        center = SelfPropelledVoronoi.circumcenter(a, b, c)
        @test center ≈ expected_center
    end

    @testset "Equilateral Triangle" begin
        s = 2.0
        p1 = SVector(0.0, 0.0)
        p2 = SVector(s, 0.0)
        p3 = SVector(s/2, s * sqrt(3.0)/2) # SVector(1.0, sqrt(3.0))
        expected_center = SVector(s/2, s / (2*sqrt(3.0))) # SVector(1.0, 1.0 / sqrt(3.0))
        center = SelfPropelledVoronoi.circumcenter(p1, p2, p3)
        @test center ≈ expected_center
    end

    @testset "Collinear Points" begin
        a = SVector(0.0, 0.0)
        b = SVector(1.0, 1.0)
        c = SVector(2.0, 2.0)
        center = SelfPropelledVoronoi.circumcenter(a, b, c)
        # For collinear points, D = 0, leading to division by zero.
        # Expect components to be Inf or NaN.
        @test isinf(center[1]) || isnan(center[1])
        @test isinf(center[2]) || isnan(center[2])
    end
end

@testset "Sort Indices Counter Clockwise" begin
    # SelfPropelledVoronoi.sort_indices_counter_clockwise is not exported.
    # StaticArrays.SVector is already imported.

    @testset "Simple Square" begin
        voronoi_center = SVector(0.5, 0.5)
        v1 = SVector(0.0, 0.0) # Index 1 in voronoi_vertices
        v2 = SVector(1.0, 0.0) # Index 2
        v3 = SVector(1.0, 1.0) # Index 3
        v4 = SVector(0.0, 1.0) # Index 4
        
        voronoi_vertices = [v1, v2, v3, v4]
        # Input order of indices (pointing to voronoi_vertices): 3, 1, 4, 2
        # Corresponding positions: v3, v1, v4, v2
        voronoi_vertex_indices_input = [3, 1, 4, 2] 
        voronoi_vertex_positions_per_particle_input = [v3, v1, v4, v2]

        expected_sorted_indices = [1, 2, 3, 4] # Sorted by angle: v1, v2, v3, v4
        expected_sorted_positions = [v1, v2, v3, v4]

        sorted_indices, sorted_positions = SelfPropelledVoronoi.sort_indices_counter_clockwise(
            voronoi_vertex_indices_input, 
            voronoi_vertex_positions_per_particle_input, 
            voronoi_vertices, 
            voronoi_center, 
            10.0, 10.0 # Lx, Ly are dummy
        )
        
        @test sorted_indices == expected_sorted_indices
        @test sorted_positions == expected_sorted_positions
    end

    @testset "More Complex Order, Different Center" begin
        voronoi_center = SVector(0.0, 0.0)
        pA = SVector(1.0, 0.0)  # Index 1
        pB = SVector(0.0, 1.0)  # Index 2
        pC = SVector(-1.0, 0.0) # Index 3
        pD = SVector(0.0, -1.0) # Index 4

        voronoi_vertices = [pA, pB, pC, pD]
        # Input order of indices: 3, 1, 4, 2
        # Corresponding positions: pC, pA, pD, pB
        voronoi_vertex_indices_input = [3, 1, 4, 2]
        voronoi_vertex_positions_per_particle_input = [pC, pA, pD, pB]

        # atan(0,1)=0 (pA), atan(1,0)=pi/2 (pB), atan(0,-1)=pi (pC), atan(-1,0)=-pi/2 (pD)
        # Sorted by angle: pD (-pi/2), pA (0), pB (pi/2), pC (pi)
        # The atan function returns values in [-pi, pi].
        # Expected sorted order: pD, pA, pB, pC
        expected_sorted_indices = [4, 1, 2, 3] 
        expected_sorted_positions = [pD, pA, pB, pC]

        sorted_indices, sorted_positions = SelfPropelledVoronoi.sort_indices_counter_clockwise(
            voronoi_vertex_indices_input,
            voronoi_vertex_positions_per_particle_input,
            voronoi_vertices,
            voronoi_center,
            10.0, 10.0 # Lx, Ly are dummy
        )
        
        @test sorted_indices == expected_sorted_indices
        @test sorted_positions == expected_sorted_positions
    end

    @testset "Empty Input" begin
        voronoi_center = SVector(0.0, 0.0)
        voronoi_vertices = SVector{2,Float64}[]
        voronoi_vertex_indices_input = Int[]
        voronoi_vertex_positions_per_particle_input = SVector{2,Float64}[]
        
        expected_sorted_indices = Int[]
        expected_sorted_positions = SVector{2,Float64}[]

        sorted_indices, sorted_positions = SelfPropelledVoronoi.sort_indices_counter_clockwise(
            voronoi_vertex_indices_input,
            voronoi_vertex_positions_per_particle_input,
            voronoi_vertices,
            voronoi_center,
            10.0, 10.0 # Lx, Ly are dummy
        )

        @test sorted_indices == expected_sorted_indices
        @test sorted_positions == expected_sorted_positions
    end

    @testset "Single Vertex" begin
        voronoi_center = SVector(0.0, 0.0)
        v = SVector(1.0, 1.0) # Index 1
        voronoi_vertices = [v]
        voronoi_vertex_indices_input = [1]
        voronoi_vertex_positions_per_particle_input = [v]

        expected_sorted_indices = [1]
        expected_sorted_positions = [v]

        sorted_indices, sorted_positions = SelfPropelledVoronoi.sort_indices_counter_clockwise(
            voronoi_vertex_indices_input,
            voronoi_vertex_positions_per_particle_input,
            voronoi_vertices,
            voronoi_center,
            10.0, 10.0 # Lx, Ly are dummy
        )

        @test sorted_indices == expected_sorted_indices
        @test sorted_positions == expected_sorted_positions
    end
end

@testset "Update Positions with PBCs" begin
    # Removed: using SelfPropelledVoronoi: ParameterStruct, ArrayStruct, VoronoiCells, SimulationBox, Output, DumpInfo # update_positions_with_pbcs! is not exported
    using Random: MersenneTwister # Already imported but good for clarity
    # SVector is already imported via StaticArrays at the top level

    # Helper to create ParameterStruct and ArrayStruct for tests
    function setup_pbc_test_structs(N, Lx, Ly, initial_positions)
        box = SelfPropelledVoronoi.SimulationBox(Lx, Ly)
        # Using dummy values for most particle and simulation parameters as they don't affect update_positions_with_pbcs! directly
        particles_data = SelfPropelledVoronoi.VoronoiCells(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N))
        dump_info = SelfPropelledVoronoi.DumpInfo(save=false)
        parameters = SelfPropelledVoronoi.ParameterStruct(N, 0.1, 100, 0.1, 0.1, 2.5, false, box, particles_data, dump_info, nothing, MersenneTwister(123))
        
        arrays = SelfPropelledVoronoi.ArrayStruct(N)
        for i in 1:N
            arrays.positions[i] = initial_positions[i]
        end
        return parameters, arrays, SelfPropelledVoronoi.Output()
    end

    @testset "No PBC images needed" begin
        N = 1; Lx = 10.0; Ly = 10.0
        initial_positions = [SVector(5.0, 5.0)]
        parameters, arrays, output = setup_pbc_test_structs(N, Lx, Ly, initial_positions)

        SelfPropelledVoronoi.update_positions_with_pbcs!(parameters, arrays, output)

        @test length(arrays.neighborlist.positions_with_pbc) == 1
        @test arrays.neighborlist.positions_with_pbc[1] == SVector(5.0, 5.0)
        @test arrays.neighborlist.position_indices == [1]
    end

    @testset "Particle near one edge (left)" begin
        N = 1; Lx = 10.0; Ly = 10.0 # pbc_layer_depth = 2.5
        initial_positions = [SVector(1.0, 5.0)] # x < 2.5
        parameters, arrays, output = setup_pbc_test_structs(N, Lx, Ly, initial_positions)

        SelfPropelledVoronoi.update_positions_with_pbcs!(parameters, arrays, output)
        
        expected_positions = [SVector(1.0, 5.0), SVector(11.0, 5.0)]
        expected_indices = [1, 1]

        @test length(arrays.neighborlist.positions_with_pbc) == 2
        # Order based on implementation: original particle, then its images
        @test arrays.neighborlist.positions_with_pbc == expected_positions
        @test arrays.neighborlist.position_indices == expected_indices
    end

    @testset "Particle near a corner (top-left)" begin
        N = 1; Lx = 10.0; Ly = 10.0 # pbc_layer_depth = 2.5
        initial_positions = [SVector(1.0, 1.0)] # x < 2.5, y < 2.5
        parameters, arrays, output = setup_pbc_test_structs(N, Lx, Ly, initial_positions)

        SelfPropelledVoronoi.update_positions_with_pbcs!(parameters, arrays, output)

        # Original first
        # Then images in order: x+Lx, (x-Lx skipped), y+Ly, (y-Ly skipped), (x+Lx,y+Ly), ...
        # Actual order from code:
        # 1. original: (1,1) idx 1
        # 2. x < pbc_layer_depth (x+Lx): (11,1) idx 1
        # 3. y < pbc_layer_depth (y+Ly): (1,11) idx 1
        # 4. x < pbc_layer_depth && y < pbc_layer_depth (x+Lx, y+Ly): (11,11) idx 1
        
        expected_positions_set = Set([SVector(1.0,1.0), SVector(11.0,1.0), SVector(1.0,11.0), SVector(11.0,11.0)])
        
        @test length(arrays.neighborlist.positions_with_pbc) == 4
        @test Set(arrays.neighborlist.positions_with_pbc) == expected_positions_set
        @test count(==(1), arrays.neighborlist.position_indices) == 4
        # Check exact order for indices if important, for positions Set is fine
        # The function adds original particle i, then all its images, then original particle i+1, etc.
        # Here N=1, so all indices must be 1.
        @test arrays.neighborlist.position_indices == [1, 1, 1, 1]
    end

    @testset "Particle near all four edges (small box)" begin
        N = 1; Lx = 4.0; Ly = 4.0 # pbc_layer_depth = 2.5. Lx-depth = 1.5, Ly-depth = 1.5
        initial_positions = [SVector(2.0, 2.0)] # x < 2.5, x > 1.5, y < 2.5, y > 1.5
        parameters, arrays, output = setup_pbc_test_structs(N, Lx, Ly, initial_positions)

        SelfPropelledVoronoi.update_positions_with_pbcs!(parameters, arrays, output)

        expected_positions_set = Set([
            SVector(2.0, 2.0),   # Original
            SVector(6.0, 2.0),   # x + Lx
            SVector(-2.0, 2.0),  # x - Lx
            SVector(2.0, 6.0),   # y + Ly
            SVector(2.0, -2.0),  # y - Ly
            SVector(6.0, 6.0),   # x + Lx, y + Ly
            SVector(6.0, -2.0),  # x + Lx, y - Ly
            SVector(-2.0, 6.0),  # x - Lx, y + Ly
            SVector(-2.0, -2.0)  # x - Lx, y - Ly
        ])
        
        @test length(arrays.neighborlist.positions_with_pbc) == 9
        @test Set(arrays.neighborlist.positions_with_pbc) == expected_positions_set
        @test count(==(1), arrays.neighborlist.position_indices) == 9
        @test all(idx -> idx == 1, arrays.neighborlist.position_indices)
    end

    @testset "Two particles, one needs images, one not" begin
        N = 2; Lx = 10.0; Ly = 10.0
        p1 = SVector(1.0, 1.0) # Needs images
        p2 = SVector(5.0, 5.0) # No images
        initial_positions = [p1, p2]
        parameters, arrays, output = setup_pbc_test_structs(N, Lx, Ly, initial_positions)

        SelfPropelledVoronoi.update_positions_with_pbcs!(parameters, arrays, output)
        
        # Particle 1 (1.0, 1.0) images: (11,1), (1,11), (11,11)
        # Particle 2 (5.0, 5.0) no images
        
        # Expected order: original p1, original p2, images for p1
        # Code structure:
        # for i in 1:N -> push original positions and indices
        #   arrays.neighborlist.positions_with_pbc = [p1, p2]
        #   arrays.neighborlist.position_indices = [1, 2]
        # for i in 1:N -> add images for particle i
        #   Images for p1 (idx 1): (11,1), (1,11), (11,11) are added
        #   Images for p2 (idx 2): none are added

        expected_positions_set = Set([
            SVector(1.0,1.0), SVector(11.0,1.0), SVector(1.0,11.0), SVector(11.0,11.0), # p1 and its images
            SVector(5.0,5.0) # p2
        ])

        @test length(arrays.neighborlist.positions_with_pbc) == 5
        @test Set(arrays.neighborlist.positions_with_pbc) == expected_positions_set
        
        # Check indices based on implementation logic:
        # Original particles first:
        @test arrays.neighborlist.positions_with_pbc[1] == p1
        @test arrays.neighborlist.positions_with_pbc[2] == p2
        @test arrays.neighborlist.position_indices[1] == 1
        @test arrays.neighborlist.position_indices[2] == 2
        
        # Images for particle 1 are appended next
        remaining_indices = arrays.neighborlist.position_indices[3:end]
        @test count(==(1), remaining_indices) == 3
        @test all(==(1), remaining_indices)
        
        # Overall counts
        @test count(==(1), arrays.neighborlist.position_indices) == 4
        @test count(==(2), arrays.neighborlist.position_indices) == 1
    end
end
