using Test
using SelfPropelledVoronoi
using StaticArrays
using Random
using JLD2

@testset "Restart Functionality" begin

    @testset "Test Case 1: Save and Load Basic Data" begin
        # Setup minimal structs
        original_params = ParameterStruct(
            N = 10,
            dt = 0.01,
            N_steps = 100,
            kBT = 0.1,
            frictionconstant = 1.0,
            box = SimulationBox(10.0, 10.0),
            particles = VoronoiCells(zeros(10), ones(10), ones(10), ones(10), ones(10), ones(10))
        )
        original_params.restart_info.save_restart = true
        original_params.restart_info.restart_filename = "test_restart_save_load.jld2"
        original_params.restart_info.restart_save_interval = 1 # Not directly used here but good for consistency

        original_arrays = ArrayStruct(10)
        original_arrays.positions .= [SVector{2, Float64}(rand(2) .* 10.0) for _ in 1:10]
        original_arrays.orientations .= rand(Float64, 10) .* 2π

        original_output = Output()
        original_output.steps_done = 50
        original_output.potential_energy = -150.0

        # Save the structs
        SelfPropelledVoronoi.save_restart_file(original_params, original_arrays, original_output)

        # Load the structs
        loaded_params, loaded_arrays, loaded_output = SelfPropelledVoronoi.load_restart_file("test_restart_save_load.jld2")

        # Compare fields
        @test loaded_params.N == original_params.N
        @test loaded_params.dt == original_params.dt
        @test loaded_params.N_steps == original_params.N_steps
        @test loaded_params.kBT == original_params.kBT
        @test loaded_params.frictionconstant == original_params.frictionconstant
        @test loaded_params.box.box_sizes == original_params.box.box_sizes
        @test loaded_params.particles.target_perimeters == original_params.particles.target_perimeters
        @test loaded_params.restart_info.restart_filename == original_params.restart_info.restart_filename

        @test loaded_arrays.positions == original_arrays.positions
        @test loaded_arrays.orientations == original_arrays.orientations
        @test isempty(loaded_arrays.forces) # Assuming forces are not saved if not computed, or are zeroed. Adjust if different.
        
        @test loaded_output.steps_done == original_output.steps_done
        @test loaded_output.potential_energy == original_output.potential_energy

        # Clean up
        rm("test_restart_save_load.jld2")
    end

    @testset "Test Case 2: Simulation Resume" begin
        N_test = 5
        initial_seed = 12345
        restart_filename = "sim_restart_test.jld2"
        
        # --- Setup for part 1: Run and save ---
        params1 = ParameterStruct(
            N = N_test,
            dt = 0.001,
            N_steps = 5, # Run for 5 steps first
            kBT = 0.05,
            frictionconstant = 1.0,
            box = SimulationBox(5.0, 5.0),
            particles = VoronoiCells(zeros(N_test), ones(N_test), ones(N_test), ones(N_test), 0.1*ones(N_test), 0.1*ones(N_test)),
            RNG = MersenneTwister(initial_seed)
        )
        params1.restart_info.save_restart = true
        params1.restart_info.restart_filename = restart_filename
        params1.restart_info.restart_save_interval = 5 # Save at step 5
        params1.verbose = false # Keep output clean for tests

        arrays1 = ArrayStruct(N_test)
        Random.seed!(params1.RNG, initial_seed) # Seed for initial positions/orientations consistent with RNG in params
        arrays1.positions .= [SVector{2, Float64}(rand(params1.RNG, 2) .* 5.0) for _ in 1:N_test]
        arrays1.orientations .= rand(params1.RNG, Float64, N_test) .* 2π
        # Copy initial state for the reference simulation later
        initial_positions_clone = deepcopy(arrays1.positions)
        initial_orientations_clone = deepcopy(arrays1.orientations)

        output1 = Output()

        # Run simulation for 5 steps
        run_simulation!(params1, arrays1, output1)
        @test output1.steps_done == 5
        
        # Store state after 5 steps
        positions_after_5_steps = deepcopy(arrays1.positions)
        orientations_after_5_steps = deepcopy(arrays1.orientations)

        # --- Setup for part 2: Load and resume ---
        @test isfile(restart_filename) # Check if restart file was created

        loaded_params, loaded_arrays, loaded_output = SelfPropelledVoronoi.load_restart_file(restart_filename)

        @test loaded_output.steps_done == 5
        @test loaded_arrays.positions ≈ positions_after_5_steps # isapprox for FP
        @test loaded_arrays.orientations ≈ orientations_after_5_steps # isapprox for FP
        @test loaded_params.RNG.seed == params1.RNG.seed # Check if RNG state is somewhat similar (JLD2 might not save perfectly identical state for MersenneTwister internals directly, but seed should be indicative if it was part of struct)

        # Continue simulation from loaded state
        loaded_params.N_steps = 10 # Continue up to 10 steps
        loaded_params.verbose = false
        # The RNG from loaded_params should be used automatically
        run_simulation!(loaded_params, loaded_arrays, loaded_output)
        @test loaded_output.steps_done == 10

        final_positions_resumed = deepcopy(loaded_arrays.positions)
        final_orientations_resumed = deepcopy(loaded_arrays.orientations)

        # --- Setup for part 3: Uninterrupted reference simulation ---
        params_ref = ParameterStruct(
            N = N_test,
            dt = 0.001,
            N_steps = 10, # Run for 10 steps straight
            kBT = 0.05,
            frictionconstant = 1.0,
            box = SimulationBox(5.0, 5.0),
            particles = VoronoiCells(zeros(N_test), ones(N_test), ones(N_test), ones(N_test), 0.1*ones(N_test), 0.1*ones(N_test)),
            RNG = MersenneTwister(initial_seed) # Crucially, use the same initial seed
        )
        params_ref.verbose = false
        # No restart saving for the reference run
        params_ref.restart_info.save_restart = false 


        arrays_ref = ArrayStruct(N_test)
        # Use the *exact same* initial positions and orientations
        arrays_ref.positions .= initial_positions_clone
        arrays_ref.orientations .= initial_orientations_clone
        
        output_ref = Output()

        run_simulation!(params_ref, arrays_ref, output_ref)
        @test output_ref.steps_done == 10

        final_positions_ref = deepcopy(arrays_ref.positions)
        final_orientations_ref = deepcopy(arrays_ref.orientations)

        # Compare final states
        # Due to potential floating point nuances, especially if RNG internal state isn't perfectly restored
        # by JLD2 for MersenneTwister (it saves fields, but not internal buffer states usually),
        # we use isapprox. For this system, it should be very close.
        @test final_positions_resumed ≈ final_positions_ref
        @test final_orientations_resumed ≈ final_orientations_ref

        # Clean up
        if isfile(restart_filename)
            rm(restart_filename)
        end
    end
end
