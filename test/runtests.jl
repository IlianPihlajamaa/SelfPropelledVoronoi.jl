using SelfPropelledVoronoi
using Test
using HDF5
using StaticArrays
import Random

# Helper function to set up mock objects for testing
function setup_test_objects(N_val=2, Lx_val=10.0, Ly_val=10.0;
                            filename_val="default_test.h5", steps_val=0,
                            save_r_val=true, save_F_val=true, save_Epot_val=true)
    
    Random.seed!(1234) # Ensure consistent random values for tests

    particles = VoronoiCells(
        [6.0 for _ in 1:N_val], [1.0 for _ in 1:N_val], # target_perimeters, target_areas
        [1.0 for _ in 1:N_val], [1.0 for _ in 1:N_val], # K_P, K_A
        [1.0 for _ in 1:N_val], [0.1 for _ in 1:N_val]  # active_force_strengths, rotational_diffusion_constants
    )
    box = SimulationBox(Lx_val, Ly_val)
    dump_info = DumpInfo(
        save=true, filename=filename_val, when_to_save_array=0:1:1000, # Adjusted for more steps
        save_r=save_r_val, save_F=save_F_val, save_Epot=save_Epot_val
    )
    
    # ParameterStruct expects a callback and RNG.
    # save_simulation_state! itself doesn't use the callback or RNG directly,
    # but they are part of the struct.
    dummy_callback(params, arrays, output) = nothing 
    rng = Random.MersenneTwister(1234)

    params = ParameterStruct(
        N_val, 0.01, 1000, 1.0, 1.0, 1.0, false, box, particles, dump_info, dummy_callback, rng
    )

    arrays = ArrayStruct(N_val)
    arrays.positions .= [SVector(rand(rng, Float64), rand(rng, Float64)) .* Lx_val for _ in 1:N_val]
    arrays.old_positions .= deepcopy(arrays.positions)
    arrays.orientations .= [rand(rng, Float64) * 2pi for _ in 1:N_val]
    arrays.old_orientations .= deepcopy(arrays.orientations)
    arrays.forces .= [SVector(rand(rng, Float64), rand(rng, Float64)) for _ in 1:N_val]
    arrays.old_forces .= deepcopy(arrays.forces)
    arrays.areas .= [rand(rng, Float64) for _ in 1:N_val]
    arrays.perimeters .= [rand(rng, Float64) * 6.0 for _ in 1:N_val] # Example perimeter
    
    output = Output()
    output.steps_done = steps_val
    output.potential_energy = rand(rng, Float64) * 100.0

    return params, arrays, output
end

@testset "SelfPropelledVoronoi.jl" begin

    @test SelfPropelledVoronoi.test_force(1) == 2

    @testset "Dump Tests" begin
        test_filename1 = "test_dump_new.h5"
        try
            rm(test_filename1, force=true) # Cleanup before test
            params1, arrays1, output1 = setup_test_objects(filename_val=test_filename1, steps_val=0)
            
            save_simulation_state!(params1, arrays1, output1)
            
            @test isfile(test_filename1)
            
            HDF5.h5open(test_filename1, "r") do file
                # Verify parameters
                @test read(file, "parameters/N") == params1.N
                @test read(file, "parameters/dt") == params1.dt
                @test read(file, "parameters/box_sizes") == params1.box.box_sizes
                @test read(file, "parameters/particles/target_perimeters") == params1.particles.target_perimeters
                @test read(file, "parameters/particles/K_A") == params1.particles.K_A

                # Verify step 0 data
                step0_g = file["0"]
                @test read(step0_g, "positions") ≈ arrays1.positions # Use ≈ for float comparisons
                @test read(step0_g, "orientations") ≈ arrays1.orientations
                @test read(step0_g, "forces") ≈ arrays1.forces
                @test read(step0_g, "potential_energy") ≈ output1.potential_energy
                @test read(step0_g, "areas") ≈ arrays1.areas
                @test read(step0_g, "perimeters") ≈ arrays1.perimeters
            end
        finally
            rm(test_filename1, force=true) # Cleanup after test
        end

        test_filename2 = "test_dump_append.h5"
        try
            rm(test_filename2, force=true)
            # Initial save (step 0)
            params2_initial, arrays2_initial, output2_initial = setup_test_objects(filename_val=test_filename2, steps_val=0)
            save_simulation_state!(params2_initial, arrays2_initial, output2_initial)

            # Prepare for append (step 100)
            params2_append, arrays2_append, output2_append = setup_test_objects(filename_val=test_filename2, steps_val=100)
            # Modify some data for step 100 to ensure it's different
            arrays2_append.positions[1] = SVector(50.0, 50.0) 
            arrays2_append.areas[1] = 20.0
            output2_append.potential_energy = -50.0

            save_simulation_state!(params2_append, arrays2_append, output2_append)

            HDF5.h5open(test_filename2, "r") do file
                # Verify parameters are still there (and from initial save)
                @test read(file, "parameters/N") == params2_initial.N 
                
                # Verify step 0 data is still there and unchanged
                step0_g = file["0"]
                @test read(step0_g, "positions") ≈ arrays2_initial.positions
                @test read(step0_g, "areas") ≈ arrays2_initial.areas
                @test read(step0_g, "potential_energy") ≈ output2_initial.potential_energy


                # Verify step 100 data
                step100_g = file["100"]
                @test read(step100_g, "positions") ≈ arrays2_append.positions
                @test read(step100_g, "areas") ≈ arrays2_append.areas
                @test read(step100_g, "potential_energy") ≈ output2_append.potential_energy
            end
        finally
            rm(test_filename2, force=true)
        end

        test_filename3 = "test_dump_flags.h5"
        try
            rm(test_filename3, force=true)
            params3, arrays3, output3 = setup_test_objects(filename_val=test_filename3, steps_val=0, 
                                                            save_r_val=false, save_F_val=false, save_Epot_val=false)
            save_simulation_state!(params3, arrays3, output3)
            
            HDF5.h5open(test_filename3, "r") do file
                step0_g = file["0"]
                @test !HDF5.exists(step0_g, "positions")
                @test !HDF5.exists(step0_g, "orientations")
                @test !HDF5.exists(step0_g, "forces")
                @test !HDF5.exists(step0_g, "potential_energy")
                @test HDF5.exists(step0_g, "areas")
                @test HDF5.exists(step0_g, "perimeters")
            end
        finally
            rm(test_filename3, force=true)
        end

        test_filename4 = "test_dump.txt" # Unsupported extension
        try
            rm(test_filename4, force=true)
            # No need to capture @warn output for this test, just check file not created
            params4, arrays4, output4 = setup_test_objects(filename_val=test_filename4)
            
            # We expect a warning, but the function should run without error and not create the file.
            # If it errors, the test framework will catch it.
            save_simulation_state!(params4, arrays4, output4) 
            
            @test !isfile(test_filename4)
        finally
            rm(test_filename4, force=true) # Should not exist, but good practice
        end
    end

    @testset "Load Tests" begin
        N_test = 3
        Lx_test = 15.0
        Ly_test = 15.0
        dt_test = 0.005
        
        # Mock parameters to be written to HDF5
        mock_params_dict = Dict(
            "N" => N_test,
            "dt" => dt_test,
            "N_steps" => 10000,
            "kBT" => 0.5,
            "frictionconstant" => 1.2,
            "box_sizes" => SVector{2, Float64}(Lx_test, Ly_test),
            "particles" => Dict(
                "target_perimeters" => [6.5 for _ in 1:N_test],
                "target_areas" => [1.1 for _ in 1:N_test],
                "K_P" => [1.2 for _ in 1:N_test],
                "K_A" => [1.3 for _ in 1:N_test],
                "active_force_strengths" => [0.8 for _ in 1:N_test],
                "rotational_diffusion_constants" => [0.05 for _ in 1:N_test]
            )
        )

        # Helper function to create a test HDF5 file
        function create_test_hdf5_file(filename, params_dict, steps_to_save; 
                                       skip_orient_for_step=nothing, 
                                       skip_forces_for_step=nothing, 
                                       skip_epot_for_step=nothing)
            h5open(filename, "w") do file
                # Save parameters
                g_params = create_group(file, "parameters")
                for (key, value) in params_dict
                    if key == "particles"
                        g_particles = create_group(g_params, "particles")
                        for (pkey, pvalue) in value
                            write(g_particles, pkey, pvalue)
                        end
                    else
                        write(g_params, key, value)
                    end
                end

                # Save data for each specified step
                N = params_dict["N"]
                for step_val in steps_to_save
                    g_step = create_group(file, string(step_val))
                    
                    # Create distinguishable data for each step
                    positions_data = [SVector{2, Float64}(rand() + step_val, rand() + step_val) for _ in 1:N]
                    areas_data = [rand() + step_val for _ in 1:N]
                    perimeters_data = [rand() * 6.0 + step_val for _ in 1:N]
                    
                    write(g_step, "positions", hcat(positions_data...)') # Store as N x 2 matrix
                    write(g_step, "areas", areas_data)
                    write(g_step, "perimeters", perimeters_data)

                    if step_val != skip_orient_for_step
                        orientations_data = [rand() * 2pi + step_val for _ in 1:N]
                        write(g_step, "orientations", orientations_data)
                    end
                    
                    if step_val != skip_forces_for_step
                        forces_data = [SVector{2, Float64}(rand() + step_val, rand() + step_val) for _ in 1:N]
                        write(g_step, "forces", hcat(forces_data...)') # Store as N x 2 matrix
                    end

                    if step_val != skip_epot_for_step
                        potential_energy_data = rand() * 10.0 + step_val
                        write(g_step, "potential_energy", potential_energy_data)
                    end
                end
            end
        end

        test_load_filename = "test_load_functions.h5"
        
        @testset "Load Simulation State" begin
            test_steps = [0, 100] # Save step 0 and step 100
            latest_step = maximum(test_steps)

            try
                create_test_hdf5_file(test_load_filename, mock_params_dict, test_steps)
                
                params_loaded, arrays_loaded, output_loaded = load_simulation_state(test_load_filename)

                # Verify ParameterStruct
                @test params_loaded.N == mock_params_dict["N"]
                @test params_loaded.dt == mock_params_dict["dt"]
                @test params_loaded.box.box_sizes == mock_params_dict["box_sizes"]
                @test params_loaded.particles.target_perimeters == mock_params_dict["particles"]["target_perimeters"]
                @test params_loaded.particles.K_A == mock_params_dict["particles"]["K_A"]
                # Add more param checks as needed

                # Verify ArrayStruct and Output (should be from latest_step)
                @test output_loaded.steps_done == latest_step
                
                # Re-generate latest step data to compare (since rand() makes it unique each time)
                # This is a bit tricky; ideally, we'd read what was written.
                # For simplicity in this test, we check if the structure is correct and values are plausible.
                # A more robust test would read from the HDF5 directly here to get expected values.
                
                h5open(test_load_filename, "r") do file_check
                    latest_step_group = file_check[string(latest_step)]
                    expected_positions_raw = read(latest_step_group, "positions")
                    expected_positions = [SVector{2, Float64}(expected_positions_raw[i,:]) for i in 1:size(expected_positions_raw,1)]
                    @test arrays_loaded.positions ≈ expected_positions

                    expected_orientations = read(latest_step_group, "orientations")
                    @test arrays_loaded.orientations ≈ expected_orientations
                    
                    expected_forces_raw = read(latest_step_group, "forces")
                    expected_forces = [SVector{2, Float64}(expected_forces_raw[i,:]) for i in 1:size(expected_forces_raw,1)]
                    @test arrays_loaded.forces ≈ expected_forces

                    expected_potential_energy = read(latest_step_group, "potential_energy")
                    @test output_loaded.potential_energy ≈ expected_potential_energy
                    
                    expected_areas = read(latest_step_group, "areas")
                    @test arrays_loaded.areas ≈ expected_areas

                    expected_perimeters = read(latest_step_group, "perimeters")
                    @test arrays_loaded.perimeters ≈ expected_perimeters
                end
                
            finally
                rm(test_load_filename, force=true)
            end

            # Test with missing optional fields for the latest step
            @testset "Load Simulation State - Missing Optional Fields" begin
                missing_fields_step = 100
                try
                    create_test_hdf5_file(test_load_filename, mock_params_dict, [0, missing_fields_step],
                                          skip_orient_for_step=missing_fields_step,
                                          skip_forces_for_step=missing_fields_step,
                                          skip_epot_for_step=missing_fields_step)
                    
                    params_loaded, arrays_loaded, output_loaded = load_simulation_state(test_load_filename)
                    
                    @test output_loaded.steps_done == missing_fields_step
                    @test all(x -> x == 0.0, arrays_loaded.orientations) # Should default to zeros
                    @test all(x -> x == SVector(0.0, 0.0), arrays_loaded.forces) # Should default to zero SVectors
                    @test output_loaded.potential_energy == 0.0 # Should default to 0.0

                    # Check that non-missing fields are still loaded correctly
                     h5open(test_load_filename, "r") do file_check
                        latest_step_group = file_check[string(missing_fields_step)]
                        expected_positions_raw = read(latest_step_group, "positions")
                        expected_positions = [SVector{2, Float64}(expected_positions_raw[i,:]) for i in 1:size(expected_positions_raw,1)]
                        @test arrays_loaded.positions ≈ expected_positions
                    end

                finally
                    rm(test_load_filename, force=true)
                end
            end
        end

        @testset "Load Trajectory" begin
            test_steps_traj = [0, 50, 100]
            try
                create_test_hdf5_file(test_load_filename, mock_params_dict, test_steps_traj)
                traj_data, params_loaded_traj = load_trajectory(test_load_filename)

                # Verify ParameterStruct
                @test params_loaded_traj.N == mock_params_dict["N"]
                @test params_loaded_traj.dt == mock_params_dict["dt"]
                # Add more param checks

                # Verify TrajectoryData
                @test traj_data.steps_saved == test_steps_traj
                @test length(traj_data.positions_trajectory) == length(test_steps_traj)
                @test length(traj_data.orientations_trajectory) == length(test_steps_traj)
                @test length(traj_data.forces_trajectory) == length(test_steps_traj)
                @test length(traj_data.potential_energy_trajectory) == length(test_steps_traj)
                @test length(traj_data.areas_trajectory) == length(test_steps_traj)
                @test length(traj_data.perimeters_trajectory) == length(test_steps_traj)

                h5open(test_load_filename, "r") do file_check
                    for (idx, step_val) in enumerate(test_steps_traj)
                        step_group = file_check[string(step_val)]
                        
                        expected_pos_raw = read(step_group, "positions")
                        expected_pos = [SVector{2, Float64}(expected_pos_raw[i,:]) for i in 1:size(expected_pos_raw,1)]
                        @test traj_data.positions_trajectory[idx] ≈ expected_pos

                        expected_orient = read(step_group, "orientations")
                        @test traj_data.orientations_trajectory[idx] ≈ expected_orient

                        expected_forces_raw = read(step_group, "forces")
                        expected_forces = [SVector{2, Float64}(expected_forces_raw[i,:]) for i in 1:size(expected_forces_raw,1)]
                        @test traj_data.forces_trajectory[idx] ≈ expected_forces
                        
                        expected_epot = read(step_group, "potential_energy")
                        @test traj_data.potential_energy_trajectory[idx] ≈ expected_epot

                        expected_areas = read(step_group, "areas")
                        @test traj_data.areas_trajectory[idx] ≈ expected_areas

                        expected_perims = read(step_group, "perimeters")
                        @test traj_data.perimeters_trajectory[idx] ≈ expected_perims
                    end
                end
            finally
                rm(test_load_filename, force=true)
            end

            # Test with missing optional fields for some steps in trajectory
            @testset "Load Trajectory - Missing Optional Fields" begin
                 test_steps_missing_traj = [0, 25, 50] # Step 25 will have missing fields
                try
                    create_test_hdf5_file(test_load_filename, mock_params_dict, test_steps_missing_traj,
                                          skip_orient_for_step=25,
                                          skip_forces_for_step=25,
                                          skip_epot_for_step=nothing) # Epot for step 25 will exist

                    traj_data_missing, _ = load_trajectory(test_load_filename)
                    
                    @test traj_data_missing.steps_saved == test_steps_missing_traj
                    
                    # Check step 0 (all fields present)
                    h5open(test_load_filename, "r") do file_check
                        step0_group = file_check["0"]
                        expected_orient_0 = read(step0_group, "orientations")
                        @test traj_data_missing.orientations_trajectory[1] ≈ expected_orient_0
                        expected_forces_raw_0 = read(step0_group, "forces")
                        expected_forces_0 = [SVector{2, Float64}(expected_forces_raw_0[i,:]) for i in 1:size(expected_forces_raw_0,1)]
                        @test traj_data_missing.forces_trajectory[1] ≈ expected_forces_0
                    end

                    # Check step 25 (idx=2) (orientations and forces missing, Epot present)
                    @test all(x -> x == 0.0, traj_data_missing.orientations_trajectory[2])
                    @test all(x -> x == SVector(0.0,0.0), traj_data_missing.forces_trajectory[2])
                    h5open(test_load_filename, "r") do file_check # Check Epot for step 25
                        step25_group = file_check["25"]
                        expected_epot_25 = read(step25_group, "potential_energy")
                        @test traj_data_missing.potential_energy_trajectory[2] ≈ expected_epot_25
                    end


                    # Check step 50 (idx=3) (all fields present again)
                     h5open(test_load_filename, "r") do file_check
                        step50_group = file_check["50"]
                        expected_orient_50 = read(step50_group, "orientations")
                        @test traj_data_missing.orientations_trajectory[3] ≈ expected_orient_50
                        expected_forces_raw_50 = read(step50_group, "forces")
                        expected_forces_50 = [SVector{2, Float64}(expected_forces_raw_50[i,:]) for i in 1:size(expected_forces_raw_50,1)]
                        @test traj_data_missing.forces_trajectory[3] ≈ expected_forces_50
                    end

                finally
                    rm(test_load_filename, force=true)
                end
            end
        end
    end
end
