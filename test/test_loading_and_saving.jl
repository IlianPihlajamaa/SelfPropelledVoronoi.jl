
function h5exists(group, name)
    # Check if a group or dataset exists in an HDF5 file
    return name in keys(group)
end

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
        N=N_val, dt=0.01, N_steps=1000, kBT=1.0, frictionconstant=1.0,
        box=box, particles=particles, dump_info=dump_info,
        callback=dummy_callback, RNG=rng
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


@testset "Loading and Saving Tests" begin
    @testset "Dump Tests" begin
        test_filename1 = "test_dump_new.h5"
        try
            rm(test_filename1, force=true) # Cleanup before test
            params1, arrays1, output1 = setup_test_objects(filename_val=test_filename1, steps_val=0)
            
            SelfPropelledVoronoi.save_simulation_state!(params1, arrays1, output1)
            
            @test isfile(test_filename1)
            
            HDF5.h5open(test_filename1, "r") do file
                # Verify parameters
                @test read(file, "parameters/N") == params1.N
                @test read(file, "parameters/dt") == params1.dt
                @test read(file, "parameters/box_sizes") == params1.box.box_sizes |> collect
                @test read(file, "parameters/particles/target_perimeters") == params1.particles.target_perimeters
                @test read(file, "parameters/particles/K_A") == params1.particles.K_A

                # Verify step 0 data
                step0_g = file["0"]
                @test read(step0_g, "positions") ≈ stack(arrays1.positions) # Use ≈ for float comparisons
                @test read(step0_g, "orientations") ≈ arrays1.orientations
                @test read(step0_g, "forces") ≈ arrays1.forces |> stack
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
            SelfPropelledVoronoi.save_simulation_state!(params2_initial, arrays2_initial, output2_initial)

            # Prepare for append (step 100)
            params2_append, arrays2_append, output2_append = setup_test_objects(filename_val=test_filename2, steps_val=100)
            # Modify some data for step 100 to ensure it's different
            arrays2_append.positions[1] = SVector(50.0, 50.0) 
            arrays2_append.areas[1] = 20.0
            output2_append.potential_energy = -50.0

            SelfPropelledVoronoi.save_simulation_state!(params2_append, arrays2_append, output2_append)

            HDF5.h5open(test_filename2, "r") do file
                # Verify parameters are still there (and from initial save)
                @test read(file, "parameters/N") == params2_initial.N 
                
                # Verify step 0 data is still there and unchanged
                step0_g = file["0"]
                @test read(step0_g, "positions") ≈ arrays2_initial.positions |> stack
                @test read(step0_g, "areas") ≈ arrays2_initial.areas
                @test read(step0_g, "potential_energy") ≈ output2_initial.potential_energy


                # Verify step 100 data
                step100_g = file["100"]
                @test read(step100_g, "positions") ≈ arrays2_append.positions |> stack
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
            SelfPropelledVoronoi.save_simulation_state!(params3, arrays3, output3)
            
            HDF5.h5open(test_filename3, "r") do file
                step0_g = file["0"]
                @test !h5exists(step0_g, "positions")
                @test !h5exists(step0_g, "orientations")
                @test !h5exists(step0_g, "forces")
                @test !h5exists(step0_g, "potential_energy")
                @test h5exists(step0_g, "areas")
                @test h5exists(step0_g, "perimeters")
            end
        finally
            rm(test_filename3, force=true)
        end

        test_filename4 = "test_dump.txt" # Unsupported extension
        try
            rm(test_filename4, force=true)
            # No need to capture @warn output for this test, just check file not created
            params4, arrays4, output4 = setup_test_objects(filename_val=test_filename4)
            
            # The function should raise an error 
            @test_throws ArgumentError SelfPropelledVoronoi.save_simulation_state!(params4, arrays4, output4) 
        finally
            @test !isfile(test_filename4) # Ensure file was not created
            rm(test_filename4, force=true) # Should not exist, but good practice
        end
    end
end