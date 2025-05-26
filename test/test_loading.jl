
              
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
              
