using SelfPropelledVoronoi
using Test

function compute_forces_finite_difference(particle_index::Int, parameters::ParameterStruct, arrays::ArrayStruct, output::Output, epsilon::Float64)
    # Get original position
    original_pos = copy(arrays.positions[particle_index])
    
    # Calculate E0 (energy at original position)
    # SelfPropelledVoronoi.voronoi_tesselation!(parameters, arrays, output) # Ensure tesselation is up-to-date
    # E0 = SelfPropelledVoronoi.compute_energy(parameters, arrays, output) # This function is in example/test.jl, need to make it accessible

    # To make compute_energy accessible, we should move it to a source file, e.g. src/AuxiliaryFunctions.jl
    # For now, let's assume it's accessible as SelfPropelledVoronoi.compute_energy

    # --- Potential Energy Function ---
    # This is a simplified view, actual energy computation might involve updating tesselation
    function potential_energy(current_arrays::ArrayStruct)
        # Critical: Ensure that any change in particle position leads to a re-tesselation 
        # before energy calculation if the energy depends on the tesselation.
        temp_arrays = deepcopy(current_arrays) # Use a deepcopy to avoid modifying the original arrays during probing
        SelfPropelledVoronoi.voronoi_tesselation!(parameters, temp_arrays, output) # Re-tesselate based on potentially modified positions
        
        # Call the exported compute_energy function
        return SelfPropelledVoronoi.compute_energy(parameters, temp_arrays, output)
    end

    # Calculate force component Fx
    # Displace in +x
    arrays.positions[particle_index] = original_pos + SVector(epsilon, 0.0)
    E_plus_x = potential_energy(arrays)
    
    # Displace in -x
    arrays.positions[particle_index] = original_pos - SVector(epsilon, 0.0)
    E_minus_x = potential_energy(arrays)
    
    Fx = -(E_plus_x - E_minus_x) / (2 * epsilon)
    
    # Reset position before calculating Fy
    arrays.positions[particle_index] = original_pos
    
    # Calculate force component Fy
    # Displace in +y
    arrays.positions[particle_index] = original_pos + SVector(0.0, epsilon)
    E_plus_y = potential_energy(arrays)
    
    # Displace in -y
    arrays.positions[particle_index] = original_pos - SVector(0.0, epsilon)
    E_minus_y = potential_energy(arrays)
    
    Fy = -(E_plus_y - E_minus_y) / (2 * epsilon)
    
    # Restore original position
    arrays.positions[particle_index] = original_pos
    # It's important to also restore the state of arrays.areas and arrays.perimeters 
    # if they were changed by potential_energy function calls.
    # A final tesselation ensures the main simulation state is correct.
    SelfPropelledVoronoi.voronoi_tesselation!(parameters, arrays, output)


    return SVector(Fx, Fy)
end

function compute_gcm_potential_energy(arrays::ArrayStruct, box_sizes::SVector{2, Float64}, N::Int)
    total_energy = 0.0
    for i in 1:N
        for j in (i+1):N # Iterate over unique pairs
            rij_vec = SelfPropelledVoronoi.compute_pair_distance_vector(arrays.positions[i], arrays.positions[j], box_sizes)
            r_sq = sum(rij_vec .* rij_vec)
            total_energy += exp(-r_sq)
        end
    end
    return total_energy
end

function compute_gcm_forces_finite_difference(particle_index::Int, N_particles::Int, box_sizes::SVector{2,Float64}, arrays::ArrayStruct, epsilon::Float64)
    original_pos = copy(arrays.positions[particle_index])
    
    # Internal function to get energy for GCM
    function gcm_energy_config(current_arrays::ArrayStruct)
        return compute_gcm_potential_energy(current_arrays, box_sizes, N_particles)
    end

    # Calculate force component Fx
    arrays.positions[particle_index] = original_pos + SVector(epsilon, 0.0)
    E_plus_x = gcm_energy_config(arrays)
    
    arrays.positions[particle_index] = original_pos - SVector(epsilon, 0.0)
    E_minus_x = gcm_energy_config(arrays)
    
    Fx = -(E_plus_x - E_minus_x) / (2 * epsilon)
    
    # Reset position
    arrays.positions[particle_index] = original_pos
    
    # Calculate force component Fy
    arrays.positions[particle_index] = original_pos + SVector(0.0, epsilon)
    E_plus_y = gcm_energy_config(arrays)
    
    arrays.positions[particle_index] = original_pos - SVector(0.0, epsilon)
    E_minus_y = gcm_energy_config(arrays)
    
    Fy = -(E_plus_y - E_minus_y) / (2 * epsilon)
    
    # Restore original position
    arrays.positions[particle_index] = original_pos

    return SVector(Fx, Fy)
end

using Random
using StaticArrays

@testset "Force Calculation" begin
    # Simulation parameters
    N = 25
    rho = 1.0
    L = sqrt(N/rho)
    Lx = L
    Ly = L
    dt = 0.001
    Nsteps = 1
    pbc_layer_depth = 2.5

    # Create SimulationBox
    box = SimulationBox(Lx, Ly)

    # Create VoronoiCells parameters
    target_perimeters = 4.0 * ones(N)
    target_areas = 1.0 * ones(N)
    K_P = 1.0 * ones(N)
    K_A = 1.0 * ones(N)
    active_force_strengths = zeros(N)
    D_r = zeros(N)
    voronoi_cells = VoronoiCells(target_perimeters, target_areas, K_P, K_A, active_force_strengths, D_r)

    # Create ParameterStruct
    kBT = 0.0
    frictionconstant = 1.0
    random_seed = 12345
    Random.seed!(random_seed)
    dump_info = DumpInfo(save=false)
    callback = (p,a,o) -> nothing
    parameter_struct = ParameterStruct(N, dt, Nsteps, kBT, frictionconstant, pbc_layer_depth, false, box, voronoi_cells, dump_info, callback, MersenneTwister(random_seed))

    # Create ArrayStruct
    arrays = ArrayStruct(N)

    # Initialize particle positions randomly
    x = rand(Float64, N) .* Lx
    y = rand(Float64, N) .* Ly
    arrays.positions .= SVector.(x, y)

    # Initialize particle orientations
    arrays.orientations .= zeros(N)

    # Create Output object
    output = Output()

    # Perform initial Voronoi tesselation
    SelfPropelledVoronoi.voronoi_tesselation!(parameter_struct, arrays, output)

    # Define a small epsilon for finite differencing
    epsilon = 1e-7

    # Define the number of particles to test
    num_particles_to_test = 5 # Test the first 5 particles

    # Call the SPV force calculation
    SelfPropelledVoronoi.compute_forces_SPV!(parameter_struct, arrays, output)

    # Loop through the selected number of particles to test forces
    for i in 1:num_particles_to_test
        # Get the i-th particle's SPV force
        F_spv = arrays.forces[i]

        # Calculate the finite difference force for the i-th particle
        # Ensure arrays.positions is not modified by compute_forces_finite_difference before this call for particle i
        # The compute_forces_finite_difference function should ideally use a deepcopy of arrays internally if it modifies positions for calculation
        # or ensure it restores the state perfectly. Given its current structure, it restores the specific particle's position.
        # However, the global 'arrays' state (like areas, perimeters updated by tesselation) might be affected by intermediate steps.
        # For a clean test, it's best if compute_forces_finite_difference takes a 'clean' arrays state or deepcopies appropriately.
        # For now, we proceed assuming compute_forces_finite_difference handles its state changes correctly and restores.
        
        # Re-calculate tesselation before finite difference to ensure a clean state for energy calculation baseline
        # This is important if the main SPV force calculation modified shared state like 'output' or 'arrays' in ways
        # that compute_forces_finite_difference doesn't fully reset for its own baseline energy calculation.
        # SelfPropelledVoronoi.voronoi_tesselation!(parameter_struct, arrays, output) # This might be redundant if compute_forces_finite_difference does it.

        F_fd = compute_forces_finite_difference(i, parameter_struct, arrays, output, epsilon)

        # Compare the forces
        @test F_spv[1] ≈ F_fd[1] atol=1e-5
        @test F_spv[2] ≈ F_fd[2] atol=1e-5
    end
end

@testset "GCM Force Calculation" begin
    # Simulation parameters
    N = 25
    rho = 1.0
    L = sqrt(N/rho)
    Lx = L
    Ly = L
    dt = 0.001
    Nsteps = 1
    pbc_layer_depth = 0.0 # GCM doesn't use Voronoi padding

    # Create SimulationBox
    box = SimulationBox(Lx, Ly)

    # Create VoronoiCells parameters (defaults for ParameterStruct)
    target_perimeters = zeros(N)
    target_areas = zeros(N)
    K_P = zeros(N)
    K_A = zeros(N)
    active_force_strengths = zeros(N)
    D_r = zeros(N)
    voronoi_cells = VoronoiCells(target_perimeters, target_areas, K_P, K_A, active_force_strengths, D_r)

    # Create ParameterStruct
    kBT = 0.0
    frictionconstant = 1.0
    random_seed = 54321 # Different seed
    Random.seed!(random_seed)
    dump_info = DumpInfo(save=false)
    callback = (p,a,o) -> nothing
    parameter_struct_gcm = ParameterStruct(N, dt, Nsteps, kBT, frictionconstant, pbc_layer_depth, false, box, voronoi_cells, dump_info, callback, MersenneTwister(random_seed))

    # Create ArrayStruct
    arrays_gcm = ArrayStruct(N)

    # Initialize particle positions randomly
    x_gcm = rand(Float64, N) .* Lx # Use different variable names to avoid potential scope issues if copy-pasting
    y_gcm = rand(Float64, N) .* Ly
    arrays_gcm.positions .= SVector.(x_gcm, y_gcm)

    # Initialize particle orientations
    arrays_gcm.orientations .= zeros(N)

    # Create Output object
    output_gcm = Output()

    # Note: No initial Voronoi tesselation needed for GCM

    # Define a small epsilon for finite differencing
    epsilon_gcm = 1e-7

    # Define the number of particles to test
    num_particles_to_test_gcm = 5 # Test the first 5 particles

    # Call the analytical GCM force calculation
    SelfPropelledVoronoi.compute_forces_GCM!(parameter_struct_gcm, arrays_gcm, output_gcm)

    # Loop through the selected number of particles to test forces
    for i in 1:num_particles_to_test_gcm
        # Get the i-th particle's analytical GCM force
        F_gcm_analytical = arrays_gcm.forces[i]

        # Calculate the finite difference GCM force for the i-th particle
        F_gcm_fd = compute_gcm_forces_finite_difference(i, parameter_struct_gcm.N, parameter_struct_gcm.box.box_sizes, arrays_gcm, epsilon_gcm)

        # Compare the forces
        @test F_gcm_analytical[1] ≈ F_gcm_fd[1] atol=1e-5
        @test F_gcm_analytical[2] ≈ F_gcm_fd[2] atol=1e-5
    end
end
