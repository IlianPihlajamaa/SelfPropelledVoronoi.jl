

function load_simulation_state(filename::String)
    h5open(filename, "r") do file
        # Read simulation parameters
        params_group = file["parameters"]
        N = read(params_group["N"])
        dt = read(params_group["dt"])
        kBT = read(params_group["kBT"])
        frictionconstant = read(params_group["frictionconstant"])
        periodic_boundary_layer_depth = read(params_group["periodic_boundary_layer_depth"])
        
        box_sizes_vec = read(params_group["box_sizes"])
        box = SimulationBox(box_sizes_vec[1], box_sizes_vec[2])

        particles_group = params_group["particles"]
        target_perimeters = read(particles_group["target_perimeters"])
        target_areas = read(particles_group["target_areas"])
        K_P = read(particles_group["K_P"])
        K_A = read(particles_group["K_A"])
        active_force_strengths = read(particles_group["active_force_strengths"])
        rotational_diffusion_constants = read(particles_group["rotational_diffusion_constants"])
        
        particles = VoronoiCells(target_perimeters, target_areas, K_P, K_A, active_force_strengths, rotational_diffusion_constants)

        # DumpInfo: Using defaults for simplicity, as not all fields are typically saved.
        # Adjust if more DumpInfo fields are saved in the HDF5.
        dump_info = DumpInfo(filename=filename) 

        # RNG and callback: Using defaults
        rng = Random.MersenneTwister()
        callback = nothing

        parameter_struct = ParameterStruct(
            N=N, dt=dt, kBT=kBT, frictionconstant=frictionconstant,
            periodic_boundary_layer_depth=periodic_boundary_layer_depth, verbose=false, box=box,
            particles=particles, dump_info=dump_info, callback=callback, RNG=rng
        )

        # Identify the most recent simulation step
        step_groups = String[]
        for name in keys(file)
            if all(isdigit, name) # Check if the group name is purely numerical
                push!(step_groups, name)
            end
        end

        if isempty(step_groups)
            error("No simulation steps found in HDF5 file: $filename")
        end

        latest_step_str = maximum(step_groups, by=x->parse(Int, x))
        latest_step = parse(Int, latest_step_str)
        step_group = file[latest_step_str]

        # Load data for the latest step
        positions_raw = read(step_group["positions"])
        positions = [SVector{2, Float64}(positions_raw[:, i]) for i in 1:N]

        orientations = Float64[]
        if "orientations" in keys(step_group)
            orientations = read(step_group["orientations"])
        else
            orientations = 2π*rand(Float64, N) # Default if not saved
        end

        forces_raw = Matrix{Float64}(undef, 0,0) # placeholder
        forces = Vector{SVector{2, Float64}}(undef, N)
        if "forces" in keys(step_group)
            forces_raw = read(step_group["forces"])
            forces = [SVector{2, Float64}(forces_raw[:, i]) for i in 1:N]
        else
            forces = [zeros(SVector{2, Float64}) for _ in 1:N] # Default if not saved
        end
        
        potential_energy = 0.0
        if "potential_energy" in keys(step_group)
            potential_energy = read(step_group["potential_energy"])
        end

        areas = read(step_group["areas"])
        perimeters = read(step_group["perimeters"])

        # Reconstruct ArrayStruct
        array_struct = ArrayStruct(N)
        array_struct.positions .= positions
        array_struct.orientations .= orientations
        array_struct.forces .= forces
        array_struct.areas .= areas
        array_struct.perimeters .= perimeters
        # old_positions, old_forces, old_orientations, random_forces remain as initialized by ArrayStruct constructor
        # neighborlist is also initialized by ArrayStruct constructor

        # Reconstruct Output
        output_struct = Output()
        output_struct.potential_energy = potential_energy
        output_struct.steps_done = latest_step
        output_struct.N_voronoi_tessellations = 0 # Not saved, default to 0

        return parameter_struct, array_struct, output_struct
    end
end

function load_trajectory(filename::String)
    h5open(filename, "r") do file
        # Read simulation parameters (similar to load_simulation_state)
        params_group = file["parameters"]
        N = read(params_group["N"])
        dt = read(params_group["dt"])
        kBT = read(params_group["kBT"])
        frictionconstant = read(params_group["frictionconstant"])
        
        box_sizes_vec = read(params_group["box_sizes"])
        box = SimulationBox(box_sizes_vec[1], box_sizes_vec[2])

        particles_group = params_group["particles"]
        target_perimeters = read(particles_group["target_perimeters"])
        target_areas = read(particles_group["target_areas"])
        K_P = read(particles_group["K_P"])
        K_A = read(particles_group["K_A"])
        active_force_strengths = read(particles_group["active_force_strengths"])
        rotational_diffusion_constants = read(particles_group["rotational_diffusion_constants"])
        
        particles = VoronoiCells(target_perimeters, target_areas, K_P, K_A, active_force_strengths, rotational_diffusion_constants)
        
        dump_info = DumpInfo(filename=filename) # Using defaults
        rng = Random.MersenneTwister()
        callback = nothing

        parameter_struct = ParameterStruct(
            N=N, dt=dt, kBT=kBT, frictionconstant=frictionconstant, 
            periodic_boundary_layer_depth=3.0, verbose=false, box=box,
            particles=particles, dump_info=dump_info, callback=callback, RNG=rng
        )

        # Initialize TrajectoryData
        trajectory_data = TrajectoryData()

        # Identify and sort step groups
        step_group_names = String[]
        for name in keys(file)
            if all(isdigit, name)
                push!(step_group_names, name)
            end
        end
        
        # Sort step groups numerically
        sort!(step_group_names, by=x->parse(Int, x))

        if isempty(step_group_names)
            @warn "No simulation step groups found in HDF5 file: $filename. Returning empty TrajectoryData."
            return trajectory_data, parameter_struct # Return empty trajectory and params
        end

        for step_str in step_group_names
            step_group = file[step_str]
            current_step = parse(Int, step_str)

            # Load and append positions
            positions_raw = read(step_group["positions"])
            loaded_positions = [SVector{2, Float64}(positions_raw[:, i]) for i in 1:N]
            push!(trajectory_data.positions_trajectory, loaded_positions)

            # Load and append orientations
            if "orientations" in keys(step_group)
                loaded_orientations = read(step_group["orientations"])
                push!(trajectory_data.orientations_trajectory, loaded_orientations)
            else
                # If not saved, append a vector of zeros (or handle as per desired behavior)
                push!(trajectory_data.orientations_trajectory, zeros(Float64, N))
            end

            # Load and append forces
            if "forces" in keys(step_group)
                forces_raw = read(step_group["forces"])
                loaded_forces = [SVector{2, Float64}(forces_raw[:,i]) for i in 1:N]
                push!(trajectory_data.forces_trajectory, loaded_forces)
            else
                # If not saved, append a vector of zero vectors
                push!(trajectory_data.forces_trajectory, [zeros(SVector{2, Float64}) for _ in 1:N])
            end

            # Load and append potential energy
            if "potential_energy" in keys(step_group)
                loaded_potential_energy = read(step_group["potential_energy"])
                push!(trajectory_data.potential_energy_trajectory, loaded_potential_energy)
            else
                push!(trajectory_data.potential_energy_trajectory, 0.0) # Default if not saved
            end
            
            # Load and append areas
            loaded_areas = read(step_group["areas"])
            push!(trajectory_data.areas_trajectory, loaded_areas)

            # Load and append perimeters
            loaded_perimeters = read(step_group["perimeters"])
            push!(trajectory_data.perimeters_trajectory, loaded_perimeters)
            
            # Append the current step number
            push!(trajectory_data.steps_saved, current_step)
        end

        return trajectory_data, parameter_struct
    end
end
