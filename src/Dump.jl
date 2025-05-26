using HDF5

# ParameterStruct, ArrayStruct, Output are expected to be defined elsewhere
# and available in the scope where this function is called.
# DumpInfo is part of ParameterStruct.

"""
    save_simulation_state!(parameters::ParameterStruct, arrays::ArrayStruct, output::Output)

Saves the current state of the simulation to an HDF5 file.

The function uses the filename specified in `parameters.dump_info.filename`. It only supports
files with `.h5` or `.hdf5` extensions.

**File Structure:**

If the HDF5 file does not exist, it will be created, and initial simulation parameters
will be saved under the `/parameters` group. These include:
- `/parameters/N`: Number of particles
- `/parameters/dt`: Timestep size
- `/parameters/N_steps`: Total number of simulation steps
- `/parameters/kBT`: Boltzmann temperature
- `/parameters/frictionconstant`: Friction constant
- `/parameters/box_sizes`: Dimensions of the simulation box (SVector)
- `/parameters/particles/target_perimeters`: Target perimeters for Voronoi cells
- `/parameters/particles/target_areas`: Target areas for Voronoi cells
- `/parameters/particles/K_P`: Perimeter stiffness constants
- `/parameters/particles/K_A`: Area stiffness constants
- `/parameters/particles/active_force_strengths`: Active force strengths
- `/parameters/particles/rotational_diffusion_constants`: Rotational diffusion constants

For both new and existing files, the simulation data for the current step (obtained from
`output.steps_done`) is saved in a new group named after the step number (e.g., `/0`,
`/1000`, etc.). If a group for the current step already exists in an existing file, it will
be overwritten with a warning.

Each step group (`/<step_number>/`) contains the following datasets, depending on the
flags in `parameters.dump_info`:
- `positions`: Vector of SVector{2, Float64} (saved if `dump_info.save_r` is true)
- `orientations`: Vector{Float64} (saved if `dump_info.save_r` is true)
- `forces`: Vector of SVector{2, Float64} (saved if `dump_info.save_F` is true)
- `potential_energy`: Float64 (saved if `dump_info.save_Epot` is true)
- `areas`: Vector{Float64} (always saved)
- `perimeters`: Vector{Float64} (always saved)

**Usage Example:**
```julia
# Assuming `params`, `arrs`, `outs` are populated ParameterStruct, ArrayStruct, and Output instances
# and params.dump_info is configured.
save_simulation_state!(params, arrs, outs)
```
The function returns `nothing`.
"""
function save_simulation_state!(parameters::ParameterStruct, arrays::ArrayStruct, output::Output)
    dump_info = parameters.dump_info
    filename = dump_info.filename
    current_step_str = string(output.steps_done)
    file_exists = isfile(filename)

    @info "Attempting to save state to $filename. Step: $current_step_str. File exists: $file_exists."

    if endswith(filename, ".h5") || endswith(filename, ".hdf5")
        if !file_exists
            @info "File $filename does not exist. Creating and writing parameters and initial step."
            HDF5.h5open(filename, "w") do file_handle
                # Save Parameters
                params_g = HDF5.create_group(file_handle, "parameters")
                params_g["N"] = parameters.N
                params_g["dt"] = parameters.dt
                params_g["N_steps"] = parameters.N_steps
                params_g["kBT"] = parameters.kBT
                params_g["frictionconstant"] = parameters.frictionconstant
                params_g["periodic_boundary_layer_depth"] = parameters.periodic_boundary_layer_depth
                # Assuming parameters.box.box_sizes is an SVector or similar directly writable by HDF5.jl
                params_g["box_sizes"] = parameters.box.box_sizes

                # Save particle-specific parameters from VoronoiCells
                # Assuming parameters.particles is of type VoronoiCells
                particles_g = HDF5.create_group(params_g, "particles")
                particles_g["target_perimeters"] = parameters.particles.target_perimeters
                particles_g["target_areas"] = parameters.particles.target_areas
                particles_g["K_P"] = parameters.particles.K_P
                particles_g["K_A"] = parameters.particles.K_A
                particles_g["active_force_strengths"] = parameters.particles.active_force_strengths
                particles_g["rotational_diffusion_constants"] = parameters.particles.rotational_diffusion_constants

                # Save Initial Step Data
                step_g = HDF5.create_group(file_handle, current_step_str)
                if dump_info.save_r
                    # Assuming arrays.positions and arrays.orientations are suitable for HDF5.jl
                    # (e.g., Vector{SVector{2, Float64}} for positions, Vector{Float64} for orientations)
                    step_g["positions"] = arrays.positions
                    step_g["orientations"] = arrays.orientations
                end
                if dump_info.save_F
                    # Assuming arrays.forces is suitable (e.g., Vector{SVector{2, Float64}})
                    step_g["forces"] = arrays.forces
                end
                if dump_info.save_Epot
                    step_g["potential_energy"] = output.potential_energy
                end
                step_g["areas"] = arrays.areas
                step_g["perimeters"] = arrays.perimeters
            end
            @info "Finished writing to new file $filename."
        else
            @info "File $filename exists. Appending data for step $current_step_str."
            HDF5.h5open(filename, "r+") do file_handle # Open in read-write mode
                # Check if group for current step already exists
                if HDF5.exists(file_handle, current_step_str)
                    @warn "Step group $current_step_str already exists in $filename. Overwriting."
                    HDF5.delete_object(file_handle, current_step_str) # Or HDF5.o_delete for older HDF5.jl versions
                end
                
                step_g = HDF5.create_group(file_handle, current_step_str)
                if dump_info.save_r
                    step_g["positions"] = arrays.positions
                    step_g["orientations"] = arrays.orientations
                end
                if dump_info.save_F
                    step_g["forces"] = arrays.forces
                end
                if dump_info.save_Epot
                    step_g["potential_energy"] = output.potential_energy
                end
                step_g["areas"] = arrays.areas
                step_g["perimeters"] = arrays.perimeters
            end
            @info "Finished appending step $current_step_str to $filename."
        end
    else
        @warn "Unsupported file extension: $filename. Only .h5 or .hdf5 are supported. Skipping save."
    end

    return nothing
end