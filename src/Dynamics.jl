function compute_forces_SPV!(parameters, arrays, output)
    # do lennard jones for now
    # println("Computing forces ", output.steps_done)
    # @time tri = Quickhull.delaunay(arrays.positions)
    # @time vor_vertices = Quickhull.voronoi_centers(tri)
    # @time vor_edges = Quickhull.voronoi_edges(tri)
    # @time vor_edge_points = [(vor_vertices[i], vor_vertices[j]) for (i, j) in vor_edges]

    # for every particle, find all the corresponding voronoi edges

    # To do this, we loop over all facets of the triangulation
    # and for each facet, we find the corresponding voronoi vertices by getting the circumcenter
    # of the triangle formed by the three vertices of the facet. We save the voronoi vertex in the list and 
    # save also the indices of the vertices of the voronoi cell in a list for all the particles

    # facets = Quickhull.facets(tri)

    # voronoi_cells = [Int[] for _ in 1:parameters.N]
    # voronoi_edges = [Int[] for _ in 1:parameters.N]




    N = parameters.N
    for i in 1:N
        Fi = SVector(0.0, 0.0)

        for j in 1:N
            if i == j
                continue
            end
            rij = compute_pair_distance_vector(arrays.positions[i], arrays.positions[j], parameters.box.box_sizes)
            r = sqrt(sum(rij.*rij))
            if r < 2^(1/6)
                # Gaussian core
                f = -2 * rij  * exp(-r^2)
                Fi += f
            end
        end
        arrays.forces[i] =  Fi
    end
    return
end


function compute_forces_GCM!(parameters, arrays, output)
    # do lennard jones for now
    N = parameters.N
    for i in 1:N
        Fi = SVector(0.0, 0.0)

        for j in 1:N
            if i == j
                continue
            end
            rij = compute_pair_distance_vector(arrays.positions[i], arrays.positions[j], parameters.box.box_sizes)
            r = sqrt(sum(rij.*rij))
            if r < 2^(1/6)
                # Gaussian core
                f = -2 * rij  * exp(-r^2)
                Fi += f
            end
        end
        arrays.forces[i] =  Fi
    end
    return
end

function do_time_step(parameters, arrays, output)
    if output.steps_done == 0
        do_time_step_Euler_Murayama(parameters, arrays, output)
    else
        do_time_step_Euler_Heun!(parameters, arrays, output)
    end
end




function do_time_step_Euler_Murayama(parameters, arrays, output)
    dt = parameters.dt
    mobility = 1 / parameters.frictionconstant
    active_force_strengths = parameters.particles.active_force_strengths
    rotational_diffusion_constants = parameters.particles.rotational_diffusion_constants

    positions = arrays.positions
    forces = arrays.forces
    orientations = arrays.orientations
    old_positions = arrays.old_positions
    old_forces = arrays.old_forces
    old_orientations = arrays.old_orientations
    box_sizes = parameters.box.box_sizes
    # compute forces
    compute_forces_SPV!(parameters, arrays, output)
    # update orientations
    for particle in 1:parameters.N
        old_orientation = orientations[particle]
        Dr = rotational_diffusion_constants[particle]
        new_orientation = old_orientation + sqrt(2*dt*Dr) * randn()
        orientations[particle] = new_orientation
    end

    # update positions
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))
        new_position = positions[particle] + dt * (forces[particle] * mobility + active_force_strengths[particle] * old_orientation_vector)
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end


    # Check if the simulation should be retesselated
    if !verify_tesselation(parameters, arrays, output)
        voronoi_tesselation!(parameters, arrays, output)
    else
        update_delauney_vertices!(parameters, arrays, output)
    end

    # update orientations
    old_forces .= forces
    old_positions .= positions
    old_orientations .= orientations

    return
end



function do_time_step_Euler_Heun!(parameters, arrays, output)
    dt = parameters.dt
    mobility = 1 / parameters.frictionconstant
    active_force_strengths = parameters.particles.active_force_strengths
    rotational_diffusion_constants = parameters.particles.rotational_diffusion_constants

    positions = arrays.positions
    forces = arrays.forces
    orientations = arrays.orientations
    old_positions = arrays.old_positions
    old_forces = arrays.old_forces
    old_orientations = arrays.old_orientations

    box_sizes = parameters.box.box_sizes
    # compute forces
    compute_forces_SPV!(parameters, arrays, output)

    # update orientations
    for particle in 1:parameters.N
        old_orientation = orientations[particle]
        Dr = rotational_diffusion_constants[particle]
        new_orientation = old_orientation + sqrt(2*dt*Dr) * randn()
        orientations[particle] = new_orientation
    end

    # predictor step
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))


        # update positions
        new_position = positions[particle] + dt * (forces[particle] * mobility + active_force_strengths[particle] * old_orientation_vector)
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end

    if !verify_tesselation(parameters, arrays, output)
        voronoi_tesselation!(parameters, arrays, output)
    else
        update_delauney_vertices!(parameters, arrays, output)
    end

    old_forces .= forces
    compute_forces_SPV!(parameters, arrays, output)

    # corrector step
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))
        new_orientation = orientations[particle]
        new_orientation_vector = SVector(cos(new_orientation), sin(new_orientation))

        old_position = old_positions[particle]
        old_force = old_forces[particle]
        new_force = forces[particle]
        
        # update positions
        new_position = old_position + dt/2 * (
            old_force * mobility + active_force_strengths[particle] * old_orientation_vector +
            new_force * mobility + active_force_strengths[particle] * new_orientation_vector
        )
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end

    if !verify_tesselation(parameters, arrays, output)
        voronoi_tesselation!(parameters, arrays, output)
    else
        update_delauney_vertices!(parameters, arrays, output)
    end
    # update orientations

    old_positions .= positions
    old_orientations .= orientations

    return
end


function run_simulation!(parameters, arrays, output)
    step = output.steps_done

    voronoi_tesselation!(parameters, arrays, output)

    # Main simulation loop
    while true
        if parameters.verbose && step % 100 == 0
            println("Step: ", step)
        end
        # invoke callback
        parameters.callback(parameters, arrays, output)

        do_time_step(parameters, arrays, output)

        # Check if the simulation should be saved
        if parameters.dump_info.save
            if step in parameters.dump_info.when_to_save_array
                save_simulation_state!(parameters, arrays, output)
            end
        end


        # Check if the simulation should be stopped
        step = step + 1
        output.steps_done = step
        if step > parameters.N_steps
            break
        end
    end

end