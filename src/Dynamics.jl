function compute_forces_SPV!(parameters, arrays, output)
    if !verify_tesselation(parameters, arrays, output)
        voronoi_tesselation!(parameters, arrays, output)
    else
        update_delauney_vertices!(parameters, arrays, output)
    end
    N = parameters.N
    voronoi_neighbors = arrays.neighborlist.voronoi_neighbors
    voronoi_vertices = arrays.neighborlist.voronoi_vertices
    voronoi_vertices_per_particle = arrays.neighborlist.voronoi_vertex_positions_per_particle
    voronoi_indices = arrays.neighborlist.voronoi_vertex_indices
    update_perimeters!(parameters, arrays, output)
    update_areas!(parameters, arrays, output)

    areas = arrays.areas
    perimeters = arrays.perimeters
    target_perimeters = parameters.particles.target_perimeters
    target_areas = parameters.particles.target_areas
    K_P = parameters.particles.K_P
    K_A = parameters.particles.K_A

    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    common_vertices = zeros(Int, 2)


    for particle_i in 1:N
        Fx = 0.0
        Fy = 0.0

        dEi_dAi = 2 * K_A[particle_i] * (areas[particle_i] - target_areas[particle_i])
        dEi_dPi = 2 * K_P[particle_i] * (perimeters[particle_i] - target_perimeters[particle_i])
        # @show particle_i, output.steps_done, dEi_dAi, dEi_dPi
        dAi_dxi = 0.0
        dAi_dyi = 0.0
        dPi_dxi = 0.0
        dPi_dyi = 0.0

        # part 1: -dEi/dri
        for j in eachindex(voronoi_indices[particle_i]) # loop over all h
            h_i_j = voronoi_vertices_per_particle[particle_i][j]
            next = j % length(voronoi_neighbors[particle_i]) + 1
            h_i_jp1 = voronoi_vertices_per_particle[particle_i][next]
            prev = j - 1 == 0 ? length(voronoi_neighbors[particle_i]) : j - 1
            h_i_jm1 = voronoi_vertices_per_particle[particle_i][prev]

            # compute the distance between h_i_j and h_i_jp1
            hijp1_m_hij = h_i_jp1 - h_i_j
            hij_m_hijm1 = h_i_j - h_i_jm1
            hijp1_m_hij_length = sqrt(hijp1_m_hij[1]^2 + hijp1_m_hij[2]^2)
            hij_m_hijm1_length = sqrt(hij_m_hijm1[1]^2 + hij_m_hijm1[2]^2)

            dAi_dhijx = (h_i_jp1[2] - h_i_jm1[2])/2
            dAi_dhijy = (h_i_jm1[1] - h_i_jp1[1])/2
            dPi_dhij = hij_m_hijm1/hij_m_hijm1_length + hijp1_m_hij/hijp1_m_hij_length
            dPi_dhijx = dPi_dhij[1]
            dPi_dhijy = dPi_dhij[2]

            i_a, i_b, i_c = arrays.neighborlist.cell_centers_that_share_a_vertex[j]

            if i_a == particle_i
                particle_l = i_b
                particle_k = i_c
            elseif i_b == particle_i
                particle_l = i_a
                particle_k = i_c
            else
                particle_l = i_a
                particle_k = i_b
            end

            rl = positions_with_pbc[particle_l]
            rk = positions_with_pbc[particle_k]
            ri = positions_with_pbc[particle_i]
            xi = ri[1]
            yi = ri[2]
            xl = rl[1]
            yl = rl[2]
            xk = rk[1]
            yk = rk[2]

            # from mathematica:
            dhijx_dxi = -(((yk - yl)*(2*xi*(xl*(-yi + yk) + xk*(yi - yl)) + xk^2*(-yi + yl) + xi^2*(-yk + yl) + (yi - yk)*(xl^2 + yi*yk - (yi + yk)*yl + yl^2)))/ (2*(xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2))
            dhijx_dyi = ((yk - yl)*(xi^2*(xk - xl) + xk^2*xl + xl*(yi - yk)^2 - xk*(xl^2 + (yi - yl)^2) + xi*(-xk^2 + xl^2 + 2*yi*yk - yk^2 - 2*yi*yl + yl^2)))/ (2*(xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2)
            dhijy_dxi = -1/2*((xk - xl)*(xk^2*(yi - yl) + xi^2*(yk - yl) - (yi - yk)*(xl^2 + yi*yk - (yi + yk)*yl + yl^2) + 2*xi*(xl*(yi - yk) + xk*(-yi + yl))))/ (xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2
            dhijy_dyi = -1/2*((xk - xl)*(xi^2*(xk - xl) + xk^2*xl + xl*(yi - yk)^2 - xk*(xl^2 + (yi - yl)^2) + xi*(-xk^2 + xl^2 + 2*yi*yk - yk^2 - 2*yi*yl + yl^2)))/ (xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2


            dAi_dxi += dAi_dhijx * dhijx_dxi + dAi_dhijy * dhijy_dxi
            dAi_dyi += dAi_dhijx * dhijx_dyi + dAi_dhijy * dhijy_dyi
            dPi_dxi += dPi_dhijx * dhijx_dxi + dPi_dhijy * dhijy_dxi
            dPi_dyi += dPi_dhijx * dhijx_dyi + dPi_dhijy * dhijy_dyi   
            # @show dAi_dxi, dAi_dyi, dPi_dxi, dPi_dyi
        end

        dEi_dxi = dEi_dAi*dAi_dxi + dEi_dPi*dPi_dxi
        dEi_dyi = dEi_dAi*dAi_dyi + dEi_dPi*dPi_dyi
        # @show dEi_dxi, dEi_dyi
        Fx -= dEi_dxi
        Fy -= dEi_dyi

        # part 2: -dEj/dri for all j
        for particle_j in voronoi_neighbors[particle_i]
            originalj = arrays.neighborlist.position_indices[particle_j]
            dEj_dAj = 2 * K_A[originalj] * (areas[originalj] - target_areas[originalj])
            dEj_dPj = 2 * K_P[originalj] * (perimeters[originalj] - target_perimeters[originalj])


            #find the two common vertices h and g that occur both in the voronoi_vertices of i and of j
            Nfound = 0
            for vertex_k in voronoi_indices[particle_i]
                if vertex_k in voronoi_indices[particle_j]
                    Nfound += 1
                    common_vertices[Nfound] = vertex_k                    
                end
            end
            if Nfound != 2
                error("Error: neighborlist does not contain two common voronoi vertices")
            end

            h = voronoi_vertices[common_vertices[1]]
            g = voronoi_vertices[common_vertices[2]]


            dAj_dxi = 0.0
            dAj_dyi = 0.0
            dPj_dxi = 0.0
            dPj_dyi = 0.0
            


            dEj_dxi = dEj_dAj*dAj_dxi + dEj_dPj*dPj_dxi
            dEj_dyi = dEj_dAj*dAj_dyi + dEj_dPj*dPj_dyi

            Fx -= dEj_dxi
            Fy -= dEj_dyi
            
        end 

        Fi = SVector(Fx , Fy)   
        arrays.forces[particle_i] = Fi

    end
    # error()
    return
end




function compute_forces_GCM!(parameters, arrays, output)
    # do Gaussian core model for now
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


    old_forces .= forces
        # Check if the simulation should be retesselated

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