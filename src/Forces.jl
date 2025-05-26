
"""
    compute_forces_SPV!(parameters, arrays, output)

Computes the forces on each particle based on the Self-Propelled Voronoi (SPV) model.
The forces arise from deviations of cell areas and perimeters from their target values,
effectively modeling area and perimeter elasticity.

The function first ensures the Voronoi tessellation is current by calling
`verify_tesselation`. If the tessellation is not valid, it's updated by
`voronoi_tesselation!` or `update_delauney_vertices!`. Subsequently, it calls
`update_perimeters!` and `update_areas!` to ensure these geometric properties
are current before force calculation.

The force on particle `i` is calculated as the negative gradient of the total potential
energy with respect to its position `r_i`. This involves complex derivatives of cell
areas and perimeters with respect to vertex positions, and then of vertex positions
with respect to particle positions.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. This provides access to:
    - `parameters.N`: Number of particles.
    - `parameters.particles`: Contains `target_areas`, `target_perimeters`, `K_A` (area spring constant), `K_P` (perimeter spring constant).
- `arrays::ArrayStruct`: The struct holding simulation arrays.
    - `arrays.neighborlist`: Used to access Voronoi neighbor information, vertex indices, vertex positions per particle, and positions with periodic boundary conditions.
    - `arrays.areas`: Current areas of Voronoi cells (read after being updated by `update_areas!`).
    - `arrays.perimeters`: Current perimeters of Voronoi cells (read after being updated by `update_perimeters!`).
    - `arrays.forces`: A vector of `SVector{2, Float64}` where the computed force for each particle will be stored. This array is modified in-place.
- `output::Output`: The simulation output struct. Passed to helper functions like `voronoi_tesselation!`, `update_perimeters!`, etc.

# Notes
- The function modifies `arrays.forces` in-place with the newly computed forces.
- The calculation involves contributions from the particle's own cell energy (`dEi/dri`) and from the energy of its neighboring cells (`dEj/dri`).
- The mathematical details of the force calculation (derivatives of geometry) can be complex and are implemented within the loops.
"""
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
            next = j % length(voronoi_vertices_per_particle[particle_i]) + 1
            prev = j == 1 ? length(voronoi_vertices_per_particle[particle_i]) : j - 1
            
            h_i_j = voronoi_vertices_per_particle[particle_i][j]
            h_i_jp1 = voronoi_vertices_per_particle[particle_i][next]
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
            # for vertex_k in voronoi_indices[particle_j]
            #     if vertex_k in voronoi_indices[particle_i]
            #         Nfound += 1
            #         common_vertices[Nfound] = vertex_k                    
            #     end
            # end
            # if Nfound != 2
            #     error("Error: neighborlist does not contain two common voronoi vertices")
            # end

            # h = voronoi_vertices[common_vertices[1]]
            # g = voronoi_vertices[common_vertices[2]]

            # # h is the first common vertex, g is the second common vertex by counterclockwise order around j

            # dAj_dhx = 0.0
            # dAj_dhy = 0.0
            # dAj_dgx = 0.0
            # dAj_dgy = 0.0
            # dPj_dhx = 0.0
            # dPj_dhy = 0.0
            # dPj_dgx = 0.0
            # dPj_dgy = 0.0

            # dhx_dxi = 0.0
            # dhx_dyi = 0.0
            # dhy_dxi = 0.0
            # dhy_dyi = 0.0
            # dgx_dxi = 0.0
            # dgx_dyi = 0.0
            # dgy_dxi = 0.0
            # dgy_dyi = 0.0
            
            # dAj_dxi = dAj_dhx*dhx_dxi + dAj_dhy*dhy_dxi + 
            #           dAj_dgx*dgx_dxi + dAj_dgy*dgy_dxi
            # dAj_dyi = dAj_dhx*dhx_dyi + dAj_dhy*dhy_dyi + 
            #           dAj_dgx*dgx_dyi + dAj_dgy*dgy_dyi
            # dPj_dxi = dPj_dhx*dhx_dxi + dPj_dhy*dhy_dxi + 
            #           dPj_dgx*dgx_dxi + dPj_dgy*dgy_dxi
            # dPj_dyi = dPj_dhx*dhx_dyi + dPj_dhy*dhy_dyi + 
            #           dPj_dgx*dgx_dyi + dPj_dgy*dgy_dyi
            


            # dEj_dxi = dEj_dAj*dAj_dxi + dEj_dPj*dPj_dxi
            # dEj_dyi = dEj_dAj*dAj_dyi + dEj_dPj*dPj_dyi

            # Fx -= dEj_dxi
            # Fy -= dEj_dyi
            
        end 

        Fi = SVector(Fx , Fy)   
        arrays.forces[particle_i] = Fi

    end
    # error()
    return
end

"""
    compute_forces_GCM!(parameters, arrays, output)

Computes the forces on each particle using a Gaussian Core Model (GCM).
The GCM is a soft, repulsive pair potential where the force between two particles `i` and `j`
is given by `F_ij = -2 * r_ij * exp(-r_ij^2)`, where `r_ij` is the distance vector between them.
This function calculates the net force on each particle by summing these pairwise forces.
Periodic boundary conditions are handled using `compute_pair_distance_vector`.

This function might serve as an alternative or simpler test case for force calculations
compared to the more complex SPV model.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. It provides:
    - `parameters.N`: The number of particles.
    - `parameters.box.box_sizes`: Dimensions of the simulation box, used by `compute_pair_distance_vector`.
- `arrays::ArrayStruct`: The struct holding simulation arrays.
    - `arrays.positions`: A vector of `SVector{2, Float64}` representing the current positions of all particles.
    - `arrays.forces`: A vector of `SVector{2, Float64}` where the computed force for each particle will be stored. This array is modified in-place.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures.

# Notes
- The function modifies `arrays.forces` in-place with the newly computed forces.
- The interaction is considered only if the distance `r` is less than `2^(1/6)` (a common cutoff for this potential, related to its inflection point).
"""
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
            # Gaussian core
            f = -2 * rij  * exp(-r^2)
            Fi += f
        end
        arrays.forces[i] =  Fi
    end
    return
end