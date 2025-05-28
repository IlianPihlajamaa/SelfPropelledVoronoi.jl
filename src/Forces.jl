
function circumcenter_derivative(ri, rj, rm)
    xi, yi = ri
    xj, yj = rj
    xm, ym = rm
    # from mathematica:
    dhx_dxi = -1/2*((yj - ym)*(2*xi*(xm*(-yi + yj) + xj*(yi - ym)) + xj^2*(-yi + ym) + xi^2*(-yj + ym) + 
    (yi - yj)*(xm^2 + yi*yj - (yi + yj)*ym + ym^2)))/ (xm*(-yi + yj) + xj*(yi - ym) + xi*(-yj + ym))^2

    dhx_dyi = ((yj - ym)*(xi^2*(xj - xm) + xj^2*xm + xm*(yi - yj)^2 - xj*(xm^2 + (yi - ym)^2) + xi*(-xj^2 + xm^2 + 2*yi*yj - yj^2 -  2*yi*ym + ym^2)))/ (2*(xm*(-yi + yj) + xj*(yi - ym) + xi*(-yj + ym))^2)

    dhy_dxi = -1/2*((xj - xm)*(xj^2*(yi - ym) + xi^2*(yj - ym) - (yi - yj)*(xm^2 + yi*yj - (yi + yj)*ym + ym^2) +  2*xi*(xm*(yi - yj) + xj*(-yi + ym))))/ (xm*(-yi + yj) + xj*(yi - ym) + xi*(-yj + ym))^2

    dhy_dyi = -1/2*((xj - xm)*(xi^2*(xj - xm) + xj^2*xm + xm*(yi - yj)^2 -  xj*(xm^2 + (yi - ym)^2) + xi*(-xj^2 + xm^2 + 2*yi*yj -  yj^2 - 2*yi*ym + ym^2)))/  (xm*(-yi + yj) + xj*(yi - ym) + xi*(-yj + ym))^2

    return dhx_dxi, dhy_dxi, dhx_dyi, dhy_dyi
end


function dAi_dhij(hi, hnext, hprev)
    # hi, hnext, hprev are the positions of the vertices in the Voronoi cell
    # returns the derivatives of the area with respect to the vertex positions
    dAi_dhijx = (hnext[2] - hprev[2])/2
    dAi_dhijy = (hprev[1] - hnext[1])/2
    return dAi_dhijx, dAi_dhijy
end

function dPi_dhij(hi, hnext, hprev)
    # hi, hnext, hprev are the positions of the vertices in the Voronoi cell
    # returns the derivatives of the perimeter with respect to the vertex positions
    hijp1_m_hij = hnext - hi
    hij_m_hijm1 = hi - hprev
    hijp1_m_hij_length = sqrt(hijp1_m_hij[1]^2 + hijp1_m_hij[2]^2)
    hij_m_hijm1_length = sqrt(hij_m_hijm1[1]^2 + hij_m_hijm1[2]^2)

    dPi_dhij = hij_m_hijm1/hij_m_hijm1_length - hijp1_m_hij/hijp1_m_hij_length
    dPi_dhijx = dPi_dhij[1]
    dPi_dhijy = dPi_dhij[2]
    return dPi_dhijx, dPi_dhijy
end





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
    N = parameters.N
    voronoi_neighbors = arrays.neighborlist.voronoi_neighbors
    voronoi_vertices = arrays.neighborlist.voronoi_vertices
    voronoi_vertices_per_particle = arrays.neighborlist.voronoi_vertex_positions_per_particle
    voronoi_indices = arrays.neighborlist.voronoi_vertex_indices

    if verify_tessellation(parameters, arrays, output) == false
        voronoi_tesselation!(parameters, arrays, output)
    end

    update_perimeters!(parameters, arrays, output)
    update_areas!(parameters, arrays, output)

    areas = arrays.areas
    perimeters = arrays.perimeters
    target_perimeters = parameters.particles.target_perimeters
    target_areas = parameters.particles.target_areas
    K_P = parameters.particles.K_P
    K_A = parameters.particles.K_A

    positions_with_pbc = arrays.neighborlist.positions_with_pbc


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
            # hijp1_m_hij = h_i_jp1 - h_i_j
            # hij_m_hijm1 = h_i_j - h_i_jm1
            # hijp1_m_hij_length = sqrt(hijp1_m_hij[1]^2 + hijp1_m_hij[2]^2)
            # hij_m_hijm1_length = sqrt(hij_m_hijm1[1]^2 + hij_m_hijm1[2]^2)

            # dAi_dhijx = (h_i_jp1[2] - h_i_jm1[2])/2
            # dAi_dhijy = (h_i_jm1[1] - h_i_jp1[1])/2
            # dPi_dhij = hij_m_hijm1/hij_m_hijm1_length + hijp1_m_hij/hijp1_m_hij_length
            # dPi_dhijx = dPi_dhij[1]
            # dPi_dhijy = dPi_dhij[2]
            dAi_dhijx, dAi_dhijy = dAi_dhij(h_i_j, h_i_jp1, h_i_jm1)
            dPi_dhijx, dPi_dhijy = dPi_dhij(h_i_j, h_i_jp1, h_i_jm1)

            i_a, i_b, i_c = arrays.neighborlist.delaunay_facet_triplets[voronoi_indices[particle_i][j]]
            if i_a == particle_i
                particle_l = i_b
                particle_k = i_c
            elseif i_b == particle_i
                particle_l = i_a
                particle_k = i_c
            elseif i_c == particle_i
                particle_l = i_a
                particle_k = i_b
            else
                error("Error: particle_i should be part of the triplet (particle_i, particle_l, particle_k) that shares the vertex h_i_j")
            end

            rl = positions_with_pbc[particle_l]
            rk = positions_with_pbc[particle_k]
            ri = positions_with_pbc[particle_i]

            dhijx_dxi, dhijy_dxi, dhijx_dyi, dhijy_dyi = circumcenter_derivative(ri, rk, rl)

            dAi_dxi += dAi_dhijx * dhijx_dxi + dAi_dhijy * dhijy_dxi
            dAi_dyi += dAi_dhijx * dhijx_dyi + dAi_dhijy * dhijy_dyi
            dPi_dxi += dPi_dhijx * dhijx_dxi + dPi_dhijy * dhijy_dxi
            dPi_dyi += dPi_dhijx * dhijx_dyi + dPi_dhijy * dhijy_dyi   
            # @show dAi_dxi, dAi_dyi, dPi_dxi, dPi_dyi
        end

        dEi_dxi = dEi_dAi*dAi_dxi + dEi_dPi*dPi_dxi
        dEi_dyi = dEi_dAi*dAi_dyi + dEi_dPi*dPi_dyi

        # we test dAi_dxi by comparing it to the finite difference approximation
        # posi = arrays.positions[particle_i]
        # posi_plus = posi + SVector(1e-6, 0.0)
        # posi_minus = posi - SVector(1e-6, 0.0)

        # arrays.positions[particle_i] = posi_plus
        # voronoi_tesselation!(parameters, arrays, output)
        # update_perimeters!(parameters, arrays, output)
        # update_areas!(parameters, arrays, output)
        # Ei_plus = compute_energy_i(particle_i, parameters, arrays, output)
        # Ai_plus = compute_area_i(particle_i, parameters, arrays, output)   
        # Pi_plus = compute_perimeter_i(particle_i, parameters, arrays, output)
        # arrays.positions[particle_i] = posi_minus
        # voronoi_tesselation!(parameters, arrays, output)
        # update_perimeters!(parameters, arrays, output)
        # update_areas!(parameters, arrays, output)
        # Ei_minus = compute_energy_i(particle_i, parameters, arrays, output)
        # Ai_minus = compute_area_i(particle_i, parameters, arrays, output)
        # Pi_minus = compute_perimeter_i(particle_i, parameters, arrays, output)
        # arrays.positions[particle_i] = posi  # restore original position
        # voronoi_tesselation!(parameters, arrays, output)
        # update_perimeters!(parameters, arrays, output)
        # update_areas!(parameters, arrays, output)
        # dEi_dxi_fd = (Ei_plus - Ei_minus)/(2*1e-6)
        # dAi_dxi_fd = (Ai_plus - Ai_minus)/(2*1e-6)
        # dPi_dxi_fd = (Pi_plus - Pi_minus)/(2*1e-6)
        # @show dEi_dxi, dEi_dxi_fd 
        # @show dAi_dxi, dAi_dxi_fd
        # @show dPi_dxi, dPi_dxi_fd
        

        # @show dEi_dxi, dEi_dyi
        Fx -= dEi_dxi
        Fy -= dEi_dyi

        # part 2: -dEj/dri for all j
        for particle_j in voronoi_neighbors[particle_i]
            originalj = arrays.neighborlist.position_indices[particle_j]
            dEj_dAj = 2 * K_A[originalj] * (areas[originalj] - target_areas[originalj])
            dEj_dPj = 2 * K_P[originalj] * (perimeters[originalj] - target_perimeters[originalj])

            # h is the first common vertex, g is the second common vertex by counterclockwise order around j
            # the positions h and g depend on the positions of the particles i, j, h, and g
            # find the two common vertices h and g that occur both in the voronoi_vertices of i and of j
            h_index = 0
            for (ih, vertex_k) in enumerate(voronoi_indices[particle_j])
                if vertex_k in voronoi_indices[particle_i]
                    h_index = ih
                    break
                end
            end
            # h_index is the index of the first common vertex in the voronoi_vertices of particle_j in counterclockwise order
            g_index = h_index + 1
            if g_index > length(voronoi_indices[particle_j])
                g_index = 1
            end
            # if h_index is 1, it is possible that this is really the second vertex in the counterclockwise order
            # if h_index is 1, we need to check if the next vertex is also in the voronoi_vertices of particle_i
            # if not, we need to set g_index to 1 and h_index to length(voronoi_indices[particle_j])
            if h_index == 1
                if !(voronoi_indices[particle_j][2] in voronoi_indices[particle_i])
                    g_index = 1
                    h_index = length(voronoi_indices[particle_j])
                end
            end
            # test if h_index and g_index are valid indices
            if !(voronoi_indices[particle_j][h_index] in voronoi_indices[particle_i] &&
                 voronoi_indices[particle_j][g_index] in voronoi_indices[particle_i])
                @show particle_i, particle_j, h_index, g_index
                @show voronoi_indices[particle_j]
                @show voronoi_indices[particle_i]
                error("Error: h_index and g_index should be valid indices in voronoi_indices[particle_i]")
            end
            

            hprev_index = h_index == 1 ? length(voronoi_indices[particle_j]) : h_index - 1
            gnext_index = g_index == length(voronoi_indices[particle_j]) ? 1 : g_index + 1

            # these are 4 subsequent voronoi vertices around particle_j
            hprev = voronoi_vertices[voronoi_indices[particle_j][hprev_index]]
            h = voronoi_vertices[voronoi_indices[particle_j][h_index]]
            g = voronoi_vertices[voronoi_indices[particle_j][g_index]]
            gnext = voronoi_vertices[voronoi_indices[particle_j][gnext_index]]

            dPj_dhx, dPj_dhy = dPi_dhij(h, g, hprev)
            dPj_dgx, dPj_dgy = dPi_dhij(g, gnext, h)
            dAj_dhx, dAj_dhy = dAi_dhij(h, g, hprev)
            dAj_dgx, dAj_dgy = dAi_dhij(g, gnext, h)

            h_triplet = arrays.neighborlist.delaunay_facet_triplets[voronoi_indices[particle_j][h_index]]
            g_triplet = arrays.neighborlist.delaunay_facet_triplets[voronoi_indices[particle_j][g_index]]
            
            # find the two particles that share the common vertex with particle_i and particle_j
            if !(particle_i in h_triplet && particle_j in h_triplet && particle_i in g_triplet && particle_j in g_triplet)
                @show particle_i, particle_j, h_triplet, g_triplet
                @show arrays.neighborlist.position_indices[collect(h_triplet)], arrays.neighborlist.position_indices[collect(g_triplet)]
                error("Error: particle_i and particle_j should be part of the triplets that share the common vertex h and g")
            end
            # the triplets are (particle_i, particle_j, particle_m) and (particle_i, particle_j, particle_n) where particle_m and particle_n are the two particles that share the common vertex with particle_i and particle_j
            particle_m_idx = findfirst(x -> x != particle_i && x != particle_j, h_triplet) 
            particle_n_idx = findfirst(x -> x != particle_i && x != particle_j, g_triplet)
            if particle_m_idx === nothing || particle_n_idx === nothing
                error("Error: could not find the two vertices that both belong to particle_i and particle_j")
            end
            particle_m = h_triplet[particle_m_idx]
            particle_n = g_triplet[particle_n_idx]
            # compute the derivatives of h and g with respect to the position of particle_i            

            dhx_dxi, dhy_dxi, dhx_dyi, dhy_dyi = circumcenter_derivative(positions_with_pbc[particle_i], positions_with_pbc[particle_j], positions_with_pbc[particle_m])

            dgx_dxi, dgy_dxi, dgx_dyi, dgy_dyi = circumcenter_derivative(positions_with_pbc[particle_i], positions_with_pbc[particle_j], positions_with_pbc[particle_n])

            
            dAj_dxi = dAj_dhx*dhx_dxi + dAj_dhy*dhy_dxi + 
                      dAj_dgx*dgx_dxi + dAj_dgy*dgy_dxi
            dAj_dyi = dAj_dhx*dhx_dyi + dAj_dhy*dhy_dyi + 
                      dAj_dgx*dgx_dyi + dAj_dgy*dgy_dyi
            dPj_dxi = dPj_dhx*dhx_dxi + dPj_dhy*dhy_dxi + 
                      dPj_dgx*dgx_dxi + dPj_dgy*dgy_dxi
            dPj_dyi = dPj_dhx*dhx_dyi + dPj_dhy*dhy_dyi + 
                      dPj_dgx*dgx_dyi + dPj_dgy*dgy_dyi
        

            dEj_dxi = dEj_dAj*dAj_dxi + dEj_dPj*dPj_dxi
            dEj_dyi = dEj_dAj*dAj_dyi + dEj_dPj*dPj_dyi

            Fx -= dEj_dxi
            Fy -= dEj_dyi

            # test dEj_dxi by comparing it to the finite difference approximation
            # posi = arrays.positions[particle_i]
            # posi_plus = posi + SVector(1e-6, 0.0)
            # posi_minus = posi - SVector(1e-6, 0.0)
            # arrays.positions[particle_i] = posi_plus
            # voronoi_tesselation!(parameters, arrays, output)
            # update_perimeters!(parameters, arrays, output)
            # update_areas!(parameters, arrays, output)
            # Ej_plus = compute_energy_i(originalj, parameters, arrays, output)
            # Aj_plus = compute_area_i(originalj, parameters, arrays, output)
            # Pj_plus = compute_perimeter_i(originalj, parameters, arrays, output)
            # arrays.positions[particle_i] = posi_minus
            # voronoi_tesselation!(parameters, arrays, output)
            # update_perimeters!(parameters, arrays, output)
            # update_areas!(parameters, arrays, output)
            # Ej_minus = compute_energy_i(originalj, parameters, arrays, output)
            # Aj_minus = compute_area_i(originalj, parameters, arrays, output)
            # Pj_minus = compute_perimeter_i(originalj, parameters, arrays, output)
            # arrays.positions[particle_i] = posi  # restore original position
            # voronoi_tesselation!(parameters, arrays, output)
            # update_perimeters!(parameters, arrays, output)
            # update_areas!(parameters, arrays, output)
            # dEj_dxi_fd = (Ej_plus - Ej_minus)/(2*1e-6)
            # dAj_dxi_fd = (Aj_plus - Aj_minus)/(2*1e-6)
            # dPj_dxi_fd = (Pj_plus - Pj_minus)/(2*1e-6)
            # @show dEj_dxi, dEj_dxi_fd
            # @show dAj_dxi, dAj_dxi_fd
            # @show dPj_dxi, dPj_dxi_fd
            
        end 
        Fi = SVector(Fx , Fy)   
        arrays.forces[particle_i] = Fi

    end
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