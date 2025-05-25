# using Combinatorics # No longer needed

function replace_or_push!(array, value, index)
    # replace the value at index if it exists, otherwise push the value
    if index <= length(array)
        array[index] = value
    elseif index == length(array) + 1
        push!(array, value)
    else
        ArgumentError("Index out of bounds")
    end
end

function update_positions_with_pbcs!(parameters, arrays, output)
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc

    positions = arrays.positions
    N = parameters.N
    Lx, Ly = parameters.box.box_sizes
    pbc_layer_depth = 2.5
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc
    N_particles = 0

    positions_with_pbc = SVector{2, Float64}[]
    pbc_position_indices = Int[]


    # first N are the real particles
    for i in 1:N
        N_particles += 1
        push!(positions_with_pbc, positions[i])
        push!(pbc_position_indices, i)
    end

    # now add the positions for the periodic boundary conditions
    for i in 1:N
        x, y = positions[i]
        if x < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y))
            push!(pbc_position_indices, i)
        end
        if x > Lx - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y))
            push!(pbc_position_indices, i)
        end
        if y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x, y + Ly))
            push!(pbc_position_indices, i)
        end
        if y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x, y - Ly))
            push!(pbc_position_indices, i)
        end

        if x < pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y + Ly))
            push!(pbc_position_indices, i)

        end
        if x < pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y - Ly))
            push!(pbc_position_indices, i)
        end
        if x > Lx - pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y + Ly))
            push!(pbc_position_indices, i)
        end
        if x > Lx - pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y - Ly))
            push!(pbc_position_indices, i)
        end
    end
    arrays.neighborlist.positions_with_pbc = positions_with_pbc
    arrays.neighborlist.position_indices = pbc_position_indices
    # resize the array to the new size
    # positions_with_pbc = resize!(positions_with_pbc, N_particles)
    # pbc_position_indices = resize!(pbc_position_indices, N_particles)
    return 
end



# for every particle, find all the corresponding voronoi edges

# To do this, we loop over all facets of the triangulation
# and for each facet, we find the corresponding delauney vertices by getting the circumcenter
# of the triangle formed by the three vertices of the facet. We save the delauney vertex in the list and 
# save also the indices of the vertices of the voronoi cell in a list for all the particle
function voronoi_tesselation!(parameters, arrays, output)

    N = parameters.N
    Lx, Ly = parameters.box.box_sizes
    update_positions_with_pbcs!(parameters, arrays, output)

    # Initialize or clear delaunay_facet_triplets
    empty!(arrays.neighborlist.delaunay_facet_triplets)

    # get the delauney triangulation
    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    N_pbc = length(positions_with_pbc)
    voronoi_vertices = SVector{2, Float64}[]
    voronoi_vertex_indices = [Int[] for _ in 1:N_pbc]
    voronoi_vertex_positions_per_particle = [SVector{2, Float64}[] for _ in 1:N_pbc]
    voronoi_neighbors = [Int[] for _ in 1:N_pbc]
    cell_centers_that_share_a_vertex = Tuple{Int, Int, Int}[]


    tri = Quickhull.delaunay(positions_with_pbc)
    delauney_facets = Quickhull.facets(tri)

    for facet in delauney_facets
        i = facet[1] # This is idx1_pbc
        j = facet[2] # This is idx2_pbc
        k = facet[3] # This is idx3_pbc

        # Populate delaunay_facet_triplets with original particle indices
        orig_idx1 = arrays.neighborlist.position_indices[i]
        orig_idx2 = arrays.neighborlist.position_indices[j]
        orig_idx3 = arrays.neighborlist.position_indices[k]
        
        sorted_triplet = NTuple{3, Int}(sort([orig_idx1, orig_idx2, orig_idx3]))
        
        if !(sorted_triplet in arrays.neighborlist.delaunay_facet_triplets)
            push!(arrays.neighborlist.delaunay_facet_triplets, sorted_triplet)
        end
    
        # add these to the voronoi neighborlist for every particle pair, checking if it is already filled
        if !(j in voronoi_neighbors[i])
            push!(voronoi_neighbors[i], j)
        end
        if !(k in voronoi_neighbors[i])
            push!(voronoi_neighbors[i], k)
        end
        if !(i in voronoi_neighbors[j])
            push!(voronoi_neighbors[j], i)
        end
        if !(k in voronoi_neighbors[j])
            push!(voronoi_neighbors[j], k)
        end
        if !(i in voronoi_neighbors[k])
            push!(voronoi_neighbors[k], i)
        end
        if !(j in voronoi_neighbors[k])
            push!(voronoi_neighbors[k], j)
        end

        # compute the voronoi vertices as the circumcenter of the facet
        # and add it to the voronoi vertices list

        voronoi_vertex_position = circumcenter(
            positions_with_pbc[i],
            positions_with_pbc[j],
            positions_with_pbc[k]
        )

        push!(voronoi_vertices, voronoi_vertex_position)
        push!(cell_centers_that_share_a_vertex, (i, j, k))
        # # add the voronoi vertex to the voronoi vertices list

        push!(voronoi_vertex_indices[i], length(voronoi_vertices))
        push!(voronoi_vertex_positions_per_particle[i], voronoi_vertex_position)
        push!(voronoi_vertex_indices[j], length(voronoi_vertices))
        push!(voronoi_vertex_positions_per_particle[j], voronoi_vertex_position)
        push!(voronoi_vertex_indices[k], length(voronoi_vertices))
        push!(voronoi_vertex_positions_per_particle[k], voronoi_vertex_position)
    end

    for particle in 1:parameters.N
        voronoi_center = arrays.positions[particle]
        # sort the voronoi vertex indices counterclockwise
        voronoi_vertex_indices_new, voronoi_vertex_positions_per_particle_new = sort_indices_counter_clockwise(voronoi_vertex_indices[particle], voronoi_vertex_positions_per_particle[particle], voronoi_vertices, voronoi_center, Lx, Ly)
        # replace the voronoi vertex indices with the new ones
        voronoi_vertex_indices[particle] = voronoi_vertex_indices_new
        voronoi_vertex_positions_per_particle[particle] = voronoi_vertex_positions_per_particle_new
    end

    arrays.neighborlist.voronoi_vertices = voronoi_vertices
    arrays.neighborlist.voronoi_neighbors = voronoi_neighbors
    arrays.neighborlist.voronoi_vertex_indices = voronoi_vertex_indices
    arrays.neighborlist.voronoi_vertex_positions_per_particle = voronoi_vertex_positions_per_particle
    arrays.neighborlist.cell_centers_that_share_a_vertex = cell_centers_that_share_a_vertex
    # The field arrays.neighborlist.delaunay_facet_triplets is already updated in the loop.
    return 
end

function sort_indices_counter_clockwise(voronoi_vertex_indices, voronoi_vertex_positions_per_particle, voronoi_vertices, voronoi_center, Lx, Ly)
    # sort the voronoi vertex indices counterclockwise
    # using the angle between the voronoi center and the voronoi vertices
    angles = zeros(Float64, length(voronoi_vertex_indices))
    for (i, voronoi_vertex_index) in enumerate(voronoi_vertex_indices)
        dx = voronoi_vertices[voronoi_vertex_index][1] - voronoi_center[1]
        dy = voronoi_vertices[voronoi_vertex_index][2] - voronoi_center[2]
        angle = atan(dy, dx)
        angles[i] = angle
    end
    sorted_indices = sortperm(angles)
    return voronoi_vertex_indices[sorted_indices], voronoi_vertex_positions_per_particle[sorted_indices]
end

function outer(a::SVector{2, T}, b::SVector{2, T}) where T
    return SVector{2, T}(a[1]*b[2] - a[2]*b[1], a[2]*b[1] - a[1]*b[2])
end

function norm2(v::SVector{2, T}) where T
    return v[1]^2 + v[2]^2
end

function circumcenter(a::SVector{2, T}, b::SVector{2, T}, c::SVector{2, T}) where T
    ax, ay = a
    bx, by = b
    cx, cy = c
    alength2 = norm2(a)
    blength2 = norm2(b)
    clength2 = norm2(c)

    D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = 1/D*(alength2 * (by - cy) + blength2 * (cy - ay) + clength2 * (ay - by))
    uy = 1/D*(alength2 * (cx - bx) + blength2 * (ax - cx) + clength2 * (bx - ax))

    return SVector{2, Float64}(ux, uy)
end



function verify_tessellation(parameters, arrays, output)
    epsilon = 1e-9
    # N = parameters.N # Not directly used in main loop, but good for context
    
    # Iterate through each pre-computed Delaunay facet triplet
    for triplet_orig_indices in arrays.neighborlist.delaunay_facet_triplets
        p1_idx, p2_idx, p3_idx = triplet_orig_indices

        # Fetch positions directly from arrays.positions using original indices
        pos1 = arrays.positions[p1_idx]
        pos2 = arrays.positions[p2_idx]
        pos3 = arrays.positions[p3_idx]

        # Calculate circumcenter and circumradius squared
        C = circumcenter(pos1, pos2, pos3)
        R_sq = norm2(pos1 - C) # Radius squared from first point of triplet to center

        # Identify particles to check (Neighbors of Neighbors - NoN)
        test_particle_indices = Set{Int}()
        for orig_particle_idx_in_triplet in triplet_orig_indices
            # Neighbors of orig_particle_idx_in_triplet are stored as PBC indices
            # arrays.neighborlist.voronoi_neighbors is indexed by original particle index (1 to N)
            # up to parameters.N. The actual list of neighbors might be longer if N_pbc > N.
            # The problem description for voronoi_neighbors_list[i] in the previous version of verify_tessellation
            # stated: "Primary particles are indexed 1 to N. Their direct entries in voronoi_neighbors_list 
            # (which is sized for N_pbc) correspond to these."
            # This implies voronoi_neighbors should be indexed up to N.
            # Let's assume voronoi_neighbors is indexed 1..N for primary particles.
            if orig_particle_idx_in_triplet <= parameters.N && orig_particle_idx_in_triplet <= length(arrays.neighborlist.voronoi_neighbors)
                for neighbor_pbc_idx in arrays.neighborlist.voronoi_neighbors[orig_particle_idx_in_triplet]
                    # Ensure neighbor_pbc_idx is valid for position_indices
                    if neighbor_pbc_idx > 0 && neighbor_pbc_idx <= length(arrays.neighborlist.position_indices)
                        original_neighbor_idx = arrays.neighborlist.position_indices[neighbor_pbc_idx]
                        push!(test_particle_indices, original_neighbor_idx)
                    else
                        # This case might indicate an issue with neighbor_pbc_idx from voronoi_neighbors list
                        # Or that position_indices is not fully populated for all pbc indices encountered.
                        # Depending on system guarantees, one might error or warn here.
                        # For now, skip if index is invalid.
                        # println("Warning: Invalid neighbor_pbc_idx $neighbor_pbc_idx encountered for triplet particle $orig_particle_idx_in_triplet")
                        continue
                    end
                end
            end
        end

        # Perform the check for each identified test particle
        for p_test_orig_idx in test_particle_indices
            # Ensure p_test_orig_idx is not one of the triplet vertices themselves
            if p_test_orig_idx == p1_idx || p_test_orig_idx == p2_idx || p_test_orig_idx == p3_idx
                continue
            end

            # Fetch the test particle's position
            # Ensure p_test_orig_idx is a valid index for arrays.positions
            if p_test_orig_idx <= 0 || p_test_orig_idx > length(arrays.positions)
                # This indicates an issue with the original_neighbor_idx obtained.
                # println("Warning: Invalid p_test_orig_idx $p_test_orig_idx derived.")
                continue
            end
            p_test = arrays.positions[p_test_orig_idx]

            # Calculate squared distance to circumcenter
            d_sq = norm2(p_test - C)

            # Check for Delaunay violation
            if d_sq < R_sq - epsilon
                return false # Violation found
            end
        end
    end

    # If all checks pass
    arrays.old_positions = deepcopy(arrays.positions)
    return true
end

function update_delauney_vertices!(parameters, arrays, output)
    error()
    return
end
