
"""
    replace_or_push!(array, value, index)

Modifies an `array` in-place by either replacing an element at a specified `index` or
pushing a new `value` onto the end of the array.

The behavior is as follows:
- If `index` is within the current bounds of the `array` (i.e., `1 <= index <= length(array)`),
  the element at `array[index]` is replaced with `value`.
- If `index` is one greater than the current length of the `array` (i.e., `index == length(array) + 1`),
  the `value` is pushed onto the end of the `array` using `push!`.
- If `index` is out of these bounds (e.g., less than 1, or greater than `length(array) + 1`),
  an `ArgumentError` is thrown.

# Arguments
- `array`: The array to be modified. This array is changed in-place.
- `value`: The value to be inserted into the `array` or pushed onto its end.
- `index::Integer`: The index at which to replace the existing element, or the position
  (equal to `length(array) + 1`) where the new element should be pushed.

# Returns
- The function does not explicitly return a value (`nothing` is implicitly returned if an `ArgumentError` is not thrown). The primary effect is the modification of the input `array`.

# Throws
- `ArgumentError`: If `index` is not within the valid range for replacement or pushing (i.e., if `index < 1` or `index > length(array) + 1`).
"""
function replace_or_push!(array, value, index)
    index < 1 && throw(ArgumentError("index must be ≥ 1 (got $index)"))

    if index <= length(array)
        array[index] = value
    elseif index == length(array) + 1
        push!(array, value)
    else
        throw(ArgumentError("index $index > length(a)+1 = $(length(array)+1)"))
    end
end

"""
    update_positions_with_pbcs!(parameters, arrays, output)

Creates an extended list of particle positions that includes the original particles
and their periodic images located within a defined boundary layer around the
simulation box. This extended list is essential for Voronoi tessellation algorithms
to correctly construct Voronoi cells for particles near the boundaries of the
periodic domain.

The function first adds all original particle positions. Then, for each original
particle, it checks if it's close to any of the box boundaries (left, right, bottom, top)
or corners. If a particle is within `pbc_layer_depth` of a boundary, its corresponding
periodic image(s) across that boundary (or corner) are added to the list.

The generated list of positions (`arrays.neighborlist.positions_with_pbc`) and a
corresponding list of original particle indices (`arrays.neighborlist.position_indices`)
are stored in the `arrays.neighborlist` structure.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct, providing:
- `arrays::ArrayStruct`: The struct holding simulation arrays.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures.

# Notes
- The fields `arrays.neighborlist.positions_with_pbc` and `arrays.neighborlist.position_indices` are modified in-place.
"""
function update_positions_with_pbcs!(parameters, arrays, output)
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc

    positions = arrays.positions
    N = parameters.N
    Lx, Ly = parameters.box.box_sizes
    pbc_layer_depth = parameters.periodic_boundary_layer_depth

    # if pbc_layer_depth > min(Lx, Ly)/2
    #     error("pbc_layer_depth must be less than half the box size in both dimensions (Lx=$Lx, Ly=$Ly), but is D_pbc=$pbc_layer_depth")
    # end
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc
    N_particles = 0

    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    pbc_position_indices = arrays.neighborlist.position_indices

    # first N are the real particles
    for i in 1:N
        N_particles += 1
        replace_or_push!(positions_with_pbc, positions[i], N_particles)
        replace_or_push!(pbc_position_indices, i, N_particles)
    end
    # now add the positions for the periodic boundary conditions
    for i in 1:N
        x, y = positions[i]
        if x < pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x + Lx, y), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if x > Lx - pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x - Lx, y), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if y < pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x, y + Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if y > Ly - pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x, y - Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end

        if x < pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x + Lx, y + Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if x < pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x + Lx, y - Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if x > Lx - pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x - Lx, y + Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
        if x > Lx - pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            replace_or_push!(positions_with_pbc, SVector(x - Lx, y - Ly), N_particles)
            replace_or_push!(pbc_position_indices, i, N_particles)
        end
    end
    arrays.neighborlist.N_positions_with_pbc = N_particles
    return 
end

"""
    voronoi_tessellation!(parameters, arrays, output)

Performs a full Voronoi tessellation of the system based on current particle positions.
This function orchestrates several steps to compute the Voronoi diagram, which is
dual to the Delaunay triangulation. The results, including lists of Voronoi neighbors,
vertex positions, and topological information, are stored in `arrays.neighborlist`.

The process involves:
1.  **Update Positions with PBC**: Calls `update_positions_with_pbcs!` to generate an extended list of particle positions, including periodic images necessary for correct tessellation across boundaries. This list is stored in `arrays.neighborlist.positions_with_pbc`.
2.  **Delaunay Triangulation**: Computes the Delaunay triangulation of the extended particle positions using `Quickhull.delaunay`.
3.  **Process Delaunay Facets**: Iterates through each facet (triangle) of the Delaunay triangulation:
    *   The three particles forming the triangle are identified as Voronoi neighbors of each other. These are added to `arrays.neighborlist.voronoi_neighbors`.
    *   The circumcenter of the Delaunay triangle is calculated using the `circumcenter` function. This circumcenter is a Voronoi vertex. It's added to a global list of `voronoi_vertices` and also to `arrays.neighborlist.voronoi_vertex_positions_per_particle` for each of the three particles forming the triangle.
    *   The index of this new Voronoi vertex is added to `arrays.neighborlist.voronoi_vertex_indices` for each of the three particles.
    *   A tuple of the three particle indices forming the triangle (cell centers sharing this Voronoi vertex) is stored in `arrays.neighborlist.cell_centers_that_share_a_vertex`.
4.  **Sort Voronoi Vertices**: For each of the original `N` particles, the collected Voronoi vertices associated with it are sorted in counter-clockwise order around the particle's actual position using `sort_indices_counter_clockwise`. This ensures a consistent representation of Voronoi cells. The sorted vertex indices and positions replace the unsorted ones in `arrays.neighborlist.voronoi_vertex_indices` and `arrays.neighborlist.voronoi_vertex_positions_per_particle` for the original particles.
5.  **Store Results**: All computed lists (global Voronoi vertices, neighbors per particle, vertex indices per particle, vertex positions per particle, and cell centers sharing a vertex) are stored in the corresponding fields of `arrays.neighborlist`.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct.
- `arrays::ArrayStruct`: The struct holding simulation arrays.
- `output::Output`: The simulation output struct. Passed to `update_positions_with_pbcs!`.

# Notes
- This function heavily modifies the `arrays.neighborlist` structure in-place.
- It relies on `Quickhull.jl` for the Delaunay triangulation and helper functions like `update_positions_with_pbcs!`, `circumcenter`, and `sort_indices_counter_clockwise`.
"""
function voronoi_tessellation!(parameters, arrays, output)
    # println("Computing Voronoi tessellation... at step ", output.steps_done)
    update_positions_with_pbcs!(parameters, arrays, output)
    N_positions_with_pbc = arrays.neighborlist.N_positions_with_pbc
    tri = Quickhull.delaunay(@view arrays.neighborlist.positions_with_pbc[1:N_positions_with_pbc])
    delaunay_facets = Quickhull.facets(tri)

    update_voronoi_vertices!(parameters, arrays, output, delaunay_facets)
    output.N_voronoi_tessellations += 1
    return 
end



"""
    norm2(v::SVector{2, T}) where T

Computes the squared Euclidean norm (also known as squared length or squared magnitude)
of a 2D vector `v`.

For a vector `v = [v1, v2]`, this is calculated as `v1^2 + v2^2`.
"""
function norm2(v::SVector{2, T}) where T
    return v[1]^2 + v[2]^2
end

"""
    circumcenter(a::SVector{2, T}, b::SVector{2, T}, c::SVector{2, T}) where T

Calculates the circumcenter of a triangle defined by three 2D vertices `a`, `b`, and `c`.
The circumcenter is the point where the perpendicular bisectors of the triangle's sides intersect,
and it is equidistant from all three vertices. This point is also the center of the
triangle's circumcircle.

The calculation is based on a standard formula for the circumcenter's coordinates (ux, uy)
derived from the coordinates of the vertices (ax, ay), (bx, by), (cx, cy):
D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
ux = (1/D) * ((ax^2 + ay^2) * (by - cy) + (bx^2 + by^2) * (cy - ay) + (cx^2 + cy^2) * (ay - by))
uy = (1/D) * ((ax^2 + ay^2) * (cx - bx) + (bx^2 + by^2) * (ax - cx) + (cx^2 + cy^2) * (bx - ax))
This function uses `norm2` which computes `ax^2 + ay^2` (squared length from origin).

# Arguments
- `a::SVector{2, T}`: The 2D coordinates of the first vertex of the triangle. `T` is its element type.
- `b::SVector{2, T}`: The 2D coordinates of the second vertex of the triangle. `T` is its element type.
- `c::SVector{2, T}`: The 2D coordinates of the third vertex of the triangle. `T` is its element type.

# Returns
- `SVector{2, Float64}`: The 2D coordinates `(ux, uy)` of the circumcenter. The coordinates are always returned as `Float64`, regardless of the input type `T`.

# Notes
- If the three points `a`, `b`, `c` are collinear, the denominator `D` will be zero, leading to `Inf` or `NaN` coordinates. In such cases, a warning is issued, and a small random displacement is added to the points to avoid undefined behavior.
"""
function circumcenter(a::SVector{2, T}, b::SVector{2, T}, c::SVector{2, T}) where T
    ax, ay = a
    bx, by = b
    cx, cy = c
    alength2 = norm2(a)
    blength2 = norm2(b)
    clength2 = norm2(c)

    D = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if abs(D) < eps(T)
        @show D
        @show a, b, c
        @warn "Points are collinear or too close together, circumcenter is undefined. You have really weird cells! Adding random displacement to save the simulations."
        disp1 = rand(SVector{2, Float64}) * sqrt(eps(T))
        disp2 = rand(SVector{2, Float64}) * sqrt(eps(T))
        disp3 = rand(SVector{2, Float64}) * sqrt(eps(T))
        return circumcenter(a + disp1, b + disp2, c + disp3)
    end
    ux = 1/D*(alength2 * (by - cy) + blength2 * (cy - ay) + clength2 * (ay - by))
    uy = 1/D*(alength2 * (cx - bx) + blength2 * (ax - cx) + clength2 * (bx - ax))

    return SVector{2, Float64}(ux, uy)
end

"""
    verify_tessellation(parameters, arrays, output)

Verifies the validity of the Voronoi tessellation by checking if all particles
are correctly positioned with respect to the Delaunay facets of the tessellation.
This function checks that for each Delaunay facet (triangle) defined by a triplet of particles,
the circumcircle of the triangle does not contain any other particle that is not one of the triangle's vertices.
# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct, providing the number of particles `N`.
- `arrays::ArrayStruct`: The struct holding simulation arrays, including positions and neighbor lists.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures.
# Returns
- `Bool`: Returns `true` if the tessellation is valid (no violations found), or `false` if any particle is found within the circumcircle of a Delaunay facet that is not one of its vertices.
# Notes
- The function assumes that the Voronoi tessellation has already been computed and that the necessary fields in `arrays.neighborlist` are populated, including `delaunay_facet_triplets`, `positions_with_pbc`, and `voronoi_neighbors`.
- The field arrays.neighborlist.positions_with_pbc is updated to include periodic boundary conditions.
- If the flag `arrays.neighborlist.check_tessellation` is set to `false`, the function immediately returns `false`, indicating that tessellation verification is not required or has been disabled.
"""
function verify_tessellation(parameters, arrays, output)
    # test if we verify the tessellation
    if !arrays.neighborlist.check_tessellation
        return false
    end

    epsilon = 1e-9
    old_used = view(arrays.neighborlist.position_indices, 1:arrays.neighborlist.N_positions_with_pbc)
    update_positions_with_pbcs!(parameters, arrays, output)
    new_used = view(arrays.neighborlist.position_indices, 1:arrays.neighborlist.N_positions_with_pbc)
    hash(old_used) == hash(new_used) || return false  

    # new ghost particles may have been added, so we need to recompute the tessellation

    # Iterate through each pre-computed Delaunay facet triplet
    for triplet_index in 1:arrays.neighborlist.N_triplets#arrays.neighborlist.delaunay_facet_triplets
        triplet = arrays.neighborlist.delaunay_facet_triplets[triplet_index]
        p1_idx, p2_idx, p3_idx = triplet

        # We only check triplets if at least one of the particles is part of the original set
        if !(p1_idx in 1:parameters.N || p2_idx in 1:parameters.N || p3_idx in 1:parameters.N)
            continue # Skip triplets that do not involve original particles
        end

        # Fetch positions directly from arrays.positions using original indices
        pos1 = arrays.neighborlist.positions_with_pbc[p1_idx]
        pos2 = arrays.neighborlist.positions_with_pbc[p2_idx]
        pos3 = arrays.neighborlist.positions_with_pbc[p3_idx]

        # Calculate circumcenter and circumradius squared
        C = circumcenter(pos1, pos2, pos3)
        R_sq = norm2(pos1 - C) # Radius squared from first point of triplet to center


        # no other particle should be inside the circumcircle of this triplet
        # Iterate through all particles in the neighbor lists of the three particles
        # we check only particles that are in the original set (1:N)

        for particle_i in triplet 
            if particle_i > parameters.N
                continue # Skip particles that are not part of the original set
            end

            # check all the neighbors of the triplet particle 
            N_neighbors = arrays.neighborlist.N_voronoi_neighbors_pp[particle_i]
            for particle_idx in 1:N_neighbors
                particle = arrays.neighborlist.voronoi_neighbors[particle_i][particle_idx]
                # Skip the particles that are part of the triplet itself
                if particle == p1_idx || particle == p2_idx || particle == p3_idx
                    continue
                end

                # Get the position of the neighbor particle
                pos_particle = arrays.neighborlist.positions_with_pbc[particle]

                # Calculate the squared distance from the circumcenter to this particle
                dist_sq = norm2(pos_particle - C)

                # Check if this particle is inside the circumcircle
                if dist_sq < R_sq - epsilon
                    # If it is, we have a violation of the tessellation
                    return false
                end
            end

        end
    end
    update_voronoi_vertices!(parameters, arrays, output)
    return true
end


function  update_voronoi_vertices!(parameters, arrays, output, facets)
    # convert facets to a vector of tuples, place in arrays.neighborlist.delaunay_facet_triplets
    N_triplets = 0 
    for facet in facets
        N_triplets += 1
        triplet = (facet[1], facet[2], facet[3])
        replace_or_push!(arrays.neighborlist.delaunay_facet_triplets, triplet, N_triplets)
    end
    arrays.neighborlist.N_triplets = N_triplets

    # set all remaining facets to (-1, -1, -1 )
    # this is not strictly necessary, but it ensures that we get hard errors if we try to access an index that is not filled
    for i in (N_triplets + 1):length(arrays.neighborlist.delaunay_facet_triplets)
        replace_or_push!(arrays.neighborlist.delaunay_facet_triplets, (-1, -1, -1), i)
    end

    # now we can update the voronoi vertices
    return update_voronoi_vertices!(parameters, arrays, output)
end


function is_in_ghost_box(pos, box_sizes, pbc_layer_depth)
    Lx, Ly = box_sizes
    x, y = pos
    if x < -pbc_layer_depth || x > Lx + pbc_layer_depth || y < -pbc_layer_depth || y > Ly + pbc_layer_depth
        return false
    end
    return true
end


function  update_voronoi_vertices!(parameters, arrays, output)
    box_sizes = parameters.box.box_sizes
    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    position_indices = arrays.neighborlist.position_indices
    N_voronoi_neighbors_pp = arrays.neighborlist.N_voronoi_neighbors_pp
    N_voronoi_vertices_pp = arrays.neighborlist.N_voronoi_vertices_pp
    N_triplets = arrays.neighborlist.N_triplets
    N_positions_with_pbc = arrays.neighborlist.N_positions_with_pbc
    
    empty!(arrays.neighborlist.voronoi_vertices)

    for i in (length(arrays.neighborlist.voronoi_neighbors)+1):N_positions_with_pbc
        replace_or_push!(arrays.neighborlist.voronoi_neighbors, Int[], i)
        replace_or_push!(arrays.neighborlist.voronoi_vertex_indices, Int[], i)
    end

    for i in 1:N_positions_with_pbc
        replace_or_push!(N_voronoi_neighbors_pp, 0, i)
        replace_or_push!(N_voronoi_vertices_pp, 0, i)
    end

    N_voronoi_vertices = 0
    for triplet_index in 1:N_triplets
        triplet = arrays.neighborlist.delaunay_facet_triplets[triplet_index]
        i = triplet[1]
        j = triplet[2]
        k = triplet[3]

        # if any pair is primal and neighbouring an image of itself print and error


        orig_i = position_indices[i]
        orig_j = position_indices[j]
        orig_k = position_indices[k]
        
        if (i == orig_i == orig_j) || (i == orig_i == orig_k) || (j == orig_j == orig_k) || (j == orig_j == orig_i) || (k == orig_k == orig_i) || (k == orig_k == orig_j) 
            @warn "Triplet $triplet_index contains particles that are images of each other: $(position_indices[i]), $(position_indices[j]), $(position_indices[k]). Cells should not neighbor themselves!"
            println("Positions: ", positions_with_pbc[i], positions_with_pbc[j], positions_with_pbc[k])
            continue
        end


        # add these to the voronoi neighborlist for every particle pair, checking if it is already filled
        @views if !(j in arrays.neighborlist.voronoi_neighbors[i][1:N_voronoi_neighbors_pp[i]])
            N_voronoi_neighbors_pp[i] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[i], j, N_voronoi_neighbors_pp[i])
        end
        @views if !(k in arrays.neighborlist.voronoi_neighbors[i][1:N_voronoi_neighbors_pp[i]])
            N_voronoi_neighbors_pp[i] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[i], k, N_voronoi_neighbors_pp[i])
        end
        @views if !(i in arrays.neighborlist.voronoi_neighbors[j][1:N_voronoi_neighbors_pp[j]])
            N_voronoi_neighbors_pp[j] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[j], i, N_voronoi_neighbors_pp[j])
        end
        @views if !(k in arrays.neighborlist.voronoi_neighbors[j][1:N_voronoi_neighbors_pp[j]])
            N_voronoi_neighbors_pp[j] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[j], k, N_voronoi_neighbors_pp[j])
        end
        @views if !(i in arrays.neighborlist.voronoi_neighbors[k][1:N_voronoi_neighbors_pp[k]])
            N_voronoi_neighbors_pp[k] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[k], i, N_voronoi_neighbors_pp[k])
        end
        @views if !(j in arrays.neighborlist.voronoi_neighbors[k][1:N_voronoi_neighbors_pp[k]])
            N_voronoi_neighbors_pp[k] += 1
            replace_or_push!(arrays.neighborlist.voronoi_neighbors[k], j, N_voronoi_neighbors_pp[k]) 
        end

        # compute the voronoi vertices as the circumcenter of the facet
        # and add it to the voronoi vertices list

        voronoi_vertex_position = circumcenter(
            positions_with_pbc[i],
            positions_with_pbc[j],
            positions_with_pbc[k]
        )

        push!(arrays.neighborlist.voronoi_vertices, voronoi_vertex_position)
        N_voronoi_vertices += 1
        # replace_or_push!(arrays.neighborlist.voronoi_vertices, voronoi_vertex_position, N_voronoi_vertices)
        # # add the voronoi vertex to the voronoi vertices list only if it falls within the extended box
        if is_in_ghost_box(voronoi_vertex_position, box_sizes, parameters.periodic_boundary_layer_depth)
  
            N_voronoi_vertices_pp[i] += 1
            replace_or_push!(arrays.neighborlist.voronoi_vertex_indices[i], N_voronoi_vertices, N_voronoi_vertices_pp[i])
            N_voronoi_vertices_pp[j] += 1
            replace_or_push!(arrays.neighborlist.voronoi_vertex_indices[j], N_voronoi_vertices, N_voronoi_vertices_pp[j])
            N_voronoi_vertices_pp[k] += 1
            replace_or_push!(arrays.neighborlist.voronoi_vertex_indices[k], N_voronoi_vertices, N_voronoi_vertices_pp[k])
        end
    end



    for particle in 1:N_positions_with_pbc
        # sort the voronoi vertex indices counterclockwise
        sort_indices_counter_clockwise!(arrays.neighborlist, particle, box_sizes)
    end
end


"""
    function sort_indices_counter_clockwise!(neighborlist, particle)

Sorts the Voronoi vertex indices and their corresponding positions for a specific
particle's Voronoi cell in a counter-clockwise (CCW) order around the cell's center.

The sorting is achieved by calculating the angle of each vertex relative to the
`voronoi_center` (the particle's position) and then sorting based on these angles.
The `atan(dy, dx)` function is used to get angles in the range `(-π, π]`, which
naturally provides an ordering for CCW sorting.

# Arguments
- `neighborlist`: The `NeighborList` struct containing Voronoi data, including vertex indices and positions.
- `particle::Int`: The index of the particle whose Voronoi cell vertices are to be sorted.
"""
function sort_indices_counter_clockwise!(neighborlist, particle, box_sizes)
    N_vertices_particle = neighborlist.N_voronoi_vertices_pp[particle]
    if N_vertices_particle == 0
        return # nothing to sort
    end
    voronoi_vertex_indices = @view neighborlist.voronoi_vertex_indices[particle][1:N_vertices_particle]
    voronoi_vertices = neighborlist.voronoi_vertices
    voronoi_center = neighborlist.positions_with_pbc[particle]
    # sort the voronoi vertex indices counterclockwise
    # using the angle between the voronoi center and the voronoi vertices
    # angles = zeros(Float64, length(voronoi_vertex_indices))
    angles = neighborlist.temp_buffer3
    if length(angles) != N_vertices_particle
        resize!(angles, N_vertices_particle)
    end

    voronoi_vertices_mat = reinterpret(reshape, Float64, voronoi_vertices)
    vor_center_x = voronoi_center[1]
    vor_center_y = voronoi_center[2]
    Lx, Ly = box_sizes
    @turbo for i in 1:N_vertices_particle
        voronoi_vertex_index = voronoi_vertex_indices[i]
        dx = voronoi_vertices_mat[1, voronoi_vertex_index] - vor_center_x
        dy = voronoi_vertices_mat[2, voronoi_vertex_index] - vor_center_y
        dx -= round(dx / Lx) * Lx # probably not necessary, but safe to do
        dy -= round(dy / Ly) * Ly
        angle = atan(dy, dx)
        angles[i] = angle
    end

    # all of this is to sort the voronoi vertex indices by angle without allocating a new array
    buf1 = neighborlist.temp_buffer1
    buf2 = neighborlist.temp_buffer2
    if length(buf1) != N_vertices_particle
        resize!(buf1, N_vertices_particle)
    end
    if length(buf2) != N_vertices_particle
        resize!(buf2, N_vertices_particle)
    end
    sortperm!(buf2, angles)
    for i in 1:N_vertices_particle
        buf1[i] = voronoi_vertex_indices[buf2[i]]
    end
    resize!(neighborlist.voronoi_vertex_indices[particle], N_vertices_particle) # Ensure the output vector is the correct size
    for i in 1:N_vertices_particle
        neighborlist.voronoi_vertex_indices[particle][i] = buf1[i]
    end
    return 
end
