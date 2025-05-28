# using Combinatorics # No longer needed

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
    if index <= length(array)
        array[index] = value
    elseif index == length(array) + 1
        push!(array, value)
    else
        ArgumentError("Index out of bounds")
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
    # add positions for to arrays.positions for pbcs
    # only considering a layer of particles with depth D_pbc
    N_particles = 0

    positions_with_pbc = SVector{2, Float64}[]
    pbc_position_indices = Int[]
    pbc_images = Tuple{Int, Int}[]


    # first N are the real particles
    for i in 1:N
        N_particles += 1
        push!(positions_with_pbc, positions[i])
        push!(pbc_position_indices, i)
        push!(pbc_images, (0, 0)) # store the original index for later reference
    end

    # now add the positions for the periodic boundary conditions
    for i in 1:N
        x, y = positions[i]
        if x < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y))
            push!(pbc_position_indices, i)
            push!(pbc_images, (1, 0)) 
        end
        if x > Lx - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y))
            push!(pbc_position_indices, i)
            push!(pbc_images, (-1, 0)) 
        end
        if y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x, y + Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (0, 1)) 
        end
        if y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x, y - Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (0, -1)) 
        end

        if x < pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y + Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (1, 1)) 

        end
        if x < pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x + Lx, y - Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (1, -1)) 
        end
        if x > Lx - pbc_layer_depth && y < pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y + Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (-1, 1)) 
        end
        if x > Lx - pbc_layer_depth && y > Ly - pbc_layer_depth
            N_particles += 1
            push!(positions_with_pbc, SVector(x - Lx, y - Ly))
            push!(pbc_position_indices, i)
            push!(pbc_images, (-1, -1)) 
        end
    end
    arrays.neighborlist.positions_with_pbc = positions_with_pbc
    arrays.neighborlist.position_indices = pbc_position_indices
    arrays.neighborlist.position_pbc_images = pbc_images
    return 
end

"""
    voronoi_tesselation!(parameters, arrays, output)

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
function voronoi_tesselation!(parameters, arrays, output)

    update_positions_with_pbcs!(parameters, arrays, output)
    tri = Quickhull.delaunay(arrays.neighborlist.positions_with_pbc)
    delauney_facets = Quickhull.facets(tri)

    update_voronoi_vertices!(parameters, arrays, output, delauney_facets)
    output.N_voronoi_tesselations += 1
    return 
end

"""
    sort_indices_counter_clockwise(voronoi_vertex_indices, voronoi_vertex_positions_per_particle, voronoi_vertices, voronoi_center, Lx, Ly)

Sorts the Voronoi vertex indices and their corresponding positions for a specific
particle's Voronoi cell in a counter-clockwise (CCW) order around the cell's center.

The sorting is achieved by calculating the angle of each vertex relative to the
`voronoi_center` (the particle's position) and then sorting based on these angles.
The `atan(dy, dx)` function is used to get angles in the range `(-π, π]`, which
naturally provides an ordering for CCW sorting.

# Arguments
- `voronoi_vertex_indices::Vector{Int}`: A vector of integer indices pointing to the global `voronoi_vertices` list. These are the vertices associated with the specific particle's cell before sorting.
- `voronoi_vertex_positions_per_particle::Vector{SVector{2, Float64}}`: A vector of 2D positions (`SVector{2, Float64}`) of the Voronoi vertices, corresponding element-wise to `voronoi_vertex_indices`. This list is also sorted along with `voronoi_vertex_indices`.
- `voronoi_vertices::Vector{SVector{2, Float64}}`: The global list containing the 2D positions of all Voronoi vertices in the system. The indices in `voronoi_vertex_indices` refer to this list.
- `voronoi_center::SVector{2, Float64}`: The 2D position of the center of the Voronoi cell, which is the position of the particle itself. This is the reference point for angle calculations.
- `Lx::Float64`: The width of the simulation box. Currently unused in this function.
- `Ly::Float64`: The height of the simulation box. Currently unused in this function.

# Returns
- `Tuple{Vector{Int}, Vector{SVector{2, Float64}}}`: A tuple containing two new vectors:
    1.  The sorted `voronoi_vertex_indices` for the particle's cell, ordered counter-clockwise.
    2.  The correspondingly sorted `voronoi_vertex_positions_per_particle`.
"""
function sort_indices_counter_clockwise(voronoi_vertex_indices, voronoi_vertex_positions_per_particle, voronoi_vertices, voronoi_center)
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


"""
    norm2(v::SVector{2, T}) where T

Computes the squared Euclidean norm (also known as squared length or squared magnitude)
of a 2D vector `v`.

For a vector `v = [v1, v2]`, this is calculated as `v1^2 + v2^2`.
This is often used in computations where the actual distance (square root) is not
needed, as it avoids a potentially costly square root operation.

# Arguments
- `v::SVector{2, T}`: The 2D input vector. `T` is its element type (e.g., `Float64`).

# Returns
- `::T`: The squared Euclidean norm of the vector `v`, equal to `v[1]^2 + v[2]^2`. The return type is the same as the element type of the input vector.
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
- The function relies on `norm2(v)` to calculate the squared magnitude of vectors `a`, `b`, and `c` (from the origin), which are used in the circumcenter formula.
- If the three points `a`, `b`, `c` are collinear, the denominator `D` will be zero, leading to `Inf` or `NaN` coordinates. This case is not explicitly handled here beyond the standard floating-point behavior.
"""
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

"""
    verify_tesselation(parameters, arrays, output)

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
- If the flag `arrays.neighborlist.check_tesselation` is set to `false`, the function immediately returns `false`, indicating that tessellation verification is not required or has been disabled.
"""
function verify_tessellation(parameters, arrays, output)
    # test if we verify the tessellation
    if !arrays.neighborlist.check_tesselation
        return false
    end

    epsilon = 1e-9
    images_hash = hash(arrays.neighborlist.position_pbc_images)
    indices_hash = hash(arrays.neighborlist.position_indices)
    update_positions_with_pbcs!(parameters, arrays, output)
    
    if images_hash != hash(arrays.neighborlist.position_pbc_images) || indices_hash != hash(arrays.neighborlist.position_indices)
        # new ghost particles may have been added, so we need to recompute the tessellation
        return false
    end
    # Iterate through each pre-computed Delaunay facet triplet
    for triplet_indices in arrays.neighborlist.delaunay_facet_triplets
        p1_idx, p2_idx, p3_idx = triplet_indices

        # Fetch positions directly from arrays.positions using original indices
        pos1 = arrays.neighborlist.positions_with_pbc[p1_idx]
        pos2 = arrays.neighborlist.positions_with_pbc[p2_idx]
        pos3 = arrays.neighborlist.positions_with_pbc[p3_idx]

        # Calculate circumcenter and circumradius squared
        C = circumcenter(pos1, pos2, pos3)
        R_sq = norm2(pos1 - C) # Radius squared from first point of triplet to center

        for particle_idx_in_triplet in triplet_indices
            # only consider particles that are within the original particle count
            # This ensures we only check particles that are part of the original set
            if particle_idx_in_triplet <= parameters.N 
                for neighbor_pbc_idx in arrays.neighborlist.voronoi_neighbors[particle_idx_in_triplet]
                    if neighbor_pbc_idx in triplet_indices
                        continue # Skip if the neighbor is one of the triplet vertices
                    end

                    # Fetch the position of the test particle
                    p_test = arrays.neighborlist.positions_with_pbc[neighbor_pbc_idx]
                    # Calculate squared distance to circumcenter
                    d_sq = norm2(p_test - C)
                    # Check for Delaunay violation
                    if d_sq < R_sq - epsilon
                        # A violation is found, return false
                        return false
                    end
                end
            end
        end
    end
    update_voronoi_vertices!(parameters, arrays, output, arrays.neighborlist.delaunay_facet_triplets)
    return true
end




"""

"""
function  update_voronoi_vertices!(parameters, arrays, output, triplets)

    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    N_pbc = length(positions_with_pbc)

    new_triplets = Tuple{Int, Int, Int}[]
    
    empty!(arrays.neighborlist.voronoi_vertices)
    empty!(arrays.neighborlist.voronoi_neighbors)
    empty!(arrays.neighborlist.voronoi_vertex_positions_per_particle)
    empty!(arrays.neighborlist.voronoi_vertex_indices)
    for _ in 1:N_pbc
        push!(arrays.neighborlist.voronoi_neighbors, Int[])
        push!(arrays.neighborlist.voronoi_vertex_positions_per_particle, SVector{2, Float64}[])
        push!(arrays.neighborlist.voronoi_vertex_indices, Int[])
    end

    for facet in triplets
        i = facet[1] # This is idx1_pbc
        j = facet[2] # This is idx2_pbc
        k = facet[3] # This is idx3_pbc

        
        triplet = (i,j,k)
        
        push!(new_triplets, triplet)
    
        # add these to the voronoi neighborlist for every particle pair, checking if it is already filled
        if !(j in arrays.neighborlist.voronoi_neighbors[i])
            push!(arrays.neighborlist.voronoi_neighbors[i], j)
        end
        if !(k in arrays.neighborlist.voronoi_neighbors[i])
            push!(arrays.neighborlist.voronoi_neighbors[i], k)
        end
        if !(i in arrays.neighborlist.voronoi_neighbors[j])
            push!(arrays.neighborlist.voronoi_neighbors[j], i)
        end
        if !(k in arrays.neighborlist.voronoi_neighbors[j])
            push!(arrays.neighborlist.voronoi_neighbors[j], k)
        end
        if !(i in arrays.neighborlist.voronoi_neighbors[k])
            push!(arrays.neighborlist.voronoi_neighbors[k], i)
        end
        if !(j in arrays.neighborlist.voronoi_neighbors[k])
            push!(arrays.neighborlist.voronoi_neighbors[k], j)
        end

        # compute the voronoi vertices as the circumcenter of the facet
        # and add it to the voronoi vertices list

        voronoi_vertex_position = circumcenter(
            positions_with_pbc[i],
            positions_with_pbc[j],
            positions_with_pbc[k]
        )

        push!(arrays.neighborlist.voronoi_vertices, voronoi_vertex_position)
        # # add the voronoi vertex to the voronoi vertices list

        push!(arrays.neighborlist.voronoi_vertex_indices[i], length(arrays.neighborlist.voronoi_vertices))
        push!(arrays.neighborlist.voronoi_vertex_positions_per_particle[i], voronoi_vertex_position)
        push!(arrays.neighborlist.voronoi_vertex_indices[j], length(arrays.neighborlist.voronoi_vertices))
        push!(arrays.neighborlist.voronoi_vertex_positions_per_particle[j], voronoi_vertex_position)
        push!(arrays.neighborlist.voronoi_vertex_indices[k], length(arrays.neighborlist.voronoi_vertices))
        push!(arrays.neighborlist.voronoi_vertex_positions_per_particle[k], voronoi_vertex_position)
    end

    for particle in 1:N_pbc
        voronoi_center = positions_with_pbc[particle]
        # sort the voronoi vertex indices counterclockwise
        voronoi_vertex_indices_new, voronoi_vertex_positions_per_particle_new = sort_indices_counter_clockwise(arrays.neighborlist.voronoi_vertex_indices[particle], arrays.neighborlist.voronoi_vertex_positions_per_particle[particle], arrays.neighborlist.voronoi_vertices, voronoi_center)
        # replace the voronoi vertex indices with the new ones
        arrays.neighborlist.voronoi_vertex_indices[particle] = voronoi_vertex_indices_new
        arrays.neighborlist.voronoi_vertex_positions_per_particle[particle] = voronoi_vertex_positions_per_particle_new
    end
    arrays.neighborlist.delaunay_facet_triplets = new_triplets
end