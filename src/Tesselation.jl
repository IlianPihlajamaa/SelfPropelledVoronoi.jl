
"""
This module is responsible for performing the Voronoi tesselation of the system.
Key functionalities include:
- Updating particle positions to account for periodic boundary conditions before tesselation (`update_positions_with_pbcs!`).
- Performing the Voronoi tesselation itself, typically by first computing the Delaunay triangulation and then deriving Voronoi cells (`voronoi_tesselation!`). This involves:
    - Identifying Voronoi neighbors for each particle.
    - Determining the vertices of each Voronoi cell.
    - Calculating geometric properties like circumcenters of Delaunay triangles, which correspond to Voronoi vertices (`circumcenter`).
- Sorting Voronoi vertices for consistent cell representation (`sort_indices_counter_clockwise`).
- Helper functions for managing and verifying the tesselation data (`verify_tesselation`, `update_delauney_vertices!`).
"""
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
    - `N::Int`: The number of original particles.
    - `box.box_sizes::SVector{2, Float64}`: Dimensions (Lx, Ly) of the simulation box.
    - `periodic_boundary_layer_depth::Float64`: The depth of the layer for including periodic images. **Note:** Currently, this parameter is overridden by a hardcoded value of `2.5` within the function.
- `arrays::ArrayStruct`: The struct holding simulation arrays.
    - `arrays.positions::Vector{SVector{2, Float64}}`: The original positions of all particles.
    - `arrays.neighborlist.positions_with_pbc::Vector{SVector{2, Float64}}`: This field is cleared and then populated in-place with the original positions and their relevant periodic images.
    - `arrays.neighborlist.position_indices::Vector{Int64}`: This field is cleared and then populated in-place with the original particle index corresponding to each position in `positions_with_pbc`.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures.

# Notes
- The fields `arrays.neighborlist.positions_with_pbc` and `arrays.neighborlist.position_indices` are modified in-place.
- The `pbc_layer_depth` used for determining which images to include is currently hardcoded to `2.5` inside the function, not taken from `parameters.periodic_boundary_layer_depth`.
"""
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
- `parameters::ParameterStruct`: The main simulation parameter struct, providing:
    - `N::Int`: The number of original (non-periodic image) particles.
    - `box.box_sizes::SVector{2, Float64}`: Dimensions (Lx, Ly) of the simulation box, used by `update_positions_with_pbcs!` and `sort_indices_counter_clockwise`.
- `arrays::ArrayStruct`: The struct holding simulation arrays.
    - `arrays.positions::Vector{SVector{2, Float64}}`: Used by `update_positions_with_pbcs!` (indirectly) and by `sort_indices_counter_clockwise` as the center for sorting vertices of original particles.
    - `arrays.neighborlist`: This substructure is extensively modified in-place. Its fields (`voronoi_vertices`, `voronoi_neighbors`, `voronoi_vertex_indices`, `voronoi_vertex_positions_per_particle`, `cell_centers_that_share_a_vertex`, `positions_with_pbc`, `position_indices`) are populated with the results of the tessellation.
- `output::Output`: The simulation output struct. Passed to `update_positions_with_pbcs!`.

# Notes
- This function heavily modifies the `arrays.neighborlist` structure in-place.
- It relies on `Quickhull.jl` for the Delaunay triangulation and helper functions like `update_positions_with_pbcs!`, `circumcenter`, and `sort_indices_counter_clockwise`.
"""
function voronoi_tesselation!(parameters, arrays, output)

    N = parameters.N
    Lx, Ly = parameters.box.box_sizes
    update_positions_with_pbcs!(parameters, arrays, output)

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
        i = facet[1]
        j = facet[2]
        k = facet[3]

    
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

# Notes
- The arguments `Lx` and `Ly` (box sizes) are included in the function signature but are not currently used in the angle calculation. This implies the sorting does not explicitly handle periodic boundary conditions for angle calculation, which might be relevant if vertices cross periodic boundaries in a way that simple `atan` would misorder. However, for typical Voronoi cells generated from Delaunay triangulations (especially with periodic images handled by `update_positions_with_pbcs!`), the vertices should generally be within a local region where direct angle calculation is sufficient.
"""
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

"""
    outer(a::SVector{2, T}, b::SVector{2, T}) where T

Computes a form of 2D vector outer product or cross product related quantity.
Given two 2D vectors `a = [a1, a2]` and `b = [b1, b2]`, the scalar component
of the 2D cross product (often referred to as the "perp dot product" or
z-component of the 3D cross product if vectors are embedded in the xy-plane)
is `a1*b2 - a2*b1`.

This function returns a 2-component `SVector{2, T}` where:
- The first component is `a[1]*b[2] - a[2]*b[1]`.
- The second component is `a[2]*b[1] - a[1]*b[2]`, which is the negative of the first component.

The specific reason for returning this two-component vector with one being the
negative of the other, rather than just the scalar value, is not immediately obvious
from the function itself and might relate to specific downstream calculations or conventions
within the broader simulation code.

# Arguments
- `a::SVector{2, T}`: The first 2D input vector. `T` is its element type.
- `b::SVector{2, T}`: The second 2D input vector. `T` is its element type.

# Returns
- `SVector{2, T}`: A 2D static vector where the first component is `a[1]*b[2] - a[2]*b[1]` and the second component is `a[2]*b[1] - a[1]*b[2]`.
"""
function outer(a::SVector{2, T}, b::SVector{2, T}) where T
    return SVector{2, T}(a[1]*b[2] - a[2]*b[1], a[2]*b[1] - a[1]*b[2])
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

Checks if the current Voronoi tessellation stored in `arrays.neighborlist` is still
considered valid. This function is intended to be used to determine if a full
re-tessellation is necessary, potentially saving computational cost if the
tessellation can be assumed to be stable under small particle displacements.

**Note:** The current implementation is a placeholder and always returns `false`.
This means that, with the current logic, the system will always perform a full
re-tessellation in every step where this check is made (e.g., in `compute_forces_SPV!`).
A more sophisticated implementation might check, for example, if any particle has
moved more than a certain threshold distance since the last tessellation.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. (Currently unused in this placeholder version).
- `arrays::ArrayStruct`: The struct holding simulation arrays, including `arrays.neighborlist` which contains the current tessellation data. (Currently unused in this placeholder version).
- `output::Output`: The simulation output struct. (Currently unused in this placeholder version).

# Returns
- `::Bool`: Currently, this function always returns `false`, indicating that the tessellation is always considered invalid and needs recomputation.
"""
function verify_tesselation(parameters, arrays, output)
    return false
end

"""
    update_delauney_vertices!(parameters, arrays, output)

Intended to update the Delaunay vertices (which correspond to Voronoi vertices)
incrementally, without performing a full re-tessellation. This could be useful
if only a few particles have moved slightly, and the overall topology of the
Delaunay triangulation is expected to remain largely unchanged.

**Note:** The current implementation is a placeholder and is not functional.
It immediately calls `error()`, indicating that this feature is not yet implemented.
A functional version would likely involve algorithms for local updates to
Delaunay triangulations.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. (Currently unused as the function is a placeholder).
- `arrays::ArrayStruct`: The struct holding simulation arrays, including `arrays.neighborlist` which would be modified if this function were implemented. (Currently unused).
- `output::Output`: The simulation output struct. (Currently unused).

# Throws
- `Error`: This function, in its current state, always throws an error, indicating it is not implemented.
"""
function update_delauney_vertices!(parameters, arrays, output)
    error()
    return
end
