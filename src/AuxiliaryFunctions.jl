
"""
apply_periodic_boundary_conditions(position, box_sizes)

Applies periodic boundary conditions to a given position, ensuring it wraps around the simulation box.
For a position coordinate `x` and a box dimension `L`, the new coordinate `x'` is `x - floor(x/L) * L`.
This maps `x` to the interval `[0, L)`. The same logic applies to all dimensions.

# Arguments
- `position`: The original position, typically an `SVector` or any `AbstractVector` representing coordinates (e.g., `[x, y]`).
- `box_sizes`: The dimensions of the simulation box, typically an `SVector` or `AbstractVector` (e.g., `[Lx, Ly]`).

# Returns
- `new_position`: The position after applying periodic boundary conditions, of the same type as `position`.
"""
function apply_periodic_boundary_conditions(position, box_sizes)
    new_position = position .- floor.(position ./ box_sizes) .* box_sizes
    return new_position
end

"""
    compute_pair_distance_vector(p1, p2, box_sizes)

Computes the shortest vector (displacement) from particle `p1` to particle `p2`
in a system with periodic boundary conditions (PBC). This is often referred to as
the minimum image convention. For a coordinate difference `dx` and box dimension `L`,
the PBC-aware difference `dx'` is `dx - round(dx/L) * L`.

# Arguments
- `p1`: Position of the first particle (e.g., an `SVector` or `AbstractVector` like `[x1, y1]`).
- `p2`: Position of the second particle (e.g., an `SVector` or `AbstractVector` like `[x2, y2]`).
- `box_sizes`: Dimensions of the simulation box (e.g., an `SVector` or `AbstractVector` like `[Lx, Ly]`).

# Returns
- `delta`: The displacement vector from `p1` to `p2` after applying the minimum image convention due to PBC. This will be of the same type as `p1` and `p2`.
"""
function compute_pair_distance_vector(p1, p2, box_sizes)
    delta = p2 - p1
    delta = delta - round.(delta ./ box_sizes) .* box_sizes
    return delta
end

"""
    update_perimeters!(parameters, arrays, output)

Calculates and updates the perimeters of all Voronoi cells.
The perimeter of each cell is computed by summing the lengths of the segments
connecting its Voronoi vertices in sequence. The results are stored in `arrays.perimeters`.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. Used here to get `N`, the number of particles.
- `arrays::ArrayStruct`: The struct holding simulation arrays. 
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures across the module.

# Notes
- This function assumes that the Voronoi vertex positions for each cell are already calculated and available in `arrays.neighborlist.voronoi_vertices`.
- The `arrays.perimeters` vector is updated in-place with the new perimeter values.
"""
function update_perimeters!(parameters, arrays, output)
    for i in 1:parameters.N
        arrays.perimeters[i] = compute_perimeter_i(i, parameters, arrays, output)
    end
end

"""
    update_areas!(parameters, arrays, output)

Calculates and updates the areas of all Voronoi cells.
The area of each cell is computed using the shoelace formula (also known as Gauss's area formula),
based on the coordinates of its Voronoi vertices. The results are stored in `arrays.areas`.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. Used here to get `N`, the number of particles.
- `arrays::ArrayStruct`: The struct holding simulation arrays.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures across the module.

# Notes
- This function assumes that the Voronoi vertex positions for each cell are already calculated and available in `arrays.neighborlist.voronoi_vertices`.
- The `arrays.areas` vector is updated in-place with the new area values. The shoelace formula calculates signed area, so the division by 2.0 yields the geometric area assuming a consistent vertex ordering.
"""
function update_areas!(parameters, arrays, output)
    for i in 1:parameters.N
        area_i = compute_area_i(i, parameters, arrays, output)
        arrays.areas[i] = area_i
    end
end


"""
    compute_energy(parameters, arrays, output)

Computes the total potential energy of the system. This energy is typically derived from
the deviations of individual Voronoi cell areas and perimeters from their target values,
penalized by corresponding spring constants.

This function first calls `update_areas!` and `update_perimeters!` to ensure that
the current areas and perimeters in `arrays` are up-to-date before calculating the energy.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. 
- `arrays::ArrayStruct`: The struct holding simulation arrays. It's passed to `update_areas!` and `update_perimeters!`, which use the current areas and perimeters.
- `output::Output`: The simulation output struct. Passed to `update_areas!` and `update_perimeters!`.

# Returns
- `potential_energy::Float64`: The total calculated potential energy of the system. This is the sum of energy contributions from each cell, where each cell's energy is `K_A * (area - target_area)^2 + K_P * (perimeter - target_perimeter)^2`.
"""
function compute_energy(parameters, arrays, output)
    N = parameters.N
    update_areas!(parameters, arrays, output)
    update_perimeters!(parameters, arrays, output)
    E = 0.0
    for i in 1:N
        E += compute_energy_i(i, parameters, arrays, output)
    end
    output.potential_energy = E
    return E
end


"""    
    compute_energy_i(i, parameters, arrays, output)

Computes the potential energy contribution of a single Voronoi cell indexed by `i`.
This function calculates the energy based on the deviation of the cell's area and perimeter
from their target values, penalized by the corresponding spring constants.
# Arguments
- `i::Int`: The index of the Voronoi cell for which to compute the energy.
- `parameters::ParameterStruct`: The main simulation parameter struct. Used to access target areas, target perimeters, and spring constants.
- `arrays::ArrayStruct`: The struct holding simulation arrays. Used to access the current areas and perimeters of the Voronoi cells.
- `output::Output`: The simulation output struct. Not directly used in this function but included for consistency in function signatures across the module.

# Returns
- `E::Float64`: The computed energy for the Voronoi cell indexed by `i`. This is calculated as:
K_A[i] * (areas[i] - target_areas[i])^2 + K_P[i] * (perimeters[i] - target_perimeters[i])^2
""" 
function compute_energy_i(i, parameters, arrays, output)
    areas = arrays.areas
    perimeters = arrays.perimeters
    target_perimeters = parameters.particles.target_perimeters
    target_areas = parameters.particles.target_areas
    K_P = parameters.particles.K_P
    K_A = parameters.particles.K_A
    E = K_A[i]*(areas[i] - target_areas[i])^2 + K_P[i]*(perimeters[i] - target_perimeters[i])^2
    return E
end


function compute_area_i(i, parameters, arrays, output)
    vor_indices = arrays.neighborlist.voronoi_vertex_indices[i]
    N_vertices_i = arrays.neighborlist.N_voronoi_vertices_pp[i]
    vor_positions = arrays.neighborlist.voronoi_vertices
    # compute the area of the voronoi cell
    area = 0.0
    for j in 1:N_vertices_i
        posj = vor_positions[vor_indices[j]]
        l = j+1
        if l > N_vertices_i
            l = 1
        end
        posl = vor_positions[vor_indices[l]]
        area += (posj[1]*posl[2] - posl[1]*posj[2])/2
    end

    return area
end


function compute_perimeter_i(i, parameters, arrays, output)
    vor_indices = arrays.neighborlist.voronoi_vertex_indices[i]
    N_vertices_i = arrays.neighborlist.N_voronoi_vertices_pp[i]
    vor_positions = arrays.neighborlist.voronoi_vertices
    # compute the perimeter of the voronoi cell
    perimeter = 0.0
    for j in 1:N_vertices_i
        posj = vor_positions[vor_indices[j]]
        l = j+1
        if l > N_vertices_i
            l = 1
        end
        posl = vor_positions[vor_indices[l]]
        # compute the distance between posj and posl
        dr = posl - posj
        perimeter += sqrt(sum(dr.*dr))
    end

    return perimeter
end


"""
    compute_msd(displacement_array)

returns [1,2,3,4,5,6,7,8,9,10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, etc]
    Computes the mean squared displacement (MSD) of particles based on their displacement vectors.
    The MSD is a measure of the average squared distance that particles have moved from their initial positions.
"""
function when_to_print_array(max_steps)
    print_array = Int[]
    step = 1
    val = 0
    while val <= max_steps
        push!(print_array, val)
        if val/step == 10
            step *= 10
        end 
        val += step
    end
    return print_array
end

function compute_eta(parameters, arrays, output, start_time, total_steps)
    elapsed_time = time() - start_time
    remaining_steps = total_steps - output.steps_done
    eta = remaining_steps * (elapsed_time / output.steps_done)
    return eta
end