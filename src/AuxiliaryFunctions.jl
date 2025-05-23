function apply_periodic_boundary_conditions(position, box_sizes)
    new_position = position - floor.(position ./ box_sizes) .* box_sizes
    return new_position
end


function compute_pair_distance_vector(p1, p2, box_sizes)
    # Compute the pair distance between two particles with periodic boundary conditions
    delta = p2 - p1
    delta = delta - round.(delta ./ box_sizes) .* box_sizes
    return delta
end


function update_perimeters!(parameters, arrays, output)
    voronoi_vertex_positions_per_particle = arrays.neighborlist.voronoi_vertex_positions_per_particle
    for particle_i in 1:parameters.N
        vor_positions = voronoi_vertex_positions_per_particle[particle_i]
        # compute the perimeter of the voronoi cell
        perimeter = 0.0
        for j in 1:length(vor_positions)
            posj = vor_positions[j]
            l = j+1
            if l > length(vor_positions)
                l = 1
            end
            posl = vor_positions[l]
            # compute the distance between posj and posl
            dr = posl - posj
            perimeter += sqrt(sum(dr.*dr))
        end
        arrays.perimeters[particle_i] = perimeter
    end
end

function update_areas!(parameters, arrays, output)
    voronoi_vertex_positions_per_particle = arrays.neighborlist.voronoi_vertex_positions_per_particle
    for particle_i in 1:parameters.N
        vor_positions = voronoi_vertex_positions_per_particle[particle_i]
        # compute the area of the voronoi cell
        area = 0.0
        for j in 1:length(vor_positions)
            posj = vor_positions[j]
            l = j+1
            if l > length(vor_positions)
                l = 1
            end
            posl = vor_positions[l]
            area += (posj[1]*posl[2] - posl[1]*posj[2])
        end
        arrays.areas[particle_i] = area / 2.0
    end
end


function compute_energy(parameters, arrays, output)
    # compute the potential energy of the system

    update_areas!(parameters, arrays, output)
    update_perimeters!(parameters, arrays, output)

    potential_energy = 0.0

    target_perimeters = parameters.particles.target_perimeters
    target_areas = parameters.particles.target_areas
    K_P = parameters.particles.K_P
    K_A = parameters.particles.K_A

    for particle in 1:parameters.N
        area = arrays.areas[particle]
        perimeter = arrays.perimeters[particle]
        # compute the potential energy of the voronoi cell
        potential_energy += K_A[particle] * (area - target_areas[particle])^2
        potential_energy += K_P[particle] * (perimeter - target_perimeters[particle])^2
    end
    return potential_energy
end


