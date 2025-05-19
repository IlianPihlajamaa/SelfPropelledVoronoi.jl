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