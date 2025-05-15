function apply_periodic_boundary_conditions(position, box_sizes)
    new_position = position - floor.(position ./ box_sizes) .* box_sizes
    return new_position
end

