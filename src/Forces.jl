"""
This module is intended to define various force laws that can be used in the
SelfPropelledVoronoi simulation. These could include inter-particle interaction
forces, forces derived from external potentials, or other custom forces relevant
to the specific system being modeled.

Currently, it may primarily contain placeholder or example force functions,
such as `test_force`.
"""
"""
    test_force(r)

An example or placeholder function that returns a force proportional to a given input `r`.
This function simply multiplies the input by two.

It is likely intended for testing purposes or as a template for more complex force calculations.

# Arguments
- `r`: An input value, typically a scalar (e.g., a distance or a generic scalar variable).

# Returns
- `2r::Number`: Twice the input value `r`.
"""
function test_force(r)
    return 2r
end
