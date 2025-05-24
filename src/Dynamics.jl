"""
This module implements the time evolution of the simulated system.
It is responsible for:
- Calculating forces acting on particles: This includes forces derived from the Voronoi tessellation (e.g., area and perimeter elasticity) via `compute_forces_SPV!`, and potentially other interaction models like the Gaussian Core Model (`compute_forces_GCM!`).
- Integrating the equations of motion: It uses numerical integration schemes to update particle positions and orientations over time. Different integrators like Euler-Murayama (`do_time_step_Euler_Murayama`) and Euler-Heun (`do_time_step_Euler_Heun!`) are implemented.
- Managing the main simulation loop: The `run_simulation!` function orchestrates the simulation steps, including force computation, integration, and data dumping.
"""
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
    if !verify_tesselation(parameters, arrays, output)
        voronoi_tesselation!(parameters, arrays, output)
    else
        update_delauney_vertices!(parameters, arrays, output)
    end
    N = parameters.N
    voronoi_neighbors = arrays.neighborlist.voronoi_neighbors
    voronoi_vertices = arrays.neighborlist.voronoi_vertices
    voronoi_vertices_per_particle = arrays.neighborlist.voronoi_vertex_positions_per_particle
    voronoi_indices = arrays.neighborlist.voronoi_vertex_indices
    update_perimeters!(parameters, arrays, output)
    update_areas!(parameters, arrays, output)

    areas = arrays.areas
    perimeters = arrays.perimeters
    target_perimeters = parameters.particles.target_perimeters
    target_areas = parameters.particles.target_areas
    K_P = parameters.particles.K_P
    K_A = parameters.particles.K_A

    positions_with_pbc = arrays.neighborlist.positions_with_pbc
    common_vertices = zeros(Int, 2)


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
            h_i_j = voronoi_vertices_per_particle[particle_i][j]
            next = j % length(voronoi_neighbors[particle_i]) + 1
            h_i_jp1 = voronoi_vertices_per_particle[particle_i][next]
            prev = j - 1 == 0 ? length(voronoi_neighbors[particle_i]) : j - 1
            h_i_jm1 = voronoi_vertices_per_particle[particle_i][prev]

            # compute the distance between h_i_j and h_i_jp1
            hijp1_m_hij = h_i_jp1 - h_i_j
            hij_m_hijm1 = h_i_j - h_i_jm1
            hijp1_m_hij_length = sqrt(hijp1_m_hij[1]^2 + hijp1_m_hij[2]^2)
            hij_m_hijm1_length = sqrt(hij_m_hijm1[1]^2 + hij_m_hijm1[2]^2)

            dAi_dhijx = (h_i_jp1[2] - h_i_jm1[2])/2
            dAi_dhijy = (h_i_jm1[1] - h_i_jp1[1])/2
            dPi_dhij = hij_m_hijm1/hij_m_hijm1_length + hijp1_m_hij/hijp1_m_hij_length
            dPi_dhijx = dPi_dhij[1]
            dPi_dhijy = dPi_dhij[2]

            i_a, i_b, i_c = arrays.neighborlist.cell_centers_that_share_a_vertex[j]

            if i_a == particle_i
                particle_l = i_b
                particle_k = i_c
            elseif i_b == particle_i
                particle_l = i_a
                particle_k = i_c
            else
                particle_l = i_a
                particle_k = i_b
            end

            rl = positions_with_pbc[particle_l]
            rk = positions_with_pbc[particle_k]
            ri = positions_with_pbc[particle_i]
            xi = ri[1]
            yi = ri[2]
            xl = rl[1]
            yl = rl[2]
            xk = rk[1]
            yk = rk[2]

            # from mathematica:
            dhijx_dxi = -(((yk - yl)*(2*xi*(xl*(-yi + yk) + xk*(yi - yl)) + xk^2*(-yi + yl) + xi^2*(-yk + yl) + (yi - yk)*(xl^2 + yi*yk - (yi + yk)*yl + yl^2)))/ (2*(xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2))
            dhijx_dyi = ((yk - yl)*(xi^2*(xk - xl) + xk^2*xl + xl*(yi - yk)^2 - xk*(xl^2 + (yi - yl)^2) + xi*(-xk^2 + xl^2 + 2*yi*yk - yk^2 - 2*yi*yl + yl^2)))/ (2*(xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2)
            dhijy_dxi = -1/2*((xk - xl)*(xk^2*(yi - yl) + xi^2*(yk - yl) - (yi - yk)*(xl^2 + yi*yk - (yi + yk)*yl + yl^2) + 2*xi*(xl*(yi - yk) + xk*(-yi + yl))))/ (xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2
            dhijy_dyi = -1/2*((xk - xl)*(xi^2*(xk - xl) + xk^2*xl + xl*(yi - yk)^2 - xk*(xl^2 + (yi - yl)^2) + xi*(-xk^2 + xl^2 + 2*yi*yk - yk^2 - 2*yi*yl + yl^2)))/ (xl*(-yi + yk) + xk*(yi - yl) + xi*(-yk + yl))^2


            dAi_dxi += dAi_dhijx * dhijx_dxi + dAi_dhijy * dhijy_dxi
            dAi_dyi += dAi_dhijx * dhijx_dyi + dAi_dhijy * dhijy_dyi
            dPi_dxi += dPi_dhijx * dhijx_dxi + dPi_dhijy * dhijy_dxi
            dPi_dyi += dPi_dhijx * dhijx_dyi + dPi_dhijy * dhijy_dyi   
            # @show dAi_dxi, dAi_dyi, dPi_dxi, dPi_dyi
        end

        dEi_dxi = dEi_dAi*dAi_dxi + dEi_dPi*dPi_dxi
        dEi_dyi = dEi_dAi*dAi_dyi + dEi_dPi*dPi_dyi
        # @show dEi_dxi, dEi_dyi
        Fx -= dEi_dxi
        Fy -= dEi_dyi

        # part 2: -dEj/dri for all j
        for particle_j in voronoi_neighbors[particle_i]
            originalj = arrays.neighborlist.position_indices[particle_j]
            dEj_dAj = 2 * K_A[originalj] * (areas[originalj] - target_areas[originalj])
            dEj_dPj = 2 * K_P[originalj] * (perimeters[originalj] - target_perimeters[originalj])


            #find the two common vertices h and g that occur both in the voronoi_vertices of i and of j
            Nfound = 0
            for vertex_k in voronoi_indices[particle_j]
                if vertex_k in voronoi_indices[particle_i]
                    Nfound += 1
                    common_vertices[Nfound] = vertex_k                    
                end
            end
            if Nfound != 2
                error("Error: neighborlist does not contain two common voronoi vertices")
            end

            h = voronoi_vertices[common_vertices[1]]
            g = voronoi_vertices[common_vertices[2]]

            # h is the first common vertex, g is the second common vertex by counterclockwise order around j

            dAj_dhx = 0.0
            dAj_dhy = 0.0
            dAj_dgx = 0.0
            dAj_dgy = 0.0
            dPj_dhx = 0.0
            dPj_dhy = 0.0
            dPj_dgx = 0.0
            dPj_dgy = 0.0

            dhx_dxi = 0.0
            dhx_dyi = 0.0
            dhy_dxi = 0.0
            dhy_dyi = 0.0
            dgx_dxi = 0.0
            dgx_dyi = 0.0
            dgy_dxi = 0.0
            dgy_dyi = 0.0
            
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
            
        end 

        Fi = SVector(Fx , Fy)   
        arrays.forces[particle_i] = Fi

    end
    # error()
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
            if r < 2^(1/6)
                # Gaussian core
                f = -2 * rij  * exp(-r^2)
                Fi += f
            end
        end
        arrays.forces[i] =  Fi
    end
    return
end

"""
    do_time_step(parameters, arrays, output)

Advances the simulation by a single time step by calling an appropriate integration scheme.

This function acts as a dispatcher for different integration methods.
Currently, it selects:
- `do_time_step_Euler_Murayama` for the very first step of the simulation (`output.steps_done == 0`).
- `do_time_step_Euler_Heun!` for all subsequent steps.

This strategy might be used, for example, to ensure proper initialization of "old" states
required by more complex integrators like Heun's method.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct, containing settings like `dt`, particle properties, etc. Passed through to the chosen integrator function.
- `arrays::ArrayStruct`: The struct holding simulation arrays (positions, forces, orientations, etc.). These are modified in-place by the chosen integrator function.
- `output::Output`: The simulation output struct. `output.steps_done` is used to determine which integration scheme to use. It's also passed through to the chosen integrator function.

# Returns
- The function itself doesn't return a value but calls one of the integrator functions (`do_time_step_Euler_Murayama` or `do_time_step_Euler_Heun!`), which in turn modify `arrays` (and potentially `output`) in-place.
"""
function do_time_step(parameters, arrays, output)
    if output.steps_done == 0
        do_time_step_Euler_Murayama(parameters, arrays, output)
    else
        do_time_step_Euler_Heun!(parameters, arrays, output)
    end
end

"""
    do_time_step_Euler_Murayama(parameters, arrays, output)

Performs a single time step using the Euler-Maruyama integration scheme, suitable for
stochastic differential equations (SDEs). This method updates particle positions and
orientations considering deterministic forces, active self-propulsion forces, and
stochastic rotational diffusion.

The function first calls `compute_forces_SPV!` to calculate the current deterministic
forces on all particles based on their positions and Voronoi cell properties.
Then, it updates orientations by adding a random component scaled by the rotational
diffusion constant. Finally, it updates positions based on the sum of deterministic
forces and active forces (directed by the orientation at the *beginning* of the step),
scaled by mobility and time step `dt`. Periodic boundary conditions are applied to
the new positions.

At the end of the step, the `old_positions`, `old_forces`, and `old_orientations`
arrays are updated to store the state for potential use in subsequent, more complex
integration schemes (though this specific integrator is often used as a simpler baseline
or for the first step).

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct, providing:
    - `dt::Float64`: The time step size.
    - `frictionconstant::Float64`: Used to calculate mobility (`1 / frictionconstant`).
    - `particles.active_force_strengths::Vector{Float64}`: Magnitudes of the active force for each particle.
    - `particles.rotational_diffusion_constants::Vector{Float64}`: Rotational diffusion rates for each particle.
    - `N::Int`: The total number of particles.
    - `box.box_sizes::SVector{2, Float64}`: Dimensions of the simulation box for periodic boundary conditions.
- `arrays::ArrayStruct`: The struct holding simulation arrays. The following fields are used and modified:
    - `positions::Vector{SVector{2, Float64}}`: Updated in-place with new particle positions.
    - `forces::Vector{SVector{2, Float64}}`: Updated in-place by `compute_forces_SPV!`.
    - `orientations::Vector{Float64}`: Updated in-place with new particle orientations.
    - `old_positions::Vector{SVector{2, Float64}}`: Updated in-place to store the new `positions`.
    - `old_forces::Vector{SVector{2, Float64}}`: Updated in-place to store the new `forces`.
    - `old_orientations::Vector{Float64}`: Updated in-place to store the new `orientations`. (Note: the position update uses orientations from *before* this step's rotational diffusion for the active force direction).
- `output::Output`: The simulation output struct. Passed to `compute_forces_SPV!`.

# Notes
- All relevant fields in `arrays` (positions, forces, orientations, old_positions, old_forces, old_orientations) are modified in-place.
"""
function do_time_step_Euler_Murayama(parameters, arrays, output)
    dt = parameters.dt
    mobility = 1 / parameters.frictionconstant
    active_force_strengths = parameters.particles.active_force_strengths
    rotational_diffusion_constants = parameters.particles.rotational_diffusion_constants

    positions = arrays.positions
    forces = arrays.forces
    orientations = arrays.orientations
    old_positions = arrays.old_positions
    old_forces = arrays.old_forces
    old_orientations = arrays.old_orientations
    box_sizes = parameters.box.box_sizes
    # compute forces
    compute_forces_SPV!(parameters, arrays, output)
    # update orientations
    for particle in 1:parameters.N
        old_orientation = orientations[particle]
        Dr = rotational_diffusion_constants[particle]
        new_orientation = old_orientation + sqrt(2*dt*Dr) * randn()
        orientations[particle] = new_orientation
    end

    # update positions
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))
        new_position = positions[particle] + dt * (forces[particle] * mobility + active_force_strengths[particle] * old_orientation_vector)
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end

    # update orientations
    old_forces .= forces
    old_positions .= positions
    old_orientations .= orientations

    return
end

"""
    do_time_step_Euler_Heun!(parameters, arrays, output)

Performs a single time step using the Euler-Heun predictor-corrector integration scheme.
This method is suitable for stochastic differential equations (SDEs) and offers better
stability and accuracy than the simpler Euler-Maruyama method for some systems.
It involves predicting a future state and then correcting it using an average of forces/drifts.

The operational steps are:
1.  **Compute forces at current positions**: `compute_forces_SPV!` is called to get forces `F(x(t_n))` based on current `arrays.positions` (state at `t_n`).
2.  **Update orientations**: Particle orientations `θ(t_n)` are updated to `θ(t_{n+1})` by adding a stochastic rotational diffusion term. `arrays.orientations` now stores `θ(t_{n+1})`. The orientations from the start of the step, `θ(t_n)`, are available in `arrays.old_orientations` (assuming they were stored from a previous step or initialization).
3.  **Predictor step**:
    -   A predicted position `x_tilde(t_{n+1})` is calculated using forces `F(x(t_n))` and active forces based on orientations `θ(t_n)` (from `arrays.old_orientations`).
    -   `arrays.positions` are updated to these predicted positions.
4.  **Store current forces**: The forces `F(x(t_n))` (used in the predictor step) are stored in `arrays.old_forces`.
5.  **Compute forces at predicted positions**: `compute_forces_SPV!` is called again to calculate forces `F(x_tilde(t_{n+1}))` based on the predicted `arrays.positions`. `arrays.forces` now stores these new forces.
6.  **Corrector step**:
    -   The final positions `x(t_{n+1})` are calculated using an average of the "old" forces/orientations and the "new" (predicted) forces/orientations.
        Specifically, it uses `0.5 * [ (F(x(t_n)) + v_active(θ(t_n))) + (F(x_tilde(t_{n+1})) + v_active(θ(t_{n+1}))) ]`.
    -   `arrays.positions` are updated to these corrected positions.
7.  **Update "old" states**: `arrays.old_positions` are updated with the corrected `arrays.positions`, and `arrays.old_orientations` are updated with `arrays.orientations` (`θ(t_{n+1})`) to prepare for the next time step.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct, providing:
    - `dt::Float64`: The time step size.
    - `frictionconstant::Float64`: Used to calculate mobility (`1 / frictionconstant`).
    - `particles.active_force_strengths::Vector{Float64}`: Magnitudes of the active force.
    - `particles.rotational_diffusion_constants::Vector{Float64}`: Rotational diffusion rates.
    - `N::Int`: The total number of particles.
    - `box.box_sizes::SVector{2, Float64}`: Dimensions for periodic boundary conditions.
- `arrays::ArrayStruct`: The struct holding simulation arrays. The following are used and modified in-place:
    - `positions::Vector{SVector{2, Float64}}`: Initially `x(t_n)`, then `x_tilde(t_{n+1})`, then corrected `x(t_{n+1})`.
    - `forces::Vector{SVector{2, Float64}}`: Initially stores `F(x(t_n))`, then `F(x_tilde(t_{n+1}))`.
    - `orientations::Vector{Float64}`: Updated to `θ(t_{n+1})` after rotational diffusion.
    - `old_positions::Vector{SVector{2, Float64}}`: Stores `x(t_n)` at the start, updated to corrected `x(t_{n+1})` at the end.
    - `old_forces::Vector{SVector{2, Float64}}`: Updated to store `F(x(t_n))`.
    - `old_orientations::Vector{Float64}`: Stores `θ(t_n)` at the start, updated to `θ(t_{n+1})` at the end.
- `output::Output`: The simulation output struct. Passed to `compute_forces_SPV!`.

# Notes
- All relevant fields in `arrays` are modified in-place throughout the function.
- This scheme aims for better accuracy by averaging drifts over the interval.
"""
function do_time_step_Euler_Heun!(parameters, arrays, output)
    dt = parameters.dt
    mobility = 1 / parameters.frictionconstant
    active_force_strengths = parameters.particles.active_force_strengths
    rotational_diffusion_constants = parameters.particles.rotational_diffusion_constants

    positions = arrays.positions
    forces = arrays.forces
    orientations = arrays.orientations
    old_positions = arrays.old_positions
    old_forces = arrays.old_forces
    old_orientations = arrays.old_orientations

    box_sizes = parameters.box.box_sizes
    # compute forces
    compute_forces_SPV!(parameters, arrays, output)

    # update orientations
    for particle in 1:parameters.N
        old_orientation = orientations[particle]
        Dr = rotational_diffusion_constants[particle]
        new_orientation = old_orientation + sqrt(2*dt*Dr) * randn()
        orientations[particle] = new_orientation
    end

    # predictor step
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))


        # update positions
        new_position = positions[particle] + dt * (forces[particle] * mobility + active_force_strengths[particle] * old_orientation_vector)
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end


    old_forces .= forces
        # Check if the simulation should be retesselated

    compute_forces_SPV!(parameters, arrays, output)

    # corrector step
    for particle in 1:parameters.N
        old_orientation = old_orientations[particle]
        old_orientation_vector = SVector(cos(old_orientation), sin(old_orientation))
        new_orientation = orientations[particle]
        new_orientation_vector = SVector(cos(new_orientation), sin(new_orientation))

        old_position = old_positions[particle]
        old_force = old_forces[particle]
        new_force = forces[particle]
        
        # update positions
        new_position = old_position + dt/2 * (
            old_force * mobility + active_force_strengths[particle] * old_orientation_vector +
            new_force * mobility + active_force_strengths[particle] * new_orientation_vector
        )
        positions[particle] = apply_periodic_boundary_conditions(new_position, box_sizes)
    end


    # update orientations

    old_positions .= positions
    old_orientations .= orientations

    return
end

"""
    run_simulation!(parameters, arrays, output)

Runs the main simulation loop for a specified number of time steps.

Before starting the loop, it performs an initial Voronoi tessellation by calling
`voronoi_tesselation!(parameters, arrays, output)` to ensure the system's geometric
properties are up-to-date.

The main loop then proceeds as follows for each step:
1.  **Verbose Output**: If `parameters.verbose` is true and the current step is a multiple of 100, the step number is printed to the console.
2.  **Callback Invocation**: The user-provided `parameters.callback` function is called, passing `parameters`, `arrays`, and `output`. This allows for custom actions or analysis during the simulation.
3.  **Time Step Advancement**: `do_time_step(parameters, arrays, output)` is called to advance the simulation state by one time step, updating particle positions, orientations, and forces.
4.  **Data Dumping**: If `parameters.dump_info.save` is true, it checks if the current `step` is in `parameters.dump_info.when_to_save_array`. If so, `save_simulation_state!(parameters, arrays, output)` is called to save the current simulation state.
5.  **Step Counting**: The local `step` counter is incremented, and `output.steps_done` is updated to reflect the completed step.
6.  **Termination**: The loop breaks if `step` exceeds `parameters.N_steps`.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. Key fields used:
    - `N_steps::Int`: The total number of simulation steps to perform.
    - `verbose::Bool`: Flag to enable/disable verbose output.
    - `dump_info::DumpInfo`: Contains parameters for saving simulation data, including `save` (flag), `when_to_save_array` (steps for saving).
    - `callback`: A user-defined function called at each step.
- `arrays::ArrayStruct`: The struct holding simulation arrays (positions, forces, orientations, etc.). This is modified in-place by the functions called within the loop (e.g., `do_time_step`, `voronoi_tesselation!`).
- `output::Output`: The simulation output struct.
    - `output.steps_done`: Used to initialize the local `step` counter and is updated at each iteration. It is also passed to and modified by functions called within the loop.

# Notes
- The `arrays` and `output` structs are modified in-place throughout the simulation by the various functions called.
- The simulation starts from the step number indicated by `output.steps_done` at the beginning of the function call.
"""
function run_simulation!(parameters, arrays, output)
    step = output.steps_done

    voronoi_tesselation!(parameters, arrays, output)

    # Main simulation loop
    while true
        if parameters.verbose && step % 100 == 0
            println("Step: ", step)
        end
        # invoke callback
        parameters.callback(parameters, arrays, output)

        do_time_step(parameters, arrays, output)

        # Check if the simulation should be saved
        if parameters.dump_info.save
            if step in parameters.dump_info.when_to_save_array
                save_simulation_state!(parameters, arrays, output)
            end
        end


        # Check if the simulation should be stopped
        step = step + 1
        output.steps_done = step
        if step > parameters.N_steps
            break
        end
    end

end