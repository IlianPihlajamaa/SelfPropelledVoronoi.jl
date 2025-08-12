

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
    run_simulation!(parameters, arrays, output, N_steps)

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
6.  **Termination**: The loop breaks if `step` exceeds `N_steps`.

# Arguments
- `parameters::ParameterStruct`: The main simulation parameter struct. Key fields used:
    - `verbose::Bool`: Flag to enable/disable verbose output.
    - `dump_info::DumpInfo`: Contains parameters for saving simulation data, including `save` (flag), `when_to_save_array` (steps for saving).
    - `callback`: A user-defined function called at each step.
- `arrays::ArrayStruct`: The struct holding simulation arrays (positions, forces, orientations, etc.). This is modified in-place by the functions called within the loop (e.g., `do_time_step`, `voronoi_tesselation!`).
- `output::Output`: The simulation output struct.
    - `output.steps_done`: Used to initialize the local `step` counter and is updated at each iteration. It is also passed to and modified by functions called within the loop.
- `N_steps::Int`: The total number of time steps to run the simulation.

# Notes
- The `arrays` and `output` structs are modified in-place throughout the simulation by the various functions called.
- The simulation starts from the step number indicated by `output.steps_done` at the beginning of the function call.
"""
function run_simulation!(parameters, arrays, output, N_steps)
    step = output.steps_done
    start_step = step
    voronoi_tesselation!(parameters, arrays, output)
    print_arr = when_to_print_array(N_steps + start_step + 10)
    # Main simulation loop
    t0 = time()
    while true
        if parameters.verbose && step in print_arr
            eta = round(compute_eta(parameters, arrays, output, t0, N_steps), digits=2)
            Epot = round(output.potential_energy, digits=2)
            elapsed = round(time() - t0, digits=2)
            println("Step: $step/$N_steps, ETA: $(eta)s, Elapsed time: $(elapsed) seconds, Epot = $(Epot)")
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
        if step >= N_steps+start_step
            break
        end
    end

end