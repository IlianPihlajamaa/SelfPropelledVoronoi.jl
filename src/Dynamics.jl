

function run_simulation!(parameters, arrays, output)
    step = output.steps_done

    # Main simulation loop
    while true

        # invoke callback
        parameters.callback(step, parameters, arrays, output)

        do_time_step!(parameters, arrays, output)

        # Check if the simulation should be saved
        if parameters.dump_info.save
            if step in parameters.dump_info.when_to_save_array
                save_simulation_state!(parameters, arrays, output)
            end
        end

        # Check if the simulation should be tesselated
        if !verify_tesselation(parameters, arrays, output)
            voronoi_tesselation!(parameters, arrays, output)
        end

        # Check if the simulation should be stopped
        step = step + 1
        output.steps_done = step
        if step > parameters.Nsteps
            break
        end
    end

end