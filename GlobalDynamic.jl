include("Parameters.jl")
include("System.jl")
include("Evolution.jl")
include("Integration.jl")
include("Division.jl")
using JSON
#using JSON3

function perfom_dynamics(
    sysGeometry::SystemGeometry,
    paramInteractions::ParamInteractions,
    paramInitialisation::ParamInitialisation,
    paramSimulation::ParamSimulation,
    paramOscillations::ParamOscillations ;
    sysState_init::SystemState,
    random_initial_state=true,
    profiling = false)

    if random_initial_state 
        sysState = initialize_system_state(sysGeometry, paramInitialisation, paramInteractions)
    else
        sysState = sysState_init
    end

    sysDynamic = initialize_system_dynamic(sysGeometry, paramInitialisation)
    sysDistances = initalize_system_distances(sysGeometry, paramInteractions, sysState)
    sysHydro = initialize_system_hydrodynamics(sysGeometry, sysState)

    N_neighbor_evaluation = 0


    # evaluate the inital force and distances
    update_cell_list!(sysGeometry, sysState, sysDistances)
    update_verlet_list!(sysGeometry, sysState, sysDistances, param_interactions);
    evaluate_force!(sysGeometry, sysState, sysDistances, sysDynamic, sysHydro, paramInteractions,paramSimulation)
    

    sysEvolution = initialize_system_evolution(sysState, sysDynamic, sysHydro)

    steps_neighbor_evaluation = []
    steps_performed = 1
    percentage = 0.01
    t0 = time()
    t=0
    ϕ = compute_packing_fraction(sysState.sizes, Lx=sysGeometry.limits[1], Ly=sysGeometry.limits[2])

    n_div_since_eval = 0

    for step in 1:paramSimulation.N_steps

        if step/paramSimulation.N_steps > percentage
            percentage += 0.1
            t_current = time()
            println("$(round(percentage*100)) % in $(round((t_current-t0)/60, digits=2)) min")
            ϕ = compute_packing_fraction(sysState.sizes, Lx=sysGeometry.limits[1], Ly=sysGeometry.limits[2])
            println("n_particles=$(sysState.Nparticles), Φ=$ϕ")
        end

        """Verify that the packing fraction is not too large"""
        if ϕ>paramSimulation.ϕmax || (time()-t0)>paramSimulation.max_duration_min*60 || maximum(sysState.displacements)>100#20*maximum(sysState.sizes)
            if ϕ>paramSimulation.ϕmax
                println("Stop because the packing fraction reached its maximum")
            elseif (time()-t0)>paramSimulation.max_duration_min*60
                println("Stop because simulation is taking too much time")
            else 
                println("Stop because the displacement was too large")
                @printvar maximum(sysState.displacements)
                #@printvar 20*maximum(sysState.sizes)
            end
            break
        end



       
        """Integrate"""
        if paramSimulation.integrator=="velocity_verlet"
            neighbor_evaluation = velocity_verlet!( sysGeometry,
            sysState, 
            sysDynamic,
            sysDistances,
            sysHydro,
            paramSimulation,
            paramInteractions)
            N_neighbor_evaluation += neighbor_evaluation
            
            
        end
        
        if paramSimulation.integrator=="euler_semi_implicit"
            implicit_euler_integration!(positions, previous_positions, velocities, previous_velocities, forces, m, dt, 0)
            potential_energy, neighbor_evaluation, maximal_displacement,maximal_growth, 
            distances,cell_to_index, index_to_cell,verlet_list,verlet_size, cutoff = evaluate_force!(n_particles, positions, sizes, displacement,
            verlet_list, forces,norm_force, U_name, cutoff, U0, diffusion,
            maximal_displacement,maximal_growth, previous_positions, previous_sizes, 
            verlet_size, cell_to_index, index_to_cell, distances,cutoff_ratio,
            phase_force, ϵ,dϕ,phases)
            
            N_neighbor_evaluation += neighbor_evaluation
        end

        if neighbor_evaluation
            push!(steps_neighbor_evaluation, step)
        end
       
        
        """Update of the phases """
        #update_phases_intrinsic!(phases, ω, dt)
        #update_phases_external!(phases, t, positions, wave_function, Dict("k"=>k, "ω"=>ω, "x0"=>x0_wave, "y0"=>y0_wave))
        update_phases_synchronization(sysState, sysDynamic, paramOscillations, paramSimulation, paramInteractions)

        """Update of the particle sizes"""
        #update_sizes_growth!(sizes, previous_sizes,growth_rate, dt)
        #update_sizes_growth_feedback!(sysState, sysDynamic, paramOscillations, paramSimulation, paramInteractions)
        update_sizes_growth_feedback_transient!(sysState, sysDynamic, paramOscillations, paramSimulation, paramInteractions)
        

        """Division"""
        
        #print_dimensions(sysState, sysDynamic, sysHydro)
        index_new_cells = divide!(sysGeometry, sysState, sysDistances, sysDynamic, sysHydro,sysEvolution,
        paramInteractions, paramInitialisation, paramOscillations, paramSimulation)
        #print_dimensions(sysState, sysDynamic, sysHydro)

        if paramSimulation.division
            n_division = length(index_new_cells)
            new_indices = sysState.Nparticles+1-n_division:sysState.Nparticles
            #n_particles+=n_division
            if n_division>0
                update_cell_list_division!(new_indices, sysGeometry, sysState, sysDistances)
                update_verlet_list_division!(index_new_cells, new_indices,  sysDistances)
                sysDynamic.maximal_displacement += paramOscillations.size_start
                n_div_since_eval += n_division
               
            end
        end

        if n_div_since_eval>1
            sysDynamic.maximal_displacement = 1e8
            # this will force the reevaluation of the neighbors
            n_div_since_eval=0
        end


        """Save the system state"""
        if (paramSimulation.saving_period!=0 && step%paramSimulation.saving_period==0)
            ϕ = compute_packing_fraction(sysState.sizes, Lx=sysGeometry.limits[1], Ly=sysGeometry.limits[2])
            sysDynamic.kinetic_energy = evaluate_kinetic_energy(sysState.velocities, paramInteractions.mass)
            update_system_evolution!(sysEvolution, sysState, sysDynamic, sysHydro)
        end

        steps_performed +=1
        t+=paramSimulation.dt
    end

    @printvar N_neighbor_evaluation
    @printvar steps_performed
    println("$(round((time()-t0)/60, digits=2)) min")
    append!(sysEvolution.parent_cell, sysState.parent_cell)

    sysEvolution.N_steps_performed = steps_performed
    sysEvolution.N_neighbors_evaluation = N_neighbor_evaluation

    directory = paramSimulation.directoryResults*"/"*paramSimulation.nameSimulation
    println("We will save results in the following directory")
    @printvar directory
    mkdir(directory)

  

    save_to_json(directory*"/SystemEvolution.json",sysEvolution)
    save_to_json(directory*"/SystemGeometry.json",sysGeometry)
    save_to_json(directory*"/ParamInitialisation.json", paramInitialisation)
    save_to_json(directory*"/ParamInteractions.json", paramInteractions)
    save_to_json(directory*"/ParamSimulation.json", paramSimulation)
    save_to_json(directory*"/ParamOscillations.json", paramOscillations)


    if !profiling
        return( sysEvolution )
    end

end


function save_to_json(file_path::String, data)
    json_str = JSON.json(data)
    open(file_path, "w") do file
        write(file, json_str)
    end
end

function load_from_json(file_path::String, T)
    json_str = read(file_path, String)
    json_dict = JSON.parse(json_str)

    attributes = fieldnames(T)
    values = [json_dict[string(attributes[k])] for k in 1:length(attributes)]

    values = [json_dict[string(attributes[k])] for k in 1:length(attributes)]

    return ( T(values...) )

end

function load_from_json_optional(file_path::String, T)
    json_str = read(file_path, String)
    json_dict = JSON.parse(json_str)
    attributes = fieldnames(T)
    values = []

    for k in 1:length(attributes)
        field_name = string(attributes[k])
        if haskey(json_dict, field_name)
            push!(values, json_dict[field_name])
        else
            # If the field is missing, assign the default value
            if field_name == "δParticlePhase"
                push!(values, default_δParticlePhase)
            else
                error("Missing required field: $field_name")
            end
        end
    end
    return ( T(values...) )

end