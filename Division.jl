### The functions in this file change the structure of the system states (number of particles)

include("System.jl")

function update_cell_list_division!(new_indices, sysGeometry::SystemGeometry, sysState::SystemState, sysDistances::SystemDistances)
    Lx, Ly = sysGeometry.limits
    nx = maximum(key[1] for key in keys(sysDistances.cell_to_index))
    ny = maximum(key[2] for key in keys(sysDistances.cell_to_index))
    rx, ry = Lx/nx, Ly/ny
 
    for k in new_indices
        #@printvar cell_to_index
        x,y = sysState.positions[1][k], sysState.positions[2][k]
        i = Int(x÷rx)+1
        j = Int(y÷ry)+1
        push!(sysDistances.cell_to_index[(i,j)], k)
        push!(sysDistances.index_to_cell, (i,j))
    end
 end

 


 function update_verlet_list_division!(indices_divided, new_indices, sysDistances::SystemDistances)
    n_new = length(indices_divided)
    for i in 1:n_new
        index_old = indices_divided[i]
        index_new = new_indices[i]
        neighbor_old = sysDistances.verlet_list[index_old]
        # neighbors of the new cell are also neighbors of the old cell
        push!(sysDistances.verlet_list, copy(neighbor_old))
        # the old cell is neighbor of the new cell
        push!(sysDistances.verlet_list[index_new], index_old)
    end
 
    # the cells that had old cells as neighbors also have the new cells as neighbors
    for i in 1:length(sysDistances.verlet_list)
        neighbors = sysDistances.verlet_list[i]
        for old_index in neighbors 
            l = findall(x->x==old_index, indices_divided)
            if length(l)>0
                new_index = new_indices[l[1]]
                if new_index < i
                    push!(sysDistances.verlet_list[i], new_index)
                elseif i<new_index
                    push!(sysDistances.verlet_list[new_index], i)
                end
 
            end
            l2 = findall(x->x==old_index, new_indices)
            if length(l2)>0
                old_index = indices_divided[l2[1]]
                if old_index < i
                    push!(sysDistances.verlet_list[i], old_index)
                elseif i<old_index
                    push!(sysDistances.verlet_list[old_index], i)
                end
            end
            
        end
    end
 end


 function divide!(sysGeometry::SystemGeometry, sysState::SystemState, sysDistances::SystemDistances, sysDynamic::SystemDynamic, sysHydro::SystemHydrodynamics,
    sysEvolution::SystemEvolution,
    paramInteractions::ParamInteractions, paramInit::ParamInitialisation, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation)
    
    #positions, previous_positions, velocities, previous_velocities, 
    #= velocities_dt2,
    forces,norm_force, sizes,previous_sizes, phases,displacement, 
    candidate_division, 
    initial_sizes, initial_phases,
    threshold, size_start,phase_force
    ; relative_distance=1/2, size_factor = 1/2, dimension=2) =#

    ### For intrinsic oscillator 
    threshold = 1
    index_new_cells = findall(>(threshold), sysState.phases)

  

    ### For external oscillator
    #= cells_above_threshold = phases .> 1-threshold
    index_new_cells = findall(cells_above_threshold.&candidate_division)
     =#
     ### If travelling wave, we should detect a maximum of the phase instead (phase(t-2)<phase(t-1) and phase(t-1)>phase(t))


    index_future_candidates = findall(<(threshold), sysState.phases)
    sysState.candidate_division[index_future_candidates] .= true


    for i in index_new_cells
        θ = 2π * rand()
        #δposition = [ sizes[i]* relative_distance *cos(θ), sizes[i]* relative_distance *sin(θ)]
        δposition = [ paramOscillations.size_start *cos(θ), paramOscillations.size_start *sin(θ)]

        if paramSimulation.division 

            for d in 1:sysGeometry.dimension
                push!(sysState.positions[d],sysState.positions[d][i] - δposition[d]  )
                push!(sysState.previous_positions[d], sysState.previous_positions[d][i])
                sysState.positions[d][i] += δposition[d] 
                push!(sysState.velocities[d], sysState.velocities[d][i])
                push!(sysState.previous_velocities[d], sysState.previous_velocities[d][i])
                push!(sysDynamic.forces[d], sysDynamic.forces[d][i])
                push!(sysDynamic.total_forces[d], sysDynamic.total_forces[d][i])
                
                push!(sysHydro.forces_hydro[d], sysHydro.forces_hydro[d][i])
                push!(sysState.velocities_dt2[d], sysState.velocities_dt2[d][i])
            end

            rescale_periodic!(sysState.positions)

            push!(sysDynamic.norm_forces, sysDynamic.norm_forces[i])

            sysState.phases[i] = sysState.phases[i]%1
            push!(sysState.phases,  sysState.phases[i])

            push!(sysDynamic.phase_force, sysDynamic.phase_force[i])
            #sizes[i] =sizes[i]*size_factor
            sysState.sizes[i] = paramOscillations.size_start

            push!(sysState.sizes, sysState.sizes[i])
            push!(sysState.previous_sizes, sysState.previous_sizes[i])

            sysState.candidate_division[i] = false
            push!(sysState.candidate_division, false)

            sysState.number_division[i]+=1
            push!(sysState.number_division, sysState.number_division[i])

            push!(sysState.parent_cell, i)

            # save which was the size before division
            if length(sysEvolution.sizes_before_division)<sysState.number_division[i]
                push!(sysEvolution.sizes_before_division, [])
            end
            push!(sysEvolution.sizes_before_division[sysState.number_division[i]],sysState.previous_sizes[i])


            if paramInteractions.add_hydrodynamics
                
                sysHydro.motility_matrix = hcat(sysHydro.motility_matrix, zeros(sysState.Nparticles))
                sysHydro.motility_matrix = vcat(sysHydro.motility_matrix, zeros(sysState.Nparticles+1)')

                sysHydro.resistance_matrix = hcat(sysHydro.resistance_matrix, zeros(sysState.Nparticles))
                sysHydro.resistance_matrix = vcat(sysHydro.resistance_matrix, zeros(sysState.Nparticles+1)')

                sysHydro.Id_plus_R_m1 = hcat(sysHydro.Id_plus_R_m1, zeros(sysState.Nparticles))
                sysHydro.Id_plus_R_m1 = vcat(sysHydro.Id_plus_R_m1, zeros(sysState.Nparticles+1)')

            end

            sysDistances.distances = hcat(sysDistances.distances, zeros(sysState.Nparticles))
            sysDistances.distances = vcat(sysDistances.distances, zeros(sysState.Nparticles+1)')

            #initial_sizes[i] = initial_sizes[i]/2
            #push!(initial_sizes, initial_sizes[i])
            #initial_phases[i] = phases[i]
            #push!(initial_phases, initial_phases[i])

            # we don't need to update the actual value of the displacemennt, because the
            # maximal_displacement is known: it's size_start
            push!(sysState.displacements, 0)
            sysState.Nparticles +=1

        else
            #@printvar i, sysState.positions[1][i], length(sysState.positions[1])
            sysState.phases[i] = sysState.phases[i]%1
            
            sysState.sizes[i] = paramOscillations.size_start
           
            sysState.candidate_division[i] = false
    
            sysState.number_division[i]+=1

            if length(sysEvolution.sizes_before_division)<sysState.number_division[i]
                push!(sysEvolution.sizes_before_division, [])
            end
            push!(sysEvolution.sizes_before_division[sysState.number_division[i]],sysState.previous_sizes[i])

        end

    end

    #sysState.Nparticles += length(index_new_cells)

    # We update maximal displacement outside of the function for efficiency
    #maximal_displacement += compute_maximal_displacement(positions, previous_positions)
    #if length(index_new_cells)>0
    #    maximal_displacement += size_start #compute_maximal_displacement(displacement)
    #end
    

    #if length(index_new_cells)>0
        #packing_fraction = compute_packing_fraction(sizes)
        #@printvar compute_packing_fraction
        #print("division of $(length(index_new_cells)) cells")
    #end

    #potential_energy, neighbor_evaluation, maximal_displacement = evaluate_force!(positions, sizes, verlet_list, 
    #forces, U, ∇U, cutoff, U0, diffusion, maximal_displacement, 
    #previous_positions, verlet_size, cell_to_index, index_to_cell, distances)

    return(index_new_cells)
end
