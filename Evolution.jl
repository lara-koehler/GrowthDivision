include("System.jl")
using LinearAlgebra
#using .System
#using .SystemState


function update_cutoff!(sysState::SystemState , sysDistances::SystemDistances, paramInteractions::ParamInteractions)
    sysDistances.cutoff = maximum(sysState.sizes)*paramInteractions.cutoff_ratio
    sysDistances.verlet_size=1.1*sysDistances.cutoff
end



function evaluate_force!(sysGeometry::SystemGeometry, 
    sysState::SystemState, 
    sysDistances::SystemDistances, 
    sysDynamic::SystemDynamic, 
    sysHydro::SystemHydrodynamics,
    paramInteractions::ParamInteractions,
    paramSimulation::ParamSimulation)
    
    #n, positions,sizes, displacement, verlet_list, forces, norm_force, 
    #U_name, cutoff, U0, diffusion, maximal_displacement, maximal_growth,
    # previous_positions, previous_sizes,
    #verlet_size, cell_to_index, index_to_cell, distances, cutoff_ratio,
    # phase_force, ϵ,dϕ, phases
    #; Lx=1.0, Ly=1.0, dimension=2)

    #n = length(positions[1])
    #dimension = length(positions)

    sysDynamic.maximal_displacement += maximum(sysState.displacements)
    sysDynamic.maximal_growth += maximum(sysState.sizes.-sysState.previous_sizes)

    if sysDynamic.maximal_displacement+sysDynamic.maximal_growth > (sysDistances.verlet_size-sysDistances.cutoff)#|| maximum(sysState.displacements)>minimum(sysState.sizes)
        
        #sysDistances .= initalize_system_distances(sysGeometry, param_interactions, sysState)

        update_cutoff!(sysState, sysDistances,paramInteractions)
        neighbor_evaluation = true
        sysDistances.distances .= initialize_pairwise_distance(sysState.positions)
        update_cell_list!(sysGeometry, sysState, sysDistances)
        update_verlet_list!(sysGeometry, sysState, sysDistances, param_interactions)

        #distances,cell_to_index, index_to_cell,verlet_list = initialize_verlet_and_cell_lists(positions, verlet_size, new_cutoff, n)
        sysDynamic.maximal_displacement = 0
        sysDynamic.maximal_growth = 0
    else 
        neighbor_evaluation = false
        #new_cutoff = cutoff
    end

    ### Reinitialize the force vectors 
    for i in 1:sysState.Nparticles
        sysDynamic.norm_forces[i]=0
        sysDynamic.phase_force[i]=0
        for d in 1:sysGeometry.dimension
            sysDynamic.forces[d][i]=0
        end
    end

    ### Reinitalize the motility matrix if add_hydrodynamics
    if paramInteractions.add_hydrodynamics
        sysHydro.motility_matrix .= zeros(sysState.Nparticles, sysState.Nparticles)
        #= for i in 1:sysState.Nparticles
            for d in 1:sysGeometry.dimension
                sysHydro.forces_hydro[d][i] = 0
            end
        end =#
    end

    sysDynamic.potential_energy = 0

    # Loop over the neighbor and update both the particles interaction and the phase coupling
    for i in 1:sysState.Nparticles
        neighbors = sysDistances.verlet_list[i]
        for j in neighbors
            ri = [sysState.positions[1][i],sysState.positions[2][i]]
            rj = [sysState.positions[1][j],sysState.positions[2][j]]
            σi = sysState.sizes[i]
            σj = sysState.sizes[j]
            dist, dx, dy = distance_periodic(ri, rj, Lx=sysGeometry.limits[1], Ly=sysGeometry.limits[2])
            if dist<sysDistances.cutoff

                #dx, dy = ri-rj
                absolute_force = ∇rU(dist, param_interactions.U0, σi, σj, param_interactions.potential_name)
                fr =  absolute_force / dist
                
                sysDynamic.potential_energy += 2*U(dist,param_interactions.U0, σi, σj , param_interactions.potential_name)
                
                sysDynamic.forces[1][i]+= dx*fr
                sysDynamic.forces[2][i]+= dy*fr
                sysDynamic.forces[1][j]+= -dx*fr
                sysDynamic.forces[2][j]+= -dy*fr

                sysDynamic.norm_forces[i] += abs(absolute_force)
                sysDynamic.norm_forces[j] += abs(absolute_force)

                #fϕ = (dist-σi-σj< min(σi,σj)/10) * sin((sysState.phases[i]-sysState.phases[j])*π)
                fϕ = (1/dist^2) * sin((sysState.phases[i]-sysState.phases[j])*π)
         
                sysDynamic.phase_force[i] += -fϕ
                sysDynamic.phase_force[j] += +fϕ

                if paramInteractions.add_hydrodynamics && dist<sysHydro.cutoff_hydro
                    sysHydro.motility_matrix[i,j] = -1/dist * 1/(8π*paramInteractions.viscosity)
                    sysHydro.motility_matrix[j,i] = -1/dist * 1/(8π*paramInteractions.viscosity)
                end
            end
        end
    end

    sysDynamic.phase_force .= sysDynamic.phase_force .* param_interactions.ϵ

   

    # Add Brownian forces
    for d in 1:sysGeometry.dimension 
        sysDynamic.forces[d].+= rand(Normal(0, 1), sysState.Nparticles) * sqrt(2*param_interactions.diffusion)
    end


    # Add hydrodynamic forces 
    if paramInteractions.add_hydrodynamics
        evaluate_hydrodynamic_forces!(sysState, sysHydro, paramInteractions, paramSimulation)
        #= for d in 1:sysGeometry.dimension 
            sysDynamic.forces[d] .+= sysHydro.forces_hydro[d]
        end =#
    end

    return(neighbor_evaluation)

end


function evaluate_hydrodynamic_forces!(sysState::SystemState, 
    #sysDistances::SystemDistances, 
    #sysDynamic::SystemDynamic, 
    sysHydro::SystemHydrodynamics,
    paramInteractions::ParamInteractions,
    paramSimulation::ParamSimulation)
    ### This assumes that the non diagonal terms of the motility matrix were already computed in the evaluate_force function

    
    ### Diagonal terms of the motility matrix
    for i in 1:sysState.Nparticles
        sysHydro.motility_matrix[i,i] = -1/(6 * π * sysState.sizes[i] * param_interactions.viscosity)
    end

    

    ### Inverting the matrix 
    sysHydro.resistance_matrix .= inv(sysHydro.motility_matrix)


    ### calculating (Id + dt/2m * R)^-1 which is useful for the integration 
    ### This might be optimized in the future to calculate this matrix from the motility matrix instead of doing two inversions
    sysHydro.Id_plus_R_m1 .= inv(Diagonal(ones(sysState.Nparticles))-paramSimulation.dt/(2*paramInteractions.mass)*sysHydro.resistance_matrix)
   
    
    #= if maximum(abs.(sysHydro.motility_matrix))>1000
        @printvar sysHydro.motility_matrix
        @printvar sysHydro.resistance_matrix
        @printvar sysHydro.Id_plus_R_m1
    end =#

    ### Calculating the force vectors
   #=  sysHydro.forces_hydro[1] = sysHydro.resistance_matrix*sysState.velocities[1]
    sysHydro.forces_hydro[2] = sysHydro.resistance_matrix*sysState.velocities[2] =#
    
end


function evaluate_total_force!(sysState::SystemState, 
    sysDynamic::SystemDynamic, 
    sysHydro::SystemHydrodynamics, 
    paramInteractions::ParamInteractions; evaluation_at_dt2=false)
    # this sums the potential force and the friction force, using the velocities 
    # at instant t or t+dt/2 depending on the stage of the intergration (cf SI of Tjhung2017SoftMatter)

    if paramInteractions.add_hydrodynamics
        # we use the hydrodynamic resistance to compute the velocity-dependent term  (this is a matrix multiplication)
        for d in 1:sysGeometry.dimension 
            if !evaluation_at_dt2
                sysHydro.forces_hydro[d] .= sysHydro.resistance_matrix*sysState.velocities[d]
                sysDynamic.total_forces[d] .= sysDynamic.forces[d] .+ sysHydro.forces_hydro[d]
                #@printvar sysHydro.resistance_matrix[1,1]
                #@printvar 6 * π * sysState.sizes[1] * param_interactions.viscosity
            else
                sysDynamic.total_forces[d] .= sysDynamic.forces[d] .+ sysHydro.resistance_matrix*sysState.velocities_dt2[d]
            end
        end
    
    else
        # we use the hydrodynamic resistance to compute the velocity-dependent term (this is a scalar multiplication)
        for d in 1:sysGeometry.dimension 
            if !evaluation_at_dt2
                sysHydro.forces_hydro[d] .= -(6π*param_interactions.viscosity .*sysState.sizes ) .*sysState.velocities[d]
                sysDynamic.total_forces[d] .= sysDynamic.forces[d] .+ sysHydro.forces_hydro[d]
            else
                sysDynamic.total_forces[d] .= sysDynamic.forces[d] .- (6π*param_interactions.viscosity .*sysState.sizes ) .* sysState.velocities_dt2[d]
            end
        end
    end

end





function update_cell_list!(sysGeometry::SystemGeometry, sysState::SystemState, sysDistances::SystemDistances)
   
   # cell_to_index is a dictionnary, the keys are the indices of the cells
   # index_to_cell is a list of tuple 
   Lx, Ly = sysGeometry.limits

   rx, ry = sysDistances.verlet_size, sysDistances.verlet_size
   nx, ny = Int(Lx÷rx)+Int(Lx%rx>1e-8), Int(Ly÷ry)+Int(Ly%ry>1e-8)
   sysDistances.nx_cell,sysDistances.ny_cell = nx, ny
   
   #= nx = maximum(key[1] for key in keys(sysDistances.cell_to_index))
   ny = maximum(key[2] for key in keys(sysDistances.cell_to_index))
   rx, ry = Lx/nx, Ly/ny =#
   n = sysState.Nparticles

   # clear the dictionnary of indices of cell to indices of particles 
   empty(sysDistances.cell_to_index)

   for i in 1:nx
        for j in 1:ny
            sysDistances.cell_to_index[(i,j)] = Set()
        end
    end



   for k in 1:n
       #i_previous, j_previous = sysDistances.index_to_cell[k]
       x,y = sysState.positions[1][k], sysState.positions[2][k]
       i = Int(x÷rx)+1
       j = Int(y÷ry)+1
       
       #if i_previous!=i || j_previous!=j
        #if haskey(sysDistances.cell_to_index, (i,j))
            push!(sysDistances.cell_to_index[(i,j)], k)
       #else
        #sysDistances.cell_to_index[(i,j)] = Set([k])
       #end
       #delete!(sysDistances.cell_to_index[(i_previous,j_previous)], k)
       sysDistances.index_to_cell[k] = (i,j)
       #end
   end
end






function update_verlet_list!(sysGeometry::SystemGeometry, 
    sysState::SystemState, 
    sysDistances::SystemDistances, 
    param_interactions::ParamInteractions)

   n = length(sysState.positions[1])
   #println("verlet")
   #nx = maximum(key[1] for key in keys(sysDistances.cell_to_index))
   #ny = maximum(key[2] for key in keys(sysDistances.cell_to_index))
   nx, ny = sysDistances.nx_cell, sysDistances.ny_cell
   #println(raw"nx, ny "*string(nx, ny))
   for i_cell in 1:nx
       for j_cell in 1:ny
           for add in [(0,0), (0,1), (1,0), (1,1), (-1,1)]
               i_next,j_next = i_cell+add[1], j_cell+add[2]
               i_next = i_next < 1 ? i_next+nx : i_next

           #for i_next in i_cell:i_cell+1
               #for j_next in j_cell:j_cell+1
               cell1 = (i_cell, j_cell)
               cell2 = ((i_next-1)%nx +1, (j_next-1)%ny +1)
               #println("cells"*string(cell1)*" "*string(cell2))
               for k1 in sysDistances.cell_to_index[cell1]
                   for k2 in sysDistances.cell_to_index[cell2]
                       p1, p2 = max(k1, k2), min(k1, k2)
                       if p1!=p2
                           #println(string(p1)*" "*string(p2)*" "*string(distances[p1,p2] < cutoff ))
                           if sysDistances.distances[p1,p2] < sysDistances.cutoff
                               push!(sysDistances.verlet_list[p1], p2)
                           else
                               delete!(sysDistances.verlet_list[p1],p2)
                           end
                       end
                   end
               end
           end
       end
   end
end




""" PHASE EVOLUTION """


function update_phases_intrinsic!(sysState::SystemState, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation)

    """ Phases are governed by an intrinsic oscillator"""
    sysState.phases.+= paramOscillations.ω * paramSimulation.dt 
    if !paramSimulation.division
        sysState.phases .= sysState.phases.%1
    end
end
function update_phases_synchronization(sysState::SystemState, sysDynamic::SystemDynamic, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation,paramInteractions::ParamInteractions)

    sysState.phases.+= paramOscillations.ω * paramSimulation.dt  .+ (max.(sysDynamic.phase_force,0) * paramSimulation.dt) 
    #.+ (paramOscillations.ω .* rand(sysState.Nparticles) .* paramSimulation.dt)
    #= for i in 1:sysState.Nparticles
        #sysState.phases[i] += paramOscillations.ω * paramSimulation.dt
        dphi =  paramOscillations.ω * paramSimulation.dt * rand()
        #@printvar  dphi
        sysState.phases[i] += dphi
    end =#

    #sysState.phases.+=   ( paramOscillations.ω * paramSimulation.dt .*(1 .+ (2 .- rand(sysState.Nparticles)) ))
    #@printvar sysState.phases
    #phases .= max.(phases, 0)

end
function update_phases_external!(phases, t, positions, wave_function, wave_parameters)
    phases.= wave!(wave_function, phases, t, positions[1], positions[2], wave_parameters)
end



""" SIZE EVOLUTION """



function update_sizes_periodic!(phases, sizes, initial_sizes, initial_phases, amplitude)
    sizes .= initial_sizes  .+ amplitude/2 .* (1 .- cos.((phases-initial_phases) * 2π ))
end

function update_sizes_growth!(sysState::SystemState, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation)
    """the sizes grow at a constant rate growth_rate"""
    sysState.previous_sizes.= sysState.sizes
    sysState.sizes .+=  paramOscillations.growth_rate*paramSimulation.dt
end

function update_sizes_growth_feedback!(sysState::SystemState, sysDynamic::SystemDynamic, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation, paramInteractions::ParamInteractions)
    """the sizes grow at a rate that depends of the forces that is feeled, 
    when the disks  feels 0 force, it grows at speed growth_rate
    the size cannot decrease below the critical size"""
    sysState.previous_sizes.= sysState.sizes
    #sizes .+=  growth_rate*dt .* max.(1 .- radius_motility.*norm_force, 0)
    sysState.sizes .+=  paramOscillations.growth_rate*paramSimulation.dt .* (1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces)
    sysState.sizes .= max.(sysState.sizes, paramOscillations.size_start)

end

function update_sizes_growth_feedback_transient!(sysState::SystemState, sysDynamic::SystemDynamic, paramOscillations::ParamOscillations, paramSimulation::ParamSimulation, paramInteractions::ParamInteractions)
    """everything like update_sizes_growth_feedback, but the growth only happens
    when the phase is above a 1-growth_ratio"""
    sysState.previous_sizes.= sysState.sizes
    #if minimum(1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces)<0.8
    #    @printvar minimum(1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces) 
    #    @printvar maximum(1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces) 
    #end
    
    sysState.sizes .+=  paramOscillations.growth_rate*paramSimulation.dt .* (1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces) .* .!((sysState.phases .< (1-paramOscillations.growth_ratio)) .& ( (1 .- paramInteractions.radius_motility.*sysDynamic.norm_forces ).> 0))
    sysState.sizes .= max.(sysState.sizes, paramOscillations.size_start)

end



function update_system_evolution!(sysEvolution::SystemEvolution, sysState::SystemState, sysDynamic::SystemDynamic, sysHydro::SystemHydrodynamics)
    push!(sysEvolution.successive_positions,deepcopy(sysState.positions))
    push!(sysEvolution.successive_velocities, deepcopy(sysState.velocities))
    push!(sysEvolution.successive_phases, deepcopy(sysState.phases))
    push!(sysEvolution.successive_sizes, deepcopy(sysState.sizes))
    push!(sysEvolution.successive_forces, deepcopy(sysDynamic.forces))
    push!(sysEvolution.successive_hydro_forces, deepcopy(sysHydro.forces_hydro))
    push!(sysEvolution.successive_phase_force, deepcopy(sysDynamic.phase_force))
    push!(sysEvolution.successive_norm_force, deepcopy(sysDynamic.norm_forces))
    push!(sysEvolution.successive_potential_energies, deepcopy(sysDynamic.potential_energy))
    push!(sysEvolution.successive_kinetic_energies, deepcopy(sysDynamic.kinetic_energy))
    push!(sysEvolution.successive_number_division, deepcopy(sysState.number_division))
            
end
    