include("System.jl")
include("Parameters.jl")
include("Evolution.jl")

function velocity_verlet!(
    sysGeometry::SystemGeometry,
    sysState::SystemState, 
    sysDynamic::SystemDynamic,
    sysDistances::SystemDistances,
    sysHydro::SystemHydrodynamics,
    paramSimulation::ParamSimulation,
    paramInteractions::ParamInteractions)
    ### this is using the SI of Tjhung2017SoftMatter paper

    #print_dimensions(sysState, sysDynamic, sysHydro, sysDistances)

    
    evaluate_total_force!(sysState, sysDynamic, sysHydro, paramInteractions ; evaluation_at_dt2=false)
    # this does not recompute the forces (i.e. the interparticle distances) but only sums the potential forces with the hydrodynamic one
    # (that result from multiplying the hydrodynamic resistance by the velocities, or just 1/ξ by the velocities)


    for d in 1:sysGeometry.dimension 
        sysState.positions[d], sysState.previous_positions[d] = sysState.positions[d] .+ (sysState.velocities[d] .* paramSimulation.dt) .+ sysDynamic.total_forces[d] .* paramSimulation.dt^2 /(2*paramInteractions.mass), sysState.positions[d]
    end

    ### the displacement is recomputed without the boundary conditions
    sysState.displacements.= sqrt.((sysState.positions[1].-sysState.previous_positions[1]).^2 + (sysState.positions[2].-sysState.previous_positions[2]).^2)

    ### then we apply the boundary conditions
    rescale_periodic!(sysState.positions, limits=sysGeometry.limits)

    ### evaluate the velocities at t+dt/2 through Taylor exapansion 
    for d in 1:sysGeometry.dimension
        sysState.velocities_dt2[d].= sysState.velocities[d] 
                    .+ sysDynamic.total_forces[d]/(2*paramInteractions.mass)*(paramSimulation.dt) 
        #push!(v_dt2, copy(velocities[d])*(1-ξ/(2m)*dt) + copy(forces[d])*dt/(2m) )
    end
    
    ### update the forces from the potential 
    neighbor_evaluation = evaluate_force!(sysGeometry, sysState, sysDistances, sysDynamic, sysHydro, paramInteractions,paramSimulation)

    ### update the total forces with the velocities at t+dt/2 
    #evaluate_total_force!(sysState, sysDynamic, sysHydro, paramInteractions; evaluation_at_dt2=true)
    

    ### update the velocities at time t+dt
    ### Here we use forces and not total_forces because the contribution to the hydrodynamic forces is included in the Id_plus_R_m1 matrix
    for d in 1:sysGeometry.dimension
        if paramInteractions.add_hydrodynamics
            sysState.velocities[d], sysState.previous_velocities[d] = sysHydro.Id_plus_R_m1 * (sysState.velocities_dt2[d] + sysDynamic.forces[d]*paramSimulation.dt/(2*paramInteractions.mass)), sysState.velocities[d]
        else
            sysState.velocities[d], sysState.previous_velocities[d] = (1 .+ (6π*param_interactions.viscosity .*sysState.sizes ) .*paramSimulation.dt/(2*paramInteractions.mass)).^(-1) .* (sysState.velocities_dt2[d] .+ sysDynamic.forces[d]*paramSimulation.dt/(2*paramInteractions.mass)), sysState.velocities[d]
        end
        #( sysState.velocities_dt2[d] .+ sysDynamic.total_forces[d]/(2*paramInteractions.mass)*(paramSimulation.dt), sysState.velocities[d]
    end
#=     @printvar paramInteractions.ξ
    @printvar sysHydro.resistance_matrix[1,1]
    @printvar sysHydro.Id_plus_R_m1[5,5]
    @printvar (1+paramInteractions.ξ*paramSimulation.dt/(2*paramInteractions.mass))^(-1)
 =#
    return(neighbor_evaluation)
end

