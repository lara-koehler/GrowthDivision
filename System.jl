include("Initialisation.jl")
include("Parameters.jl")
include("SideFunctions.jl")

macro printvar(var)
    ex = esc(var)
    varname = string(var)
    quote
        println(string($varname), " = ", $ex)
        $ex
    end
end

#module System 

#using .Parameters

#export SystemGeometry, SystemState, SystemDynamic, SystemDistances, initialize_system_state

struct SystemGeometry
    dimension::Int64
    limits::Vector{Float64}
end

Base.@kwdef mutable struct SystemState
    Nparticles::Int64
    positions::Vector{Vector{Float64}}  # Vector of 2D positions
    previous_positions::Vector{Vector{Float64}}
    displacements::Vector{Float64}
    velocities::Vector{Vector{Float64}}#,Vector{Float64}} # Vector of 2D velocities
    velocities_dt2::Vector{Vector{Float64}}
    previous_velocities::Vector{Vector{Float64}}#,Vector{Float64}} # Vector of 2D velocities
    sizes::Vector{Float64}  # Vector of particle sizes
    previous_sizes::Vector{Float64}  # Vector of particle sizes
    phases::Vector{Float64} # Vector of particle phases
    candidate_division::Vector{Int64} #indices of the particle that could divide
    number_division::Vector{Int64}
    parent_cell::Vector{Int64}
end

Base.@kwdef mutable struct SystemDynamic
    forces::Vector{Vector{Float64}}  # Vector of 2D positions
    norm_forces::Vector{Float64} # the absolute values of the forces
    total_forces::Vector{Vector{Float64}} # the other is only the potential forces
    maximal_displacement::Float64
    maximal_growth::Float64
    phase_force::Vector{Float64}    
    potential_energy::Float64
    kinetic_energy::Float64
end

Base.@kwdef mutable struct SystemDistances
    distances::Matrix{Float64}
    nx_cell::Int64 # the x-dimension of the cell table 
    ny_cell::Int64 # the y-dimension of the cell table 
    cell_to_index::Dict{Tuple{Int64,Int64}, Set{Int64}}
    index_to_cell::Vector{Tuple{Int64, Int64}}
    verlet_list::Vector{Set{Int64}}
    verlet_size::Float64
    cutoff::Float64

end

Base.@kwdef mutable struct SystemHydrodynamics
    motility_matrix::Matrix{Float64}
    resistance_matrix::Matrix{Float64}
    forces_hydro::Vector{Vector{Float64}}
    Id_plus_R_m1::Matrix{Float64}
    cutoff_hydro::Float64
end


Base.@kwdef mutable struct SystemEvolution
    successive_positions::Vector{Vector{Vector{Float64}}}
    successive_velocities::Vector{Vector{Vector{Float64}}}
    successive_phases::Vector{Vector{Float64}}
    successive_sizes::Vector{Vector{Float64}}
    successive_forces::Vector{Vector{Vector{Float64}}}
    successive_hydro_forces::Vector{Vector{Vector{Float64}}}
    successive_norm_force::Vector{Vector{Float64}}
    successive_phase_force::Vector{Vector{Float64}}
    successive_kinetic_energies::Vector{Float64}
    successive_potential_energies::Vector{Float64}
    successive_number_division::Vector{Vector{Int64}}
    parent_cell::Vector{Int64} # this is useless, but we keep it to be able to read the old data
    sizes_before_division::Vector{Vector{Float64}}
    N_steps_performed::Int64
    N_neighbors_evaluation::Int64

end



function initialize_system_state(sysGeometry::SystemGeometry, 
    param_initalisation::ParamInitialisation, param_interactions::ParamInteractions)
    Lx, Ly = sysGeometry.limits

    # initalize the positions
    positions = position_lattice_noisy(param_initalisation.Nparticles, Lx=Lx, Ly=Ly,  η=param_initalisation.δParticlePosition)

    # initalize the velocities
    velocities = initialize_velocities(param_initalisation.Nparticles, param_initalisation.kinetic_energy, param_interactions.mass)

    # initalize the phases 
    if param_initalisation.random_initial_phase
        #phases = initialize_random_phases(param_initalisation.Nparticles)
        phases = initialize_random_phases_gaussian_positive(param_initalisation.Nparticles, param_initalisation.ParticlePhase, param_initalisation.δParticlePhase)
    else
        phases = initialize_uniform_phases(param_initalisation.Nparticles).+param_initalisation.ParticlePhase
    end

    #initalize the sizes
    sizes = initialize_polydisperse_sizes(param_initalisation.Nparticles, param_initalisation.ParticleSize, param_initalisation.δParticleSize)

    candidate_division = initialize_candidate_division(param_initalisation.Nparticles)

    number_divison = zeros(param_initalisation.Nparticles)
    parent_cell = [k for k in 1:param_initalisation.Nparticles]

    sysState = SystemState(Nparticles=param_initalisation.Nparticles, 
    positions=positions, 
    previous_positions=deepcopy(positions),
    displacements=zeros(param_initalisation.Nparticles),
    velocities=velocities,
    velocities_dt2=deepcopy(velocities),
    previous_velocities=deepcopy(velocities),
    sizes=sizes,
    previous_sizes=deepcopy(sizes),
    phases=phases,
    candidate_division=candidate_division,
    number_division = number_divison,
    parent_cell=parent_cell)

    return(sysState)
end



function initialize_empty_system_state()
    sysState = SystemState(Nparticles=0, 
    positions=[[]], 
    previous_positions=[[]],
    displacements=[],
    velocities=[[]],
    velocities_dt2=[[]],
    previous_velocities=[[]],
    sizes=[],
    previous_sizes=[],
    phases=[],
    candidate_division=[],
    number_division = [],
    parent_cell=[])
    return(sysState)
end

function rescale_phases!(phases, Δφ_ratio ; average_phase=nothing )
    if average_phase == nothing 
        average_phase = mean(phases)
    end

    phases .= (phases .- average_phase) .* Δφ_ratio .+ average_phase
    for i in 1:length(phases)
        if phases[i]<0
            phases[i]=-phases[i]
        elseif phases[i]>1
            phases[i]=2-phases[i]
        end
    end
    return(phases)
end



function initialize_state_from_evolution(sysEvolution::SystemEvolution, index; Δφ_ratio = 1, average_phase=nothing)

    phases = deepcopy(sysEvolution.successive_phases[index])
    if abs(Δφ_ratio-1)>1e-3
        phases .= rescale_phases!(phases, Δφ_ratio, average_phase=average_phase)
    end

    return( SystemState(Nparticles=length(sysEvolution.successive_positions[index][1]),
    positions=deepcopy(sysEvolution.successive_positions[index]),
    previous_positions=deepcopy(sysEvolution.successive_positions[index]),
    displacements=zeros(length(sysEvolution.successive_positions[index][1])),
    velocities=deepcopy(sysEvolution.successive_velocities[index]),
    velocities_dt2 = deepcopy(sysEvolution.successive_velocities[index]),
    previous_velocities=deepcopy(sysEvolution.successive_velocities[index]),
    sizes=deepcopy(sysEvolution.successive_sizes[index]),
    previous_sizes = deepcopy(sysEvolution.successive_sizes[index]),
    phases = phases,
    candidate_division=initialize_candidate_division(length(sysEvolution.successive_positions[index][1])), 
    number_division=deepcopy(sysEvolution.successive_number_division[index]), 
    parent_cell=deepcopy(sysEvolution.parent_cell)))

end


function initialize_system_dynamic(sysGeometry::SystemGeometry, param_initalisation::ParamInitialisation)
    forces = initialize_forces(param_initalisation.Nparticles ; sysGeometry.dimension)
    norm_forces = zeros(param_initalisation.Nparticles)
    total_forces = initialize_forces(param_initalisation.Nparticles ; sysGeometry.dimension)
    maximal_displacement = 0
    maximal_growth = 0
    phase_force = zeros(param_initalisation.Nparticles)

    return( SystemDynamic(forces=forces, norm_forces=norm_forces, total_forces=total_forces,
    maximal_displacement=maximal_displacement,
    maximal_growth=maximal_growth,
    phase_force=phase_force,
    potential_energy=0,
    kinetic_energy=0)
 )

end

function initialize_system_dynamic_empty()
    return( SystemDynamic(forces=[[]], norm_forces=[], total_forces=[[]],
    maximal_displacement=0,
    maximal_growth=0,
    phase_force=[],
    potential_energy=0,
    kinetic_energy=0))
end

function initalize_system_distances(sysGeometry::SystemGeometry, 
    param_interactions::ParamInteractions,
    sysState::SystemState)
   
    distances = initialize_pairwise_distance(sysState.positions)
    cutoff = maximum(sysState.sizes)*param_interactions.cutoff_ratio
    verlet_size = 1.1*cutoff
    nx, ny, cell_to_index, index_to_cell =  initialize_cell_list(sysState.positions, verlet_size, verlet_size)
    verlet_list = initialize_verlet_list(sysState.Nparticles) 
   
    #cutoff = maximum(sysState.sizes)*param_interactions.cutoff_ratio
    #update_cell_list!(sysState.positions, cell_to_index, index_to_cell)
    #update_verlet_list!(sysState.positions, distances, verlet_list,cell_to_index, param_interactions.cutoff)

    
    return( SystemDistances(distances=distances, 
    nx_cell=nx,
    ny_cell=ny,
    cell_to_index=cell_to_index, 
    index_to_cell=index_to_cell, 
    verlet_list=verlet_list, 
    verlet_size=verlet_size, 
    cutoff=cutoff ))
    
end


function initialize_system_hydrodynamics(sysGeometry::SystemGeometry, sysState::SystemState)
    motility_matrix=zeros(sysState.Nparticles, sysState.Nparticles)
    resistance_matrix=zeros(sysState.Nparticles, sysState.Nparticles)
    forces_hydro = initialize_forces(sysState.Nparticles ; sysGeometry.dimension)
    Id_plus_R_m1 = zeros(sysState.Nparticles, sysState.Nparticles)
    return( SystemHydrodynamics(motility_matrix=motility_matrix, 
    resistance_matrix=resistance_matrix, forces_hydro=forces_hydro, Id_plus_R_m1=Id_plus_R_m1, cutoff_hydro=1.0))
end


function initialize_system_evolution(sysState::SystemState, sysDynamic::SystemDynamic, sysHydro::SystemHydrodynamics)
    return( SystemEvolution(
        successive_positions=[deepcopy(sysState.positions)],
        successive_velocities=[deepcopy(sysState.velocities)],
        successive_phases = [deepcopy(sysState.phases)],
        successive_sizes = [deepcopy(sysState.sizes)],
        successive_forces = [deepcopy(sysDynamic.forces)], 
        successive_hydro_forces = [deepcopy(sysHydro.forces_hydro)] ,
        successive_norm_force = [deepcopy(sysDynamic.norm_forces)] ,
        successive_phase_force = [deepcopy(sysDynamic.phase_force)], 
        successive_kinetic_energies = [deepcopy(sysDynamic.kinetic_energy)],
        successive_potential_energies = [deepcopy(sysDynamic.potential_energy)], 
        successive_number_division = [deepcopy(sysState.number_division)],
        parent_cell = [],
        sizes_before_division = [[]],
        N_steps_performed=0,
        N_neighbors_evaluation=0
        ))
end

function initialize_system_evolution_empty()
    return( SystemEvolution(

    successive_positions = [[[]]],
    successive_velocities = [[[]]],
    successive_phases = [[]],
    successive_sizes = [[]],
    successive_forces= [[[]]],
    successive_hydro_forces= [[[]]],
    successive_norm_force= [[]],
    successive_phase_force= [[]],
    successive_kinetic_energies = [],
    successive_potential_energies= [],
    successive_number_division= [[]],
    parent_cell= [], # this is useless, but we keep it to be able to read the old data
    sizes_before_division= [[]],
    N_steps_performed = 0,
    N_neighbors_evaluation = 0))

end



function print_dimensions(sysState::SystemState, sysDynamic::SystemDynamic, sysHydro::SystemHydrodynamics, sysDistances::SystemDistances)

    @printvar length(sysState.positions[1])
    @printvar length(sysState.previous_positions[1])
    @printvar length(sysState.velocities[1])
    @printvar length(sysDynamic.forces[1])
    @printvar length(sysHydro.forces_hydro[1])
    @printvar length(sysState.displacements)
    @printvar size(sysHydro.motility_matrix)
    @printvar size(sysDistances.distances)
    
end