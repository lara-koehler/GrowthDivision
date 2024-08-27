const default_δParticlePhase = 0;
using Parameters

Base.@kwdef struct ParamInitialisation
    Nparticles::Int64
    ParticleSize::Float64
    δParticleSize::Float64
    δParticlePosition::Float64
    kinetic_energy::Float64
    random_initial_phase::Bool
    ParticlePhase::Float64 # if the particles all have the same phase 
    δParticlePhase::Union{Float64, Nothing} = default_δParticlePhase
end


Base.@kwdef struct ParamSimulation 
    dt::Float64
    N_steps::Int64
    ϕmax::Float64
    integrator::String
    max_duration_min::Float64
    saving_period::Int64
    division::Bool
    directoryResults::String
    nameSimulation::String
end

Base.@kwdef struct ParamInteractions
    U0::Float64
    potential_name::String # "soft" or "LJ" potentials are coded
    mass::Float64
    ξ::Float64 #this parameter is irrelevant if add_hydrodynamics=true
    cutoff_ratio::Float64 # cutoff for forces evaluation = maximal_size*cutoff_ratio
    diffusion::Float64
    radius_motility::Float64
    ϵ::Float64
    add_hydrodynamics::Bool
    viscosity::Float64

end

Base.@kwdef struct ParamOscillations
    amplitude::Float64 
    ω::Float64
    k::Float64
    x0_wave::Float64 
    y0_wave::Float64
    #wave_function = param_oscillations["wave_function"]
    growth_rate::Float64 
    size_start::Float64 
    growth_ratio::Float64
end

