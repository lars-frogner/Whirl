abstract type EnergyComponent{D} end

"""
    get_specific_energies(energy_component::EnergyComponent) -> Vector{Float}

Return the specific energies of the particles in the energy component.
"""
function get_specific_energies end

"""
    initderivatives(energy_component::EnergyComponent{ED}) -> ED

Creates initial derivatives for the energy component.
"""
initderivatives(::EnergyComponent{Nothing}) = nothing

"""
    initderivatives!(
        derivatives::ED,
        ::EnergyComponent{ED},
    )

Resets the given derivatives for the energy component.
"""
function initderivatives!(::Nothing, ::EnergyComponent{Nothing}) end

"""
    updatederivatives!(
        energy_derivatives::ED,
        (i, j)::Tuple{Unsigned,Unsigned},
        vᵣ::Number,
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
        energy_component::EnergyComponent{ED},
    )

Update energy component derivatives for particles `i` and `j`.
"""
function updatederivatives!(
    ::Nothing,
    ::Tuple{Unsigned,Unsigned},
    ::Number,
    ::Tuple{Number,Number,Number},
    ::EnergyComponent{Nothing},
) end

"""
    evolvevariables!(
        energy_component::EnergyComponent{ED},
        energy_derivatives::ED,
        Δt::Number,
    )

Evolve energy component with the given derivatives and time step.
"""
function evolvevariables!(::EnergyComponent{Nothing}, ::Nothing, ::Number) end

include("energy/standard.jl")
include("energy/isothermal.jl")
