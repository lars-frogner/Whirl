abstract type EnergyComponent{D} end

"""
    get_specific_energies(energy_component::EnergyComponent) -> Vector{Float}

Return the specific energies of the particles in the energy component.
"""
function get_specific_energies end

"""
    get_specific_energies(
        energy_component::EnergyComponent,
        (i, j)::Tuple{Unsigned,Unsigned},
    ) -> (Float, Float)

Return the specific energies of the particles with the given indices.
"""
function get_specific_energies(
    energy_component::EnergyComponent,
    (i, j)::Tuple{Unsigned,Unsigned},
)
    specific_energies = get_specific_energies(energy_component)
    specific_energies[i], specific_energies[j]
end

"""
    get_specific_energy(energy_component::EnergyComponent, i::Unsigned) -> Float

Return the specific energy of the particle with the given index.
"""
function get_specific_energy(energy_component::EnergyComponent, i::Unsigned)
    specific_energies = get_specific_energies(energy_component)
    specific_energies[i]
end

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
    updatederivative!(
        energy_derivatives::ED,
        i::Unsigned,
        vᵣ::Number,
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
        energy_component::EnergyComponent{ED},
    )

Update energy component derivative for particle `i`.
"""
function updatederivative!(
    ::Nothing,
    ::Unsigned,
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
