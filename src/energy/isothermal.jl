export IsothermalEnergy

mutable struct IsothermalEnergy <: EnergyComponent{Nothing}
    specific_energies::Vector{Float}

    """
        IsothermalEnergy(
            positions::AbstractVector{SVectorF{N}},
            mass_densities::AbstractVector{Float},
            energy_distribution::InitialEnergyDistribution{N},
            eos::EquationOfState,
        )
    """
    function IsothermalEnergy(
        positions::AbstractVector{SVectorF{N}},
        mass_densities::AbstractVector{Float},
        energy_distribution::InitialEnergyDistribution{N},
        eos::EquationOfState,
    ) where {N}
        specific_energies = energy_distribution(positions, mass_densities, eos)
        new(specific_energies)
    end
end

get_specific_energies(isothermal_energy::IsothermalEnergy) =
    isothermal_energy.specific_energies
initderivatives(::IsothermalEnergy) = nothing
function initderivatives!(::Nothing, ::IsothermalEnergy) end
