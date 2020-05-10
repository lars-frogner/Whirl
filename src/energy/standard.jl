export StandardEnergy

mutable struct StandardEnergy <: EnergyComponent{Vector{Float}}
    specific_energies::Vector{Float}

    """
        StandardEnergy(
            positions::AbstractVector{SVectorF{N}},
            mass_densities::AbstractVector{Float},
            energy_distribution::InitialEnergyDistribution{N},
            eos::EquationOfState,
        )
    """
    function StandardEnergy(
        positions::AbstractVector{SVectorF{N}},
        mass_densities::AbstractVector{Float},
        energy_distribution::InitialEnergyDistribution{N},
        eos::EquationOfState,
    ) where {N}
        specific_energies = energy_distribution(positions, mass_densities, eos)
        new(specific_energies)
    end
end

get_specific_energies(standard_energy::StandardEnergy) =
    standard_energy.specific_energies

initderivatives(standard_energy::StandardEnergy) =
    zeros(length(standard_energy.specific_energies))

initderivatives!(dudt::AbstractVector{Float}, ::StandardEnergy) =
    fill!(dudt, 0.0)

function updatederivatives!(
    dudt::AbstractVector{Float},
    (i, j)::Tuple{Unsigned,Unsigned},
    vᵣ::Number,
    (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
    ::StandardEnergy,
)
    Δdudtᵢ = (Δaᵢ_pressure + Δa_diffusion) * vᵣ
    Δdudtⱼ = (Δaⱼ_pressure + Δa_diffusion) * vᵣ
    dudt[i] += Δdudtᵢ
    dudt[j] += Δdudtⱼ
end

evolvevariables!(
    standard_energy::StandardEnergy,
    dudt::AbstractVector{Float},
    Δt::Number,
) = standard_energy.specific_energies .+= Δt * dudt
