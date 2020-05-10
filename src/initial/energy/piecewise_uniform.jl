export PiecewiseUniformEnergyDistribution

"Initial piecewise uniform energy distribution."
struct PiecewiseUniformEnergyDistribution <: InitialEnergyDistribution{1}
    interior_boundary_positions::Vector{Float}
    specific_energies::Vector{Float}

    function PiecewiseUniformEnergyDistribution(
        interior_boundary_positions::AbstractVector{Float},
        specific_energies::AbstractVector{Float},
    )
        @checkeq(
            length(specific_energies),
            length(interior_boundary_positions) + 1
        )
        @check all(diff(interior_boundary_positions) .> 0)
        @check all(specific_energies .> 0)
        new(interior_boundary_positions, specific_energies)
    end
end

function (distribution::PiecewiseUniformEnergyDistribution)(
    positions::AbstractVector{SVectorF{1}},
    mass_densities::AbstractVector{Float},
    eos::EquationOfState,
)
    @checkeq(length(positions), length(mass_densities))
    [specific_energy_at(distribution, r[]) for r in positions]
end

"""
    specific_energy_at(distr::PiecewiseUniformEnergyDistribution, r) -> Float

Return the energy density of the piecewise uniform energy distribution
at the given position `r`.
"""
specific_energy_at(distribution::PiecewiseUniformEnergyDistribution, r) =
    distribution.specific_energies[1+searchsortedlast(
        distribution.interior_boundary_positions,
        r,
    )]
