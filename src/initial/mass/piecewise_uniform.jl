export PiecewiseUniformMassDistribution

"Initial 1D piecewise uniform mass distribution."
struct PiecewiseUniformMassDistribution <: InitialMassDistribution{1}
    number_of_particles::UInt
    boundary_positions::Vector{Float}
    mass_densities::Vector{Float}

    """
        PiecewiseUniformMassDistribution(
            exterior_boundary_positions::Bounds{1,<:Number},
            interior_boundary_positions::AbstractVector{Float},
            mass_densities::AbstractVector{Float},
            number_of_particles::Integer,
        )

    Create a new piecewise uniform mass distribution with the given mass densities
    between the given boundary positions, with the given number of particles.
    """
    function PiecewiseUniformMassDistribution(
        exterior_boundary_positions::Bounds{1,<:Number},
        interior_boundary_positions::AbstractVector{Float},
        mass_densities::AbstractVector{Float},
        number_of_particles::Integer,
    )
        @check number_of_particles > 0
        @checkeq(
            length(mass_densities),
            length(interior_boundary_positions) + 1
        )
        boundary_positions = [
            lower(exterior_boundary_positions)[]
            interior_boundary_positions
            upper(exterior_boundary_positions)[]
        ]
        @check all(diff(boundary_positions) .> 0)
        @check all(mass_densities .> 0)
        new(number_of_particles, boundary_positions, mass_densities)
    end
end

function (distribution::PiecewiseUniformMassDistribution)()
    m =
        sum(
            distribution.mass_densities .*
            diff(distribution.boundary_positions),
        ) / distribution.number_of_particles

    positions::Vector{SVectorF{1}} = []
    mass_densities::Vector{Float} = []
    for (ρ, lower_bound, upper_bound) in zip(
        distribution.mass_densities,
        distribution.boundary_positions[1:end-1],
        distribution.boundary_positions[2:end],
    )
        extent = upper_bound - lower_bound
        local_number_of_particles = UInt(round(ρ * extent / m))
        cell_extent = extent / local_number_of_particles
        start = lower_bound + 0.5 * cell_extent
        append!(
            positions,
            map(
                x -> SA_F64[x],
                range(
                    start,
                    step = cell_extent,
                    length = local_number_of_particles,
                ),
            ),
        )
        append!(mass_densities, fill(ρ, local_number_of_particles))
    end

    MassDistributionData(m, positions, mass_densities)
end
