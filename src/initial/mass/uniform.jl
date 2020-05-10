export UniformMassDistribution

"Initial uniform mass distribution."
struct UniformMassDistribution{N} <: InitialMassDistribution{N}
    shape::SVector{N,Int}
    lower_bounds::SVectorF{N}
    cell_extent::Float
    ρ::Float

    """
        UniformMassDistribution(
            bounds::Bounds{N,<:Number},
            ρ::Number,
            number_of_particles::Integer,
        )

    Create a new uniform mass distribution with the given mass density `ρ`
    between the given bounds, with the given number of particles.
    """
    function UniformMassDistribution(
        bounds::Bounds{N,<:Number},
        ρ::Number,
        number_of_particles::Integer,
    ) where {N}
        @check ρ > 0
        @check number_of_particles > 0
        grid_extents = extents(bounds)
        volume = prod(grid_extents)
        cell_extent = (volume / number_of_particles)^(1 / N)
        shape =
            max.(
                ones(SVector{N,Int}),
                SVector{N,Int}(round.(grid_extents / cell_extent)),
            )
        new_extents = shape * cell_extent
        extensions = new_extents - grid_extents
        new_bounds = extended(bounds, extensions)
        new_lower_bounds = lower(new_bounds)
        new_cell_extent = (prod(new_extents) / prod(shape))^(1 / N)
        new{N}(shape, new_lower_bounds, new_cell_extent, ρ)
    end
end

Base.show(io::IO, ::Type{UniformMassDistribution{N}}) where {N} =
    print(io, "UniformMassDistribution")

function (distribution::UniformMassDistribution{N})() where {N}
    number_of_particles = prod(distribution.shape)
    m = distribution.ρ * distribution.cell_extent

    coordinates = [
        range(
            distribution.lower_bounds[dim] + 0.5 * distribution.cell_extent,
            step = distribution.cell_extent,
            length = distribution.shape[dim],
        ) for dim = 1:N
    ]
    indices = CartesianIndices(Tuple(distribution.shape))
    positions = [
        SVectorF{N}([coordinates[dim][indices[i][dim]] for dim = 1:N]) for i = 1:number_of_particles
    ]

    mass_densities = fill(distribution.ρ, number_of_particles)

    MassDistributionData(m, positions, mass_densities)
end
