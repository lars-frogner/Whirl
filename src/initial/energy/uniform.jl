export UniformEnergyDistribution

"Initial uniform energy distribution."
struct UniformEnergyDistribution{N} <: InitialEnergyDistribution{N}
    "Specific energy."
    u::Float

    UniformEnergyDistribution{N}(u::Number) where {N} = @check(u > 0) && new(u)
end

Base.show(io::IO, ::Type{UniformEnergyDistribution{N}}) where {N} =
    print(io, "UniformEnergyDistribution")

function (distribution::UniformEnergyDistribution{N})(
    positions::AbstractVector{SVectorF{N}},
    mass_densities::AbstractVector{Float},
    eos::EquationOfState,
) where {N}
    @checkeq(length(positions), length(mass_densities))
    fill(distribution.u, length(positions))
end
