export StaticVelocityDistribution

"Initial velocity distribution with all particles at rest."
struct StaticVelocityDistribution{N} <: InitialVelocityDistribution{N} end

Base.show(io::IO, ::Type{StaticVelocityDistribution{N}}) where {N} =
    print(io, "StaticVelocityDistribution")

function (::StaticVelocityDistribution{N})(
    positions::AbstractVector{SVectorF{N}},
    mass_densities::AbstractVector{Float},
) where {N}
    @checkeq(length(positions), length(mass_densities))
    [zeros(SVectorF{N}) for _ in eachindex(positions)]
end
