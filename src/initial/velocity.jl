abstract type InitialVelocityDistribution{N} end

"""
    (::InitialVelocityDistribution{N})(
        positions::AbstractVector{SVectorF{N}},
        mass_densities::AbstractVector{Float},
    ) -> Vector{SVectorF{N}}

Compute the fluid particle velocities corresponding to this initial velocity
distribution.
"""
function (::InitialVelocityDistribution) end

include("velocity/static.jl")
