abstract type InitialMassDistribution{N} end

struct MassDistributionData{N}
    "Mass of the fluid particles."
    m::Float
    "Positions of the fluid particles."
    positions::Vector{SVectorF{N}}
    "Fluid particle mass densities."
    mass_densities::Vector{Float}

    function MassDistributionData(
        m::Number,
        positions::AbstractVector{SVectorF{N}},
        mass_densities::AbstractVector{Float},
    ) where {N}
        new{N}(m, positions, mass_densities)
    end
end

"""
    (::InitialMassDistribution{N})() -> MassDistributionData{N}

Compute the mass distribution data corresponding to this initial mass
distribition.
"""
function (::InitialMassDistribution) end

include("mass/empty.jl")
include("mass/uniform.jl")
include("mass/piecewise_uniform.jl")
