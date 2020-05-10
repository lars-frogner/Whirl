abstract type InitialEnergyDistribution{N} end

"""
    (::InitialEnergyDistribution{N})(
        positions::Vector{SVectorF{N}},
        mass_densities::Vector{Float},
        eos::EquationOfState,
    ) -> Vector{Float}

Compute the fluid particle specific energies corresponding to this initial energy
distribution.
"""
function (::InitialEnergyDistribution) end

include("energy/uniform.jl")
include("energy/piecewise_uniform.jl")
include("energy/isothermal.jl")
