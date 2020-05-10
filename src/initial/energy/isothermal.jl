export IsothermalEnergyDistribution

"Initial isothermal energy distribution."
struct IsothermalEnergyDistribution{N} <: InitialEnergyDistribution{N}
    "Temperature."
    T::Float

    IsothermalEnergyDistribution{N}(T::Number) where {N} =
        @check(T > 0) && new(T)
end

Base.show(io::IO, ::Type{IsothermalEnergyDistribution{N}}) where {N} =
    print(io, "IsothermalEnergyDistribution")

function (distribution::IsothermalEnergyDistribution{N})(
    positions::AbstractVector{SVectorF{N}},
    mass_densities::AbstractVector{Float},
    eos::EquationOfState,
) where {N}
    @checkeq(length(positions), length(mass_densities))
    fill(specific_energy_from(eos, distribution.T), length(positions))
end
