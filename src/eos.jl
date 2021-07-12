export energy_density_from,
    gas_pressure_from, temperature_from, mass_density_from, specific_energy_from

abstract type EquationOfState end

"""
    energy_density_from(eos::EquationOfState, ρ::Number, u::Number) -> Float

Compute the internal energy density for the given mass density `ρ` and specific
energy `u`."
"""
function energy_density_from end

"""
    gas_pressure_from(eos::EquationOfState, ρ::Number, u::Number) -> Float

Compute the gas pressure for the given mass density `ρ` and internal specific
energy `u`."
"""
function gas_pressure_from end

"""
    temperature_from(eos::EquationOfState, u::Number) -> Float

Compute the temperature for the given internal specific energy `u`."
"""
function temperature_from end

"""
    mass_density_from(eos::EquationOfState, P::Number, u::Number) -> Float

Compute the mass density for the given gas pressure `P` and specific energy `u`."
"""
function mass_density_from end

"""
    specific_energy_from(eos::EquationOfState, T::Number) -> Float

Compute the internal specific energy for the given temperature `T`."
"""
function specific_energy_from end

"""
    specific_energy_from(eos::EquationOfState, ρ::Number, P::Number) -> Float

Compute the internal specific energy for the given mass density `ρ` and
gas pressure `P`."
"""
function specific_energy_from end

"""
    soundspeed²_from(eos::EquationOfState, u::Number) -> Float

Compute the squared sound speed from the given internal specific energy `u`."
"""
function soundspeed²_from end

include("eos/ideal_gas.jl")
