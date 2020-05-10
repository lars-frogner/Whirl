export SI

"Constants written in the same system of physical units."
struct UnitSystem
    name::Symbol
    "Boltzmann constant."
    k::Float
    "Atomic mass unit."
    AMU::Float

    UnitSystem(; name::Symbol, k::Number, AMU::Number) = new(name, k, AMU)
end

"Constants written in SI units"
const SI = UnitSystem(
    name = :SI,
    k = 1.38064852e-23, # [m^2 kg/(s^2 K)]
    AMU = 1.66054e-27,  # [kg]
)
