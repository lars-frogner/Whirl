export AdiabaticIdealGasEOS, MONATOMIC, DIATOMIC, PHYSICAL, NORMALIZED

abstract type Atomicity end
struct Monatomic <: Atomicity end
struct Diatomic <: Atomicity end
const MONATOMIC = Monatomic()
const DIATOMIC = Diatomic()

abstract type EOSType end
struct Physical <: EOSType end
struct Normalized <: EOSType end
const PHYSICAL = Physical()
const NORMALIZED = Normalized()

degrees_of_freedom(::Dim{N}, ::Monatomic) where {N} = N
degrees_of_freedom(::Dim{N}, ::Diatomic) where {N} = N + 2

"Equation of state for an ideal adiabatic gas."
struct AdiabaticIdealGasEOS <: EquationOfState
    "Adiabatic index."
    γ::Float
    "Pressure from internal energy density."
    P_from_e::Float
    "Internal energy density from mass density and temperature."
    e_from_ρT::Float

    """
        AdiabaticIdealGasEOS(
            μ::Number,
            ::Physical;
            dim::Dim = DIM3,
            atomicity::Atomicity = MONATOMIC,
            units::UnitSystem = SI,
        )

    Create a new ideal equation of state for an adiabatic gas with the given
    mean molecular mass `μ`, dimensionality, atomicity and unit system.
    """
    function AdiabaticIdealGasEOS(
        μ::Number,
        ::Physical;
        dim::Dim = DIM3,
        atomicity::Atomicity = MONATOMIC,
        units::UnitSystem = SI,
    )
        @check μ > 0
        γ = 1 + 2 / degrees_of_freedom(dim, atomicity)
        P_from_e = γ - 1
        e_from_ρT = units.k / (μ * P_from_e)
        new(γ, P_from_e, e_from_ρT)
    end

    """
        AdiabaticIdealGasEOS(γ::Number, ::Normalized)

    Create a new ideal equation of state, normalized so `P = ρT`, for an
    adiabatic gas with the given adiabatic index `γ`.
    """
    function AdiabaticIdealGasEOS(γ::Number, ::Normalized)
        @check γ > 0
        P_from_e = γ - 1
        e_from_ρT = 1 / P_from_e
        new(γ, P_from_e, e_from_ρT)
    end
end

energy_density_from(eos::AdiabaticIdealGasEOS, ρ::Number, u::Number) =
    @dbgassert(ρ ≥ 0) && @dbgassert(u ≥ 0) && ρ * u

gas_pressure_from(eos::AdiabaticIdealGasEOS, ρ::Number, u::Number) =
    eos.P_from_e * energy_density_from(eos, ρ, u)

temperature_from(eos::AdiabaticIdealGasEOS, u::Number) =
    @dbgassert(u ≥ 0) && u / eos.e_from_ρT

mass_density_from(eos::AdiabaticIdealGasEOS, P::Number, T::Number) =
    @dbgassert(P ≥ 0) &&
    @dbgassert(T > 0) && P / (eos.e_from_ρT * eos.P_from_e * T)

specific_energy_from(eos::AdiabaticIdealGasEOS, T::Number) =
    @dbgassert(T ≥ 0) && eos.e_from_ρT * T

specific_energy_from(eos::AdiabaticIdealGasEOS, ρ::Number, P::Number) =
    @dbgassert(P ≥ 0) && @dbgassert(ρ > 0) && P / (eos.P_from_e * ρ)

soundspeed²_from(eos::AdiabaticIdealGasEOS, u::Number) =
    @dbgassert(u ≥ 0) && eos.γ * eos.P_from_e * u
