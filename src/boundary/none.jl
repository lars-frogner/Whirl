export NoBoundaries

struct NoBoundaries{N} <: BoundaryComponent{N} end

Base.show(io::IO, ::Type{NoBoundaries{N}}) where {N} = print(io, "NoBoundaries")

get_positions(::NoBoundaries{N}) where {N} = SVectorF{N}[]
get_velocities(::NoBoundaries{N}) where {N} = SVectorF{N}[]
get_pressures(::NoBoundaries) = Float[]

function compute_mass_density(
    ::NoBoundaries,
    ::EnergyComponent,
    ::EquationOfState,
    ::Unsigned,
    ::Unsigned,
)
    nothing
end

function updateboundaries!(
    ::NoBoundaries,
    ::AbstractVector{SVectorF{N}},
    ::Velocities{N},
    ::MassComponent{MD,<:Kernel{N}},
    ::Pressures,
) where {N,MD} end

function evolveboundaries!(::NoBoundaries, ::Number) end
