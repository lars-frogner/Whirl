mutable struct Velocities{N}
    velocities::Vector{SVectorF{N}}

    function Velocities(
        positions::AbstractVector{SVectorF{N}},
        mass_densities::AbstractVector{Float},
        velocity_distribution::InitialVelocityDistribution{N},
    ) where {N}
        velocities = velocity_distribution(positions, mass_densities)
        new{N}(velocities)
    end
end

"""
    get_velocities(velocities::Velocities{N})
        -> Vector{SVectorF{N}}

Return the velocities of the particles.
"""
get_velocities(velocities::Velocities) = velocities.velocities

"""
    get_velocities(velocities::Velocities{N}, (i, j)::Tuple{Unsigned,Unsigned})
        -> (SVectorF{N}, SVectorF{N})

Return the velocities of the particles with the given indices.
"""
function get_velocities(
    velocities::Velocities,
    (i, j)::Tuple{Unsigned,Unsigned},
)
    velocities.velocities[i], velocities.velocities[j]
end

"""
    get_velocity(velocities::Velocities{N}, i::Unsigned) -> SVectorF{N}

Return the velocity of the particle with the given index.
"""
function get_velocity(velocities::Velocities, i::Unsigned)
    velocities.velocities[i]
end

"""
    initderivatives(velocities::Velocities{N})
        -> Vector{SVectorF{N}}

Creates an initial array of accelerations for the velocity component.
"""
initderivatives(velocities::Velocities{N}) where {N} =
    [zeros(SVectorF{N}) for _ in eachindex(velocities.velocities)]

"""
    initderivatives!(
        accelerations::Vector{SVectorF{N}},
        ::Velocities{N},
    )

Resets the given array of accelerations for the velocity component.
"""
function initderivatives!(
    accelerations::AbstractVector{SVectorF{N}},
    ::Velocities{N},
) where {N}
    map!(_ -> zeros(SVectorF{N}), accelerations, eachindex(accelerations))
end

"""
    updatederivatives!(
        accelerations::AbstractVector{SVectorF{N}},
        (i, j)::Tuple{Unsigned,Unsigned},
        r̂ᵢⱼ::SVectorF{N},
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
        ::Velocities{N},
    )

Updates the accelerations of particles `i` and `j`.
"""
function updatederivatives!(
    accelerations::AbstractVector{SVectorF{N}},
    (i, j)::Tuple{Unsigned,Unsigned},
    r̂ᵢⱼ::SVectorF{N},
    Δa::Tuple{Number,Number,Number},
    ::Velocities{N},
) where {N}
    Δaᵢⱼ = compute_Δaᵢⱼ(Δa, r̂ᵢⱼ)
    accelerations[i] += Δaᵢⱼ
    accelerations[j] -= Δaᵢⱼ
end

"""
    updatederivative!(
        accelerations::AbstractVector{SVectorF{N}},
        i::Unsigned,
        r̂ᵢⱼ::SVectorF{N},
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
        ::Velocities{N},
    )

Updates the acceleration of particle `i`.
"""
function updatederivative!(
    accelerations::AbstractVector{SVectorF{N}},
    i::Unsigned,
    r̂ᵢⱼ::SVectorF{N},
    Δa::Tuple{Number,Number,Number},
    ::Velocities{N},
) where {N}
    accelerations[i] += compute_Δaᵢⱼ(Δa, r̂ᵢⱼ)
end

compute_Δaᵢⱼ(
    (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion)::Tuple{Number,Number,Number},
    r̂ᵢⱼ::SVectorF{N},
) where {N} = (Δaᵢ_pressure + Δaⱼ_pressure + Δa_diffusion) * r̂ᵢⱼ

"""
    evolvevariables!(
        velocities::Velocities{N},
        accelerations::AbstractVector{SVectorF{N}},
        Δt::Number,
    )

Evolves the velocities with the given accelerations and time step.
"""
evolvevariables!(
    velocities::Velocities{N},
    accelerations::AbstractVector{SVectorF{N}},
    Δt::Number,
) where {N} = velocities.velocities .+= Δt * accelerations

"""
    evolvepositions!(
        positions::AbstractVector{SVectorF{N}},
        velocities::Velocities{N},
        Δt::Number,
    )

Evolves the positions with the given velocities and time step.
"""
evolvepositions!(
    positions::AbstractVector{SVectorF{N}},
    velocities::Velocities{N},
    Δt::Number,
) where {N} = positions .+= Δt * velocities.velocities
