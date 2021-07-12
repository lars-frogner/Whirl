abstract type BoundaryComponent{N} end

"""
    get_positions(boundary_component::BoundaryComponent{N})
        -> Vector{SVectorF{N}}

Return the positions of the boundary particles.
"""
function get_positions end

"""
    get_velocities(boundary_component::BoundaryComponent{N})
        -> Vector{SVectorF{N}}

Return the velocities of the boundary particles.
"""
function get_velocities end

"""
    get_velocity(boundary_component::BoundaryComponent{N}, i::Unsigned)
        -> SVectorF{N}

Return the velocity of the boundary particle with the given index.
"""
function get_velocity(boundary_component::BoundaryComponent, i::Unsigned)
    velocities = get_boundary_velocities(boundary_component)
    velocities[i]
end

"""
    get_pressures(boundary_component::BoundaryComponent)
        -> Vector{Float}

Return the virtual pressures of the boundary particles.
"""
function get_pressures end

"""
    get_pressure(boundary_component::BoundaryComponent, i::Unsigned)
        -> Float

Return the virtual pressure of the boundary particle with the given index.
"""
function get_pressure(boundary_component::BoundaryComponent, i::Unsigned)
    pressures = get_boundary_pressures(boundary_component)
    pressures[i]
end

"""
    compute_boundary_mass_density(
        boundary_component::BoundaryComponent,
        energy_component::EnergyComponent,
        eos::EquationOfState,
        idx_fluid::Unsigned,
        idx_boundary::Unsigned,
    )

Compute the virtual mass density of the boundary particle with index `idx_boundary`
due to its interaction with the fluid particle with index `idx_fluid`.
"""
function compute_mass_density end

"""
    updateboundaries!(
        boundary_component::BoundaryComponent{N},
        positions::AbstractVector{SVectorF{N}},
        velocities::Velocities{N},
        mass_component::MassComponent{MD,<:Kernel{N}},
        pressures::Pressures,
    )

Update the state of the boundary particles according to the given fluid particles.
"""
function updateboundaries! end

"""
    evolveboundaries!(boundary_component::BoundaryComponent, Δt::Number)

Evolves the boundary particles with the given time step.
"""
function evolveboundaries! end

include("boundary/none.jl")
include("boundary/wall.jl")
