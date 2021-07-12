export WallBoundaries, FreeSlip, NoSlip

abstract type SlipComponent{N} end

struct FreeSlip{N} <: SlipComponent{N}
    FreeSlip(::Dim{N}, ::Integer) where {N} = new{N}()
end

Base.show(io::IO, ::Type{FreeSlip{N}}) where {N} = print(io, "FreeSlip")

mutable struct NoSlip{N} <: SlipComponent{N}
    virtual_velocities::Vector{SVectorF{N}}

    NoSlip(::Dim{N}, number_of_particles::Integer) where {N} =
        new{N}([zeros(SVectorF{N}) for _ = 1:number_of_particles])
end

Base.show(io::IO, ::Type{NoSlip{N}}) where {N} = print(io, "NoSlip")

mutable struct WallBoundaries{N,S} <: BoundaryComponent{N}
    positions::Vector{SVectorF{N}}
    velocities::Vector{SVectorF{N}}
    acceleration::SVectorF{N}
    virtual_pressures::Vector{Float}
    slip_component::S

    function WallBoundaries(
        positions::Vector{SVectorF{N}},
        velocities::Vector{SVectorF{N}},
        acceleration::SVectorF{N},
        slip_component::Type{S} = FreeSlip,
    ) where {N,S<:SlipComponent{N}}
        number_of_particles = length(positions)
        @checkeq(length(velocities), number_of_particles)
        new{N,S}(
            positions,
            velocities,
            acceleration,
            fill(0.0, number_of_particles),
            slip_component(Dim{N}(), number_of_particles),
        )
    end
end

Base.show(io::IO, ::Type{WallBoundaries{N,S}}) where {N,S} =
    print(io, "WallBoundaries{$(S)}")

get_positions(wall::WallBoundaries) = wall.positions
get_velocities(wall::WallBoundaries) = wall.velocities
get_pressures(wall::WallBoundaries) = wall.virtual_pressures

function compute_mass_density(
    wall::WallBoundaries,
    energy_component::EnergyComponent,
    eos::EquationOfState,
    idx_fluid::Unsigned,
    idx_boundary::Unsigned,
)
    mass_density_from(
        eos,
        wall.virtual_pressures[idx_boundary],
        get_specific_energy(energy_component, idx_fluid),
    )
end

function updateboundaries!(
    wall::WallBoundaries{N},
    positions::AbstractVector{SVectorF{N}},
    velocities::Velocities{N},
    standard_mass::StandardMass{ConstantWidthKernel{N}},
    pressures::Pressures,
) where {N}
    update_wall_properties!(
        wall,
        positions,
        get_velocities(velocities),
        get_mass_densities(standard_mass),
        get_pressures(pressures),
        get_kernel(standard_mass),
    )
end

function evolveboundaries!(wall::WallBoundaries, Δt::Number)
    wall.positions .+= Δt * wall.velocities
end

function update_wall_properties!(
    wall::WallBoundaries{N},
    fluid_positions::AbstractVector{SVectorF{N}},
    fluid_velocities::AbstractVector{SVectorF{N}},
    fluid_mass_densities::AbstractVector{Float},
    fluid_pressures::AbstractVector{Float},
    W::ConstantWidthKernel{N},
) where {N}
    h = get_kernel_width(W)
    h⁻¹ = get_inverse_kernel_width(W)
    for (i, (rᵢ, vᵢ)) in enumerate(zip(wall.positions, wall.velocities))
        ΣPⱼWᵢⱼ::Float = 0.0
        ΣρⱼrᵢⱼWᵢⱼ = zeros(SVectorF{N})
        ΣvⱼWᵢⱼ = zeros(SVectorF{N})
        ΣWᵢⱼ::Float = 0.0
        for (rⱼ, vⱼ, ρⱼ, Pⱼ) in zip(
            fluid_positions,
            fluid_velocities,
            fluid_mass_densities,
            fluid_pressures,
        )
            rᵢⱼ = rⱼ - rᵢ
            qᵢⱼ = norm(rᵢⱼ) * h⁻¹
            if nonzeroat(W, qᵢⱼ)
                Wᵢⱼ = W(qᵢⱼ, h)
                ΣPⱼWᵢⱼ += Pⱼ * Wᵢⱼ
                ΣρⱼrᵢⱼWᵢⱼ .+= (ρⱼ * Wᵢⱼ) * rᵢⱼ
                accumulate_virtual_velocity(
                    ΣvⱼWᵢⱼ,
                    vⱼ,
                    Wᵢⱼ,
                    wall.slip_component,
                )
                ΣWᵢⱼ += W_rᵢⱼ_h
            end
        end
        ΣWᵢⱼ⁻¹ = 1 / ΣWᵢⱼ
        wall.virtual_pressures[i] =
            (ΣPⱼWᵢⱼ - wall.acceleration ⋅ ΣρⱼrᵢⱼWᵢⱼ) * ΣWᵢⱼ⁻¹
        set_virtual_velocity(wall.slip_component, vᵢ, ΣvⱼWᵢⱼ, ΣWᵢⱼ⁻¹, i)
    end
end

function accumulate_virtual_velocity(
    ::SVectorF{N},
    ::SVectorF{N},
    ::Number,
    ::FreeSlip{N},
) where {N} end

function set_virtual_velocity(
    ::FreeSlip{N},
    ::SVectorF{N},
    ::SVectorF{N},
    ::Number,
    ::Integer,
) where {N} end

function accumulate_virtual_velocity(
    ΣvⱼWᵢⱼ::SVectorF{N},
    vⱼ::SVectorF{N},
    Wᵢⱼ::Number,
    no_slip::NoSlip{N},
) where {N}
    ΣvⱼWᵢⱼ .+= vⱼ * Wᵢⱼ
end

function set_virtual_velocity(
    no_slip::NoSlip{N},
    vᵢ::SVectorF{N},
    ΣvⱼWᵢⱼ::SVectorF{N},
    ΣWᵢⱼ⁻¹::Number,
    i::Integer,
) where {N}
    ṽᵢ = ΣvⱼWᵢⱼ * ΣWᵢⱼ⁻¹
    no_slip.virtual_velocities[i] = 2 * vᵢ - ṽᵢ
end
