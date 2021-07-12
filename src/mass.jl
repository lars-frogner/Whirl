abstract type MassComponent{D,K} end

"""
    get_kernel(mass_component::MassComponent{MD,K}) -> K

Return the smoothing kernel of the mass component.
"""
function get_kernel end

"""
    get_particle_mass(mass_component::MassComponent) -> Float

Return the mass of the particles in the mass component.
"""
function get_particle_mass end

"""
    get_mass_densities(mass_component::MassComponent) -> Vector{Float}

Return the mass densities of the particles in the mass component.
"""
function get_mass_densities end

"""
    get_kernel_widths(mass_component::MassComponent) -> Vector{Float}

Return the kernel widths for the particles in the mass component.
"""
function get_kernel_widths end

"""
    get_mass_densities(
        mass_component::MassComponent,
        (i, j)::Tuple{Unsigned,Unsigned},
    ) -> (Float, Float)

Return the mass densities of the particles with the given indices.
"""
function get_mass_densities(
    mass_component::MassComponent,
    (i, j)::Tuple{Unsigned,Unsigned},
)
    mass_densities = get_mass_densities(mass_component)
    mass_densities[i], mass_densities[j]
end

"""
    get_mass_density(mass_component::MassComponent, i::Unsigned) -> Float

Return the mass density of the particle with the given index.
"""
function get_mass_density(mass_component::MassComponent, i::Unsigned)
    mass_densities = get_mass_densities(mass_component)
    mass_densities[i]
end

"""
    get_kernel_widths(
        mass_component::MassComponent,
        (i, j)::Tuple{Unsigned,Unsigned},
    ) -> (Float, Float)

Return the kernel widths of the particles with the given indices.
"""
function get_kernel_widths(
    mass_component::MassComponent,
    (i, j)::Tuple{Unsigned,Unsigned},
)
    kernel_widths = get_kernel_widths(mass_component)
    kernel_widths[i], kernel_widths[j]
end

"""
    get_kernel_width(mass_component::MassComponent, i::Unsigned) -> Float

Return the kernel width of the particles with the given index.
"""
function get_kernel_width(mass_component::MassComponent, i::Unsigned)
    kernel_widths = get_kernel_widths(mass_component)
    kernel_widths[i]
end

"""
    initderivatives(mass_component::MassComponent{MD}) -> MD

Creates initial derivatives for the mass component.
"""
initderivatives(::MassComponent{Nothing}) = nothing

"""
    initderivatives!(
        derivatives::MD,
        ::MassComponent{MD},
    )

Resets the given derivatives for the mass component.
"""
function initderivatives!(::Nothing, ::MassComponent{Nothing}) end

"""
    updatederivatives!(
        mass_derivatives::MD,
        (i, j)::Tuple{Unsigned,Unsigned},
        vᵣ::Number,
        m_ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
        mass_component::MassComponent{MD},
    )

Update mass component derivatives for particles `i` and `j`.
"""
function updatederivatives!(
    ::Nothing,
    ::Tuple{Unsigned,Unsigned},
    ::Number,
    ::Union{Number,Tuple{Number,Number}},
    ::MassComponent{Nothing},
) end

"""
    updatederivative!(
        mass_derivatives::MD,
        i::Unsigned,
        vᵣ::Number,
        m_ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
        mass_component::MassComponent{MD},
    )

Update mass component derivative for particle `i`.
"""
function updatederivative!(
    ::Nothing,
    ::Unsigned,
    ::Number,
    ::Union{Number,Tuple{Number,Number}},
    ::Number,
    ::MassComponent{Nothing},
) end

"""
    evolvevariables!(
        mass_component::MassComponent{MD},
        mass_derivatives::MD,
        Δt::Number,
    )

Evolve mass component with the given derivatives and time step.
"""
function evolvevariables!(::MassComponent{Nothing}, ::Nothing, ::Number) end

"""
    normalize_distance(
        ::MassComponent,
        r::Number,
        (hᵢ, hⱼ)::Tuple{Number,Number},
    ) -> (Float, Float)

Normalize the distance by the two kernel widths.
"""
normalize_distance(::MassComponent, r::Number, (hᵢ, hⱼ)::Tuple{Number,Number}) =
    r / hᵢ, r / hⱼ

"""
    normalize_distance(::MassComponent, r::Number, h::Number) -> Float

Normalize the distance by the kernel width.
"""
normalize_distance(::MassComponent, r::Number, h::Number) = r / h

"""
    estimate_initial_kernel_width(
        ρ::Number,
        m::Number,
        η::Number,
        ::Dim{N},
    ) -> Vector{Float}

Estimate initial widths of the SPH smoothing kernels.
"""
function estimate_initial_kernel_width(
    ρ::Number,
    m::Number,
    η::Number,
    ::Dim{N},
) where {N}
    η * (m / ρ)^(1 / N)
end

"""
    updatemasses!(
        mass_component::MassComponent,
        positions::AbstractVector{SVectorF{N}},
    )

Update the state of the mass component according to the given particle
positions.
"""
function updatemasses!(::MassComponent, ::AbstractVector{SVectorF{N}}) where {N} end

function compute_m_ddr_W_r_h(
    ::Union{Unsigned,Tuple{Unsigned,Unsigned}},
    ddr_W_r_h::Number,
    mass_component::MassComponent,
)
    get_particle_mass(mass_component) * ddr_W_r_h
end

function compute_m_ddr_W_r_h(
    ::Tuple{Unsigned,Unsigned},
    (ddr_W_r_hᵢ, ddr_W_r_hⱼ)::Tuple{Number,Number},
    mass_component::MassComponent,
)
    m = get_particle_mass(mass_component)
    m * ddr_W_r_hᵢ, m * ddr_W_r_hⱼ
end

include("mass/standard.jl")
include("mass/adaptive.jl")
