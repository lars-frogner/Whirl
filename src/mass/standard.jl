export StandardMass

mutable struct StandardMass{K} <: MassComponent{Vector{Float},K}
    "Smoothing kernel."
    W::K
    "Mass of the fluid particles."
    m::Float
    mass_densities::Vector{Float}

    """
        StandardMass(
            W::Kernel{N},
            mass_distribution::InitialMassDistribution{N};
            η::Number = 1.2,
        ) -> (StandardMass, Vector{SVectorF{N}})
    """
    function StandardMass(
        W::K,
        mass_distribution::InitialMassDistribution{N};
        η::Number = 1.2,
    ) where {N,K<:Kernel{N}}
        @check !isa(W, ConstantWidthKernel)
        @check η > 0
        data = mass_distribution()
        m = data.m
        mass_densities = data.mass_densities
        positions = data.positions
        @check data.m > 0
        @checkeq(length(mass_densities), length(positions))
        h = @show estimate_initial_kernel_width(
            mean(mass_densities),
            m,
            η,
            Dim{N}(),
        )
        @check(h > 0)
        constant_width_W = ConstantWidthKernel(W, h)
        new{ConstantWidthKernel{N,K}}(constant_width_W, m, mass_densities),
        positions
    end
end

get_kernel(standard_mass::StandardMass) = standard_mass.W
get_particle_mass(standard_mass::StandardMass) = standard_mass.m
get_mass_densities(standard_mass::StandardMass) = standard_mass.mass_densities
get_kernel_widths(standard_mass::StandardMass) =
    fill(get_kernel_width(standard_mass), length(standard_mass.mass_densities))

"""
    get_kernel_widths(standard_mass::StandardMass, ::Tuple{Unsigned,Unsigned})
        -> Float

Return the constant kernel width of the standard mass component.
"""
get_kernel_widths(standard_mass::StandardMass, ::Tuple{Unsigned,Unsigned}) =
    get_kernel_width(standard_mass)

"""
    get_kernel_width(standard_mass::StandardMass) -> Float

Return the constant kernel width of the standard mass component.
"""
get_kernel_width(standard_mass::StandardMass) =
    get_kernel_width(standard_mass.W)

"""
    normalize_distance(
        standard_mass::StandardMass,
        r::Number,
        ::Number,
    ) -> Float

Normalize the distance by the constant kernel width of the standard mass component.
"""
normalize_distance(standard_mass::StandardMass, r::Number, ::Number) =
    r * get_inverse_kernel_width(standard_mass.W)

initderivatives(standard_mass::StandardMass) =
    zeros(length(standard_mass.mass_densities))

initderivatives!(dρdt::AbstractVector{Float}, ::StandardMass) = fill!(dρdt, 0.0)

function updatederivatives!(
    dρdt::AbstractVector{Float},
    (i, j)::Tuple{Unsigned,Unsigned},
    vᵣ::Number,
    m_ddr_W_r_h::Number,
    ::StandardMass,
)
    Δdρdt = m_ddr_W_r_h * vᵣ
    dρdt[i] += Δdρdt
    dρdt[j] += Δdρdt
end

function updatederivatives!(
    ::AbstractVector{Float},
    ::Tuple{Unsigned,Unsigned},
    ::Number,
    ::Tuple{Number,Number},
    ::StandardMass,
)
    throw(InvalidStateException)
end

evolvevariables!(
    standard_mass::StandardMass,
    dρdt::AbstractVector{Float},
    Δt::Number,
) = standard_mass.mass_densities .+= Δt * dρdt
