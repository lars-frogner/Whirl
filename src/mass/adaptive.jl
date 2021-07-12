export AdaptiveMass

mutable struct AdaptiveMass{K} <: MassComponent{Nothing,K}
    "Smoothing kernel."
    W::K
    "Scaling for the kernel widths (typically between 1.2 and 1.5)."
    η::Float
    "Tolerance for Newton's method."
    ε::Float
    "Mass of the fluid particles."
    m::Float
    kernel_widths::Vector{Float}
    mass_densities::Vector{Float}
    correction_scales::Vector{Float}

    """
        AdaptiveMass(
            W::Kernel{N},
            mass_distribution::InitialMassDistribution{N};
            η::Number = 1.2,
            ε::Number = 1e-3,
        ) -> (AdaptiveMass, Vector{SVectorF{N}})
    """
    function AdaptiveMass(
        W::K,
        mass_distribution::InitialMassDistribution{N};
        η::Number = 1.2,
        ε::Number = 1e-3,
    ) where {N,K<:Kernel{N}}
        @check η > 0
        @check ε > 0
        data = mass_distribution()
        m = data.m
        mass_densities = data.mass_densities
        positions = data.positions
        @check data.m > 0
        @checkeq(length(mass_densities), length(positions))
        kernel_widths =
            estimate_initial_kernel_width.(mass_densities, m, η, Ref(Dim{N}()))
        correction_scales = fill(1.0, length(data.positions))
        update_particle_volumes!(
            kernel_widths,
            mass_densities,
            correction_scales,
            positions,
            m,
            W,
            η,
            ε,
        )
        new{K}(W, η, ε, m, kernel_widths, mass_densities, correction_scales),
        positions
    end
end

get_kernel(adaptive_mass::AdaptiveMass) = adaptive_mass.W
get_particle_mass(adaptive_mass::AdaptiveMass) = adaptive_mass.m
get_mass_densities(adaptive_mass::AdaptiveMass) = adaptive_mass.mass_densities
get_kernel_widths(adaptive_mass::AdaptiveMass) = adaptive_mass.kernel_widths

updatemasses!(
    adaptive_mass::AdaptiveMass{<:Kernel{N}},
    positions::AbstractVector{SVectorF{N}},
) where {N} = update_particle_volumes!(adaptive_mass, positions)

function compute_m_ddr_W_r_h(
    ::Union{Unsigned,Tuple{Unsigned,Unsigned}},
    ::Number,
    adaptive_mass::AdaptiveMass,
)
    throw(InvalidStateException)
end

function compute_m_ddr_W_r_h(
    (i, j)::Tuple{Unsigned,Unsigned},
    (ddr_W_r_hᵢ, ddr_W_r_hⱼ)::Tuple{Number,Number},
    adaptive_mass::AdaptiveMass,
)
    m = adaptive_mass.m
    Ωᵢ, Ωⱼ =
        adaptive_mass.correction_scales[i], adaptive_mass.correction_scales[j]
    m * ddr_W_r_hᵢ / Ωᵢ, m * ddr_W_r_hⱼ / Ωⱼ
end

"""
    update_particle_volumes!(
        adaptive_mass::AdaptiveMass,
        positions::AbstractVector{SVectorF{N}},
    )

Update the particle volumes of the mass component according to the given
particle positions.
"""
function update_particle_volumes!(
    adaptive_mass::AdaptiveMass{<:Kernel{N}},
    positions::AbstractVector{SVectorF{N}},
) where {N}
    update_particle_volumes!(
        adaptive_mass.kernel_widths,
        adaptive_mass.mass_densities,
        adaptive_mass.correction_scales,
        positions,
        adaptive_mass.m,
        adaptive_mass.W,
        adaptive_mass.η,
        adaptive_mass.ε,
    )
end

function update_particle_volumes!(
    kernel_widths::AbstractVector{Float},
    mass_densities::AbstractVector{Float},
    correction_scales::AbstractVector{Float},
    positions::AbstractVector{SVectorF{N}},
    m::Number,
    W::Kernel{N},
    η::Number,
    ε::Number,
) where {N}
    MAX_ITER::UInt = 10000

    @dbgasserteq(length(mass_densities), length(positions))

    for (i, (rᵢ, hᵢ)) in enumerate(zip(positions, kernel_widths))
        h, ρ, Ω, n_iter =
            solve_h_ρ_Ω_newton(hᵢ, rᵢ, positions, m, W, η, ε, MAX_ITER)
        if n_iter >= MAX_ITER
            println(
                stderr,
                "warning: particle volume computation reached max number of iterations ($(MAX_ITER))",
            )
        end
        kernel_widths[i] = h
        mass_densities[i] = ρ
        correction_scales[i] = Ω
    end
end

function solve_h_ρ_Ω_newton(
    hᵢ::Number,
    rᵢ::SVectorF{N},
    positions::AbstractVector{SVectorF{N}},
    m::Number,
    W::Kernel{N},
    η::Number,
    ε::Number,
    max_iter::Unsigned,
) where {N}
    @dbgassert hᵢ > 0
    h = hᵢ
    ρ = 0.0
    Ω = 0.0
    n_iter = 0
    while n_iter < max_iter
        ρ, Ω = compute_ρ_Ω_for_h(h, rᵢ, positions, m, W, η)
        h_perturb = clamp(h_perturbance(h, ρ, Ω, m, η, N), -0.99, 0.99)
        h *= 1 + h_perturb

        n_iter += 1
        if abs(h_perturb) < ε
            break
        end
    end
    h, ρ, Ω, n_iter
end

h_perturbance(
    h::Number,
    ρ::Number,
    Ω::Number,
    m::Number,
    η::Number,
    N::Integer,
) = (m * (η / h)^N - ρ) / (N * ρ * Ω)

function compute_ρ_Ω_for_h(
    h::Number,
    rᵢ::SVectorF{N},
    positions::AbstractVector{SVectorF{N}},
    m::Number,
    W::Kernel{N},
    η::Number,
) where {N}
    ρ = 0.0
    Ω = 0.0
    for rⱼ in positions
        q = norm(rⱼ - rᵢ) / h
        if nonzeroat(W, q)
            ρ += W(q, h)
            Ω += ddh(W, q, h)
        end
    end
    ρ *= m
    Ω = 1 + m * Ω * h / (N * ρ)
    ρ, Ω
end
