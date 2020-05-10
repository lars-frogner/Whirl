mutable struct TimeDerivatives{N,MD,ED}
    accelerations::Vector{SVectorF{N}}
    mass_derivatives::MD
    energy_derivatives::ED

    function TimeDerivatives(
        velocities::Velocities{N},
        mass_component::MassComponent{MD,<:Kernel{N}},
        energy_component::EnergyComponent{ED},
    ) where {N,MD,ED}
        new{N,MD,ED}(
            initderivatives(velocities),
            initderivatives(mass_component),
            initderivatives(energy_component),
        )
    end
end

"""
    initderivatives!(
        derivatives::TimeDerivatives{N,MD,ED},
        velocities::Velocities{N},
        mass_component::MassComponent{MD,<:Kernel{N}},
        energy_component::EnergyComponent{ED},
    )

Resets the given fluid derivatives.
"""
function initderivatives!(
    derivatives::TimeDerivatives{N,MD,ED},
    velocities::Velocities{N},
    mass_component::MassComponent{MD,<:Kernel{N}},
    energy_component::EnergyComponent{ED},
) where {N,MD,ED}
    initderivatives!(derivatives.accelerations, velocities)
    initderivatives!(derivatives.mass_derivatives, mass_component)
    initderivatives!(derivatives.energy_derivatives, energy_component)
end

"""
    updatefluid!(
        time_derivatives::TimeDerivatives{N,MD,ED},
        mass_component::MassComponent{MD,<:Kernel{N}},
        energy_component::EnergyComponent{ED},
        diffusion_component::DiffusionComponent,
        pressures::Pressures,
        velocities::Velocities{N},
        positions::Vector{SVectorF{N}},
        eos::EquationOfState,
    )

Update the state of the fluid variables and derivatives according to the
current particle positions.
"""
function updatefluid!(
    time_derivatives::TimeDerivatives{N,MD,ED},
    mass_component::MassComponent{MD,<:Kernel{N}},
    energy_component::EnergyComponent{ED},
    diffusion_component::DiffusionComponent,
    pressures::Pressures,
    velocities::Velocities{N},
    positions::Vector{SVectorF{N}},
    eos::EquationOfState,
) where {N,MD,ED}
    updatemasses!(mass_component, positions)
    updatepressures!(pressures, mass_component, energy_component, eos)
    updatediffusion!(diffusion_component, energy_component, eos)

    initderivatives!(
        time_derivatives,
        velocities,
        mass_component,
        energy_component,
    )

    n::UInt = length(positions)
    for i = 1:n-1, j = i:n
        rᵢⱼ = positions[j] - positions[i]
        r = norm(rᵢⱼ)
        hᵢ_hⱼ = get_kernel_widths(mass_component, (i, j))
        qᵢ_qⱼ = normalize_distance(mass_component, r, hᵢ_hⱼ)

        if r > 0 && nonzeroat(get_kernel(mass_component), qᵢ_qⱼ)
            r̂ᵢⱼ = rᵢⱼ / r
            updatederivatives!(
                time_derivatives,
                (i, j),
                (rᵢⱼ, r̂ᵢⱼ, r),
                hᵢ_hⱼ,
                qᵢ_qⱼ,
                velocities,
                pressures,
                mass_component,
                energy_component,
                diffusion_component,
            )
        end
    end
end

"""
    updatederivatives!(
        time_derivatives::TimeDerivatives{N,MD,ED},
        (i, j)::Tuple{Unsigned,Unsigned},
        (rᵢⱼ, r̂ᵢⱼ, r)::Tuple{SVectorF{N},SVectorF{N},Number},
        hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
        qᵢ_qⱼ::Union{Number,Tuple{Number,Number}},
        velocities::Velocities{N},
        pressures::Pressures,
        mass_component::MassComponent{MD,<:Kernel{N}},
        energy_component::EnergyComponent{ED},
        diffusion_component::DiffusionComponent,
    )

Update time derivatives for particles `i` and `j`.
"""
function updatederivatives!(
    time_derivatives::TimeDerivatives{N,MD,ED},
    i_j::Tuple{Unsigned,Unsigned},
    (rᵢⱼ, r̂ᵢⱼ, r)::Tuple{SVectorF{N},SVectorF{N},Number},
    hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
    qᵢ_qⱼ::Union{Number,Tuple{Number,Number}},
    velocities::Velocities{N},
    pressures::Pressures,
    mass_component::MassComponent{MD,<:Kernel{N}},
    energy_component::EnergyComponent{ED},
    diffusion_component::DiffusionComponent,
) where {N,MD,ED}
    vᵢ, vⱼ = get_velocities(velocities, i_j)
    vᵢⱼ = vⱼ - vᵢ
    vᵣ = vᵢⱼ ⋅ r̂ᵢⱼ

    ρᵢ_ρⱼ = get_mass_densities(mass_component, i_j)

    ddr_W_r_hᵢ_hⱼ = ddr(get_kernel(mass_component), qᵢ_qⱼ, hᵢ_hⱼ)

    m_ddr_W_r_hᵢ_hⱼ = compute_m_ddr_W_r_h(i_j, ddr_W_r_hᵢ_hⱼ, mass_component)

    Δaᵢ_pressure, Δaⱼ_pressure =
        compute_Δa_pressure(i_j, ρᵢ_ρⱼ, m_ddr_W_r_hᵢ_hⱼ, pressures)

    Δa_diffusion = compute_Δa_diffusion(
        r,
        vᵣ,
        hᵢ_hⱼ,
        ρᵢ_ρⱼ,
        ddr_W_r_hᵢ_hⱼ,
        mass_component,
        diffusion_component,
    )

    updatederivatives!(
        time_derivatives.mass_derivatives,
        i_j,
        vᵣ,
        m_ddr_W_r_hᵢ_hⱼ,
        mass_component,
    )
    updatederivatives!(
        time_derivatives.energy_derivatives,
        i_j,
        vᵣ,
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion),
        energy_component,
    )
    updatederivatives!(
        time_derivatives.accelerations,
        i_j,
        r̂ᵢⱼ,
        (Δaᵢ_pressure, Δaⱼ_pressure, Δa_diffusion),
        velocities,
    )
end

"""
    evolvevariables!(
        velocities::Velocities{N},
        mass_component::MassComponent{MD},
        energy_component::EnergyComponent{ED},
        time_derivatives::TimeDerivatives{N,MD,ED},
        Δt::Number,
    )

Evolve the components with the given derivatives and time step.
"""
function evolvevariables!(
    velocities::Velocities{N},
    mass_component::MassComponent{MD},
    energy_component::EnergyComponent{ED},
    time_derivatives::TimeDerivatives{N,MD,ED},
    Δt::Number,
) where {N,MD,ED}
    evolvevariables!(velocities, time_derivatives.accelerations, Δt)
    evolvevariables!(mass_component, time_derivatives.mass_derivatives, Δt)
    evolvevariables!(energy_component, time_derivatives.energy_derivatives, Δt)
end

"""
    evolvevariables!(
        velocities::Velocities{N},
        mass_component::MassComponent{MD},
        energy_component::EnergyComponent{ED},
        Δt::Number,
        time_derivatives_and_weights::Tuple{TimeDerivatives{N,MD,ED},Float}...,
    )

Evolve the components with the given time step and weighed sets of derivatives.
"""
function evolvevariables!(
    velocities::Velocities{N},
    mass_component::MassComponent{MD},
    energy_component::EnergyComponent{ED},
    Δt::Number,
    time_derivatives_and_weights::Tuple{TimeDerivatives{N,MD,ED},Float}...,
) where {N,MD,ED}
    for (time_derivatives, weight) in time_derivatives_and_weights
        weighted_Δt = weight * Δt
        evolvevariables!(
            velocities,
            mass_component,
            energy_component,
            time_derivatives,
            weighted_Δt,
        )
    end
end

"""
    evolvepositions!(
        positions::Vector{SVectorF{N}},
        Δt::Number,
        velocities_and_weights::Tuple{Velocities{N},Number}...,
    )

Evolves the positions with the given time step weighed sets of velocities.
"""
function evolvepositions!(
    positions::Vector{SVectorF{N}},
    Δt::Number,
    velocities_and_weights::Tuple{Velocities{N},Number}...,
) where {N}
    for (velocities, weight) in velocities_and_weights
        weighted_Δt = weight * Δt
        evolvepositions!(positions, velocities, weighted_Δt)
    end
end
