mutable struct Pressures
    pressures::Vector{Float}

    Pressures(
        mass_component::MassComponent,
        energy_component::EnergyComponent,
        eos::EquationOfState,
    ) = new([
        gas_pressure_from(eos, ρ, u)
        for
        (ρ, u) in zip(
            get_mass_densities(mass_component),
            get_specific_energies(energy_component),
        )
    ])
end

"""
    get_pressures(pressures::Pressures)
        -> Vector{Float}

Return the pressures of the particles.
"""
get_pressures(pressures::Pressures) = pressures.pressures

"""
    get_pressures(pressures::Pressures, (i, j)::Tuple{Unsigned,Unsigned})
        -> (Float, Float)

Return the pressures of the particles with the given indices.
"""
function get_pressures(pressures::Pressures, (i, j)::Tuple{Unsigned,Unsigned})
    pressures.pressures[i], pressures.pressures[j]
end

"""
    get_pressure(pressures::Pressures, i::Unsigned) -> (Float, Float)

Return the pressure of the particle with the given index.
"""
function get_pressure(pressures::Pressures, i::Unsigned)
    pressures.pressures[i]
end

"""
    updatepressures!(
        pressures::Pressures,
        mass_component::MassComponent,
        energy_component::EnergyComponent,
        eos::EquationOfState,
    )

Update the state of the pressure component according to the given mass
and energy component.
"""
updatepressures!(
    pressures::Pressures,
    mass_component::MassComponent,
    energy_component::EnergyComponent,
    eos::EquationOfState,
) = map!(
    (ρ, u) -> gas_pressure_from(eos, ρ, u),
    pressures.pressures,
    get_mass_densities(mass_component),
    get_specific_energies(energy_component),
)

"""
    compute_Δa_pressure(
        (ρᵢ, ρⱼ)::Tuple{Number,Number},
        (Pᵢ, Pⱼ)::Tuple{Number,Number},
        m_ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
    ) -> (Float, Float)

Compute pairwise particle accelerations due to pressure.
"""
function compute_Δa_pressure(
    (ρᵢ, ρⱼ)::Tuple{Number,Number},
    (Pᵢ, Pⱼ)::Tuple{Number,Number},
    m_ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
)
    compute_Δa_pressure(m_ddr_W_r_hᵢ_hⱼ, (Pᵢ / ρᵢ^2, Pⱼ / ρⱼ^2))
end

"""
    compute_Δa_pressure(
        m_ddr_W_r_h::Number,
        (Pᵢ_over_ρᵢ², Pⱼ_over_ρⱼ²)::Tuple{Number,Number},
    ) -> (Float, Float)

Compute pairwise particle accelerations due to pressure.
"""
function compute_Δa_pressure(
    m_ddr_W_r_h::Number,
    (Pᵢ_over_ρᵢ², Pⱼ_over_ρⱼ²)::Tuple{Number,Number},
)
    m_ddr_W_r_h * Pᵢ_over_ρᵢ², m_ddr_W_r_h * Pⱼ_over_ρⱼ²
end

"""
    compute_Δa_pressure(
        (m_ddr_W_r_hᵢ, m_ddr_W_r_hⱼ)::Tuple{Number,Number},
        (Pᵢ_over_ρᵢ², Pⱼ_over_ρⱼ²)::Tuple{Number,Number},
    ) -> (Float, Float)

Compute pairwise particle accelerations due to pressure.
"""
function compute_Δa_pressure(
    (m_ddr_W_r_hᵢ, m_ddr_W_r_hⱼ)::Tuple{Number,Number},
    (Pᵢ_over_ρᵢ², Pⱼ_over_ρⱼ²)::Tuple{Number,Number},
)
    m_ddr_W_r_hᵢ * Pᵢ_over_ρᵢ², m_ddr_W_r_hⱼ * Pⱼ_over_ρⱼ²
end
