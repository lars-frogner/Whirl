export StandardDiffusion

mutable struct StandardDiffusion <: DiffusionComponent
    "Scale for linear term."
    α::Float
    "Scale for quadratic term."
    β::Float
    "Denominator term for avoiding division by zero."
    ε::Float
    "Squared characteristic sound speed."
    cs²::Float

    function StandardDiffusion(; α = 1.0, β = 2.0, ε = 1e-3)
        @check α ≥ 0
        @check β ≥ 0
        @check ε ≥ 0
        new(α, β, ε, 0.0)
    end
end

function updatediffusion!(
    standard_diffusion::StandardDiffusion,
    energy_component::EnergyComponent,
    eos::EquationOfState,
)
    standard_diffusion.cs² =
        soundspeed²_from(eos, mean(get_specific_energies(energy_component)))
end

function compute_Δa_diffusion(
    r::Number,
    vᵣ::Number,
    hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
    ρᵢ_ρⱼ::Tuple{Number,Number},
    ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
    mass_component::MassComponent,
    standard_diffusion::StandardDiffusion,
)
    m = get_particle_mass(mass_component)
    α = standard_diffusion.α
    β = standard_diffusion.β
    ε = standard_diffusion.ε
    cs² = standard_diffusion.cs²

    hᵢⱼ = mean(hᵢ_hⱼ)
    ρᵢⱼ = mean(ρᵢ_ρⱼ)
    ddr_W_r_hᵢⱼ = mean(ddr_W_r_hᵢ_hⱼ)
    μᵢⱼ = vᵣ * r * hᵢⱼ / (r^2 + ε * hᵢⱼ^2)
    Πᵢⱼ = μᵢⱼ < 0 ? (β * μᵢⱼ - α * cs²) * μᵢⱼ / ρᵢⱼ : 0.0
    m * Πᵢⱼ * ddr_W_r_hᵢⱼ
end
