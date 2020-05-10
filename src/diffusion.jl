abstract type DiffusionComponent end

"""
    updatediffusion!(
        diffusion_component::DiffusionComponent
        energy_component::EnergyComponent,
        eos::EquationOfState,
    )

Update the state of the diffusion component according to the given energy
component.
"""
function updatediffusion!(
    ::DiffusionComponent,
    ::EnergyComponent,
    ::EquationOfState,
) end

"""
    compute_Δa_diffusion(
        r::Number,
        vᵣ::Number,
        hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
        ρᵢ_ρⱼ::Tuple{Number,Number},
        ddr_W_r_hᵢ_hⱼ::Union{Number,Tuple{Number,Number}},
        mass_component::MassComponent,
        diffusion_component::DiffusionComponent,
    ) -> Float

Compute pairwise particle acceleration due to diffusion.
"""
compute_Δa_diffusion(
    ::Number,
    ::Number,
    ::Union{Number,Tuple{Number,Number}},
    ::Tuple{Number,Number},
    ::Union{Number,Tuple{Number,Number}},
    ::MassComponent,
    ::DiffusionComponent,
) = 0.0

include("diffusion/none.jl")
include("diffusion/standard.jl")
