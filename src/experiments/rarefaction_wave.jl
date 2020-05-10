export rarefaction_wave_builder

function rarefaction_wave_builder(;
    width::Number = 1.0,
    ρ::Number = 1.0,
    u::Number = 2.0,
    number_of_particles::Integer = 200,
    γ::Number = 1.4,
)
    exterior_boundary_positions = Bounds(-0.5 * width, 0.5 * width)
    interior_boundary_positions = Float[]
    mass_densities = [ρ]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        number_of_particles,
    )
    energy_distribution = UniformEnergyDistribution{1}(u)

    SimulationBuilder(
        mass_distribution,
        energy_distribution,
        eos = AdiabaticIdealGasEOS(γ, NORMALIZED),
        velocity_distribution = StaticVelocityDistribution{1}(),
    )
end
