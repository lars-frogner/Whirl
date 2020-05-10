export sod_shock_tube_builder

function sod_shock_tube_builder(;
    number_of_particles::Integer = 200,
    γ::Number = 1.4,
)
    exterior_boundary_positions = Bounds(0.0, 1.0)
    interior_boundary_positions = [0.5]
    mass_densities = [1.0, 0.125]
    pressures = [1.0, 0.1]

    eos = AdiabaticIdealGasEOS(γ, NORMALIZED)
    specific_energies =
        specific_energy_from.(Ref(eos), mass_densities, pressures)

    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        number_of_particles,
    )
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )

    SimulationBuilder(
        mass_distribution,
        energy_distribution,
        eos = eos,
        velocity_distribution = StaticVelocityDistribution{1}(),
    )
end
