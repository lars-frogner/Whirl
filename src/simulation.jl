export SimulationBuilder, set!, build

mutable struct SimulationBuilder{N}
    unit_system::UnitSystem
    kernel::Kernel{N}
    eos::EquationOfState
    mass_distribution::InitialMassDistribution{N}
    velocity_distribution::InitialVelocityDistribution{N}
    energy_distribution::InitialEnergyDistribution{N}
    mass_component_type::Tuple{Type,Any}
    energy_component_type::Tuple{Type,Any}
    diffusion_component_type::Tuple{Type,Any}
    stepper_type::Tuple{Type,Any}

    function SimulationBuilder(
        mass_distribution::InitialMassDistribution{N},
        energy_distribution::InitialEnergyDistribution{N};
        unit_system::UnitSystem = SI,
        kernel::Kernel{N} = TabulatedKernel(CubicSplineKernel{N}()),
        eos::EquationOfState = AdiabaticIdealGasEOS(1.4, NORMALIZED),
        velocity_distribution::InitialVelocityDistribution{N} = StaticVelocityDistribution{
            N,
        }(),
    ) where {N}
        new{N}(
            unit_system,
            kernel,
            eos,
            mass_distribution,
            velocity_distribution,
            energy_distribution,
            (AdaptiveMass, []),
            (StandardEnergy, []),
            (StandardDiffusion, []),
            (EulerStepper, []),
        )
    end
end

function Base.show(io::IO, builder::SimulationBuilder{N}) where {N}
    print(
        io,
        """
        SimulationBuilder{
            $(N)D
            $(typeof(builder.mass_distribution))
            $(typeof(builder.velocity_distribution))
            $(typeof(builder.energy_distribution))
            $(typeof(builder.eos))
            $(typeof(builder.kernel))
            $(builder.mass_component_type[1])
            $(builder.energy_component_type[1])
            $(builder.diffusion_component_type[1])
            $(builder.stepper_type[1])
        }""",
    )
end

(Base.:+)(builder::SimulationBuilder, property) = set!(builder, property)

function set!(builder::SimulationBuilder{N}, unit_system::UnitSystem) where {N}
    builder.unit_system = unit_system
    builder
end

function set!(builder::SimulationBuilder{N}, kernel::Kernel{N}) where {N}
    builder.kernel = kernel
    builder
end

function set!(builder::SimulationBuilder{N}, eos::EquationOfState) where {N}
    builder.eos = eos
    builder
end

function set!(
    builder::SimulationBuilder{N},
    mass_distribution::InitialMassDistribution{N},
) where {N}
    builder.mass_distribution = mass_distribution
    builder
end

function set!(
    builder::SimulationBuilder{N},
    velocity_distribution::InitialVelocityDistribution{N},
) where {N}
    builder.velocity_distribution = velocity_distribution
    builder
end

function set!(
    builder::SimulationBuilder{N},
    energy_distribution::InitialEnergyDistribution{N},
) where {N}
    builder.energy_distribution = energy_distribution
    builder
end

function set!(
    builder::SimulationBuilder{N},
    mass_component_type::Type{<:MassComponent};
    kwargs...,
) where {N}
    builder.mass_component_type = (mass_component_type, kwargs)
    builder
end

function set!(
    builder::SimulationBuilder{N},
    energy_component_type::Type{<:EnergyComponent};
    kwargs...,
) where {N}
    builder.energy_component_type = (energy_component_type, kwargs)
    builder
end

function set!(
    builder::SimulationBuilder{N},
    diffusion_component_type::Type{<:DiffusionComponent};
    kwargs...,
) where {N}
    builder.diffusion_component_type = (diffusion_component_type, kwargs)
    builder
end

function set!(
    builder::SimulationBuilder{N},
    stepper_type::Type{<:Stepper{N}};
    kwargs...,
) where {N}
    builder.stepper_type = (stepper_type, kwargs)
    builder
end

function build(builder::SimulationBuilder{N}) where {N}
    mass_component, positions = builder.mass_component_type[1](
        builder.kernel,
        builder.mass_distribution;
        builder.mass_component_type[2]...,
    )
    mass_densities = get_mass_densities(mass_component)
    energy_component = builder.energy_component_type[1](
        positions,
        mass_densities,
        builder.energy_distribution,
        builder.eos;
        builder.energy_component_type[2]...,
    )
    diffusion_component =
        builder.diffusion_component_type[1](builder.diffusion_component_type[2]...)
    velocities =
        Velocities(positions, mass_densities, builder.velocity_distribution)
    stepper = builder.stepper_type[1](
        positions,
        velocities,
        mass_component,
        energy_component,
        diffusion_component,
        builder.eos;
        builder.stepper_type[2]...,
    )
    stepper
end
