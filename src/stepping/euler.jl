export EulerStepper

mutable struct EulerStepper{N,EOS,K,MD,M,ED,E,D,B} <: Stepper{N}
    eos::EOS
    positions::Vector{SVectorF{N}}
    velocities::Velocities{N}
    mass_component::M
    energy_component::E
    diffusion_component::D
    boundary_component::B
    pressures::Pressures
    time_derivatives::TimeDerivatives{N,MD,ED}

    function EulerStepper(
        positions::AbstractVector{SVectorF{N}},
        velocities::Velocities{N},
        mass_component::M,
        energy_component::E,
        diffusion_component::D,
        boundary_component::B,
        eos::EOS,
    ) where {
        N,
        EOS<:EquationOfState,
        K<:Kernel{N},
        MD,
        M<:MassComponent{MD,K},
        ED,
        E<:EnergyComponent{ED},
        D<:DiffusionComponent,
        B<:BoundaryComponent{N},
    }
        pressures = Pressures(mass_component, energy_component, eos)
        time_derivatives =
            TimeDerivatives(velocities, mass_component, energy_component)
        updatefluid!(
            time_derivatives,
            mass_component,
            energy_component,
            diffusion_component,
            pressures,
            velocities,
            positions,
            eos,
            boundary_component,
        )
        new{N,EOS,K,MD,M,ED,E,D,B}(
            eos,
            positions,
            velocities,
            mass_component,
            energy_component,
            diffusion_component,
            boundary_component,
            pressures,
            time_derivatives,
        )
    end
end

Base.show(
    io::IO,
    ::Type{EulerStepper{N,EOS,K,MD,M,ED,E,D,B}},
) where {N,EOS,K,MD,M,ED,E,D,B} = print(io, "EulerStepper")

function Base.show(
    io::IO,
    ::EulerStepper{N,EOS,K,MD,M,ED,E,D,B},
) where {N,EOS,K,MD,M,ED,E,D,B}
    print(
        io,
        """
        EulerStepper{
            $(N)D
            $(EOS)
            $(K)
            $(M)
            $(E)
            $(D)
            $(B)
        }""",
    )
end

get_positions(stepper::EulerStepper) = stepper.positions
get_velocities(stepper::EulerStepper) = get_velocities(stepper.velocities)
get_mass_densities(stepper::EulerStepper) =
    get_mass_densities(stepper.mass_component)
get_kernel_widths(stepper::EulerStepper) =
    get_kernel_widths(stepper.mass_component)
get_specific_energies(stepper::EulerStepper) =
    get_specific_energies(stepper.energy_component)
get_pressures(stepper::EulerStepper) = get_pressures(stepper.pressures)

function step!(stepper::EulerStepper, Δt::Number)
    evolvevariables!(
        stepper.velocities,
        stepper.mass_component,
        stepper.energy_component,
        stepper.time_derivatives,
        Δt,
    )
    evolvepositions!(stepper.positions, stepper.velocities, Δt)
    evolveboundaries!(stepper.boundary_component, Δt)
    updatefluid!(
        stepper.time_derivatives,
        stepper.mass_component,
        stepper.energy_component,
        stepper.diffusion_component,
        stepper.pressures,
        stepper.velocities,
        stepper.positions,
        stepper.eos,
        stepper.boundary_component,
    )
end
