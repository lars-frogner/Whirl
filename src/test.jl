using Test

function test_geometry()
    println(extents(Bounds(SA_F64[0, -1, 0], SA_F64[1, 2, 3])))
    true
end

function test_ideal_gas_eos()
    physical_eos = @show AdiabaticIdealGasEOS(
        0.1,
        PHYSICAL;
        dim = DIM1,
        atomicity = MONATOMIC,
        units = SI,
    )
    norm_eos = @show AdiabaticIdealGasEOS(2.1, NORMALIZED)
    true
end

function test_gaussian_kernel()
    W1 = @show GaussianKernel{1}()
    W2 = @show GaussianKernel{2}()
    W3 = @show GaussianKernel{3}()
    q, h = 0.7, 0.2
    @show W1(q, h)
    @show W2(q, h)
    @show W3(q, h)
    true
end

function test_initial_piecewise_uniform_mass_distr()
    @show PiecewiseUniformMassDistribution(Bounds(0, 3), [1, 2], [1, 2, 1], 10)()
    true
end

function test_initial_static_velocity_distr()
    @show StaticVelocityDistribution{2}()(
        fill(zeros(SVectorF{2}), 10),
        fill(1.0, 10),
    )
    true
end

function test_initial_uniform_energy_distr()
    @show UniformEnergyDistribution{2}(0.3)(
        fill(zeros(SVectorF{2}), 10),
        fill(1.0, 10),
        AdiabaticIdealGasEOS(2.0, NORMALIZED),
    )
    true
end

function test_initial_piecewise_uniform_energy_distr()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )
    data = mass_distribution()
    @show energy_distribution(
        data.positions,
        data.mass_densities,
        AdiabaticIdealGasEOS(2.1, NORMALIZED),
    )
    true
end

function test_initial_isothermal_energy_distr()
    @show IsothermalEnergyDistribution{2}(0.3)(
        fill(zeros(SVectorF{2}), 10),
        fill(1.0, 10),
        AdiabaticIdealGasEOS(2.5, NORMALIZED),
    )
    true
end

function test_velocity()
    positions = fill(zeros(SVectorF{2}), 10)
    mass_densities = fill(1.0, 10)
    velocity_distribution = StaticVelocityDistribution{2}()
    velocities = @show Velocities(
        positions,
        mass_densities,
        velocity_distribution,
    )
    accelerations = @show initderivatives(velocities)
    accelerations[2] = SA_F64[10, 5]
    initderivatives!(accelerations, velocities)
    @show accelerations
    true
end

function test_adaptive_mass_component()
    mass_distribution =
        PiecewiseUniformMassDistribution(Bounds(0, 3), [1, 2], [1, 2, 1], 10)
    kernel = GaussianKernel{1}()
    mass_component, positions = @show AdaptiveMass(mass_distribution, kernel)
    @show get_kernel_widths(mass_component)
    updatemasses!(mass_component, positions, kernel)
    @show get_kernel_widths(mass_component)
    true
end

function test_standard_energy_component()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    mass_component, positions =
        AdaptiveMass(mass_distribution, GaussianKernel{1}())
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )
    energy_component = @show StandardEnergy(
        positions,
        get_mass_densities(mass_component),
        energy_distribution,
        AdiabaticIdealGasEOS(2.5, NORMALIZED),
    )
    dudt = @show initderivatives(energy_component)
    dudt[2] = 5.0
    initderivatives!(dudt, energy_component)
    @show dudt
    true
end

function test_isothermal_energy_component()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    mass_component, positions =
        AdaptiveMass(mass_distribution, GaussianKernel{1}())
    energy_distribution = IsothermalEnergyDistribution{1}(4.5)
    energy_component = @show IsothermalEnergy(
        positions,
        get_mass_densities(mass_component),
        energy_distribution,
        AdiabaticIdealGasEOS(2.5, NORMALIZED),
    )
    true
end

function test_pressure()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    eos = AdiabaticIdealGasEOS(2.5, NORMALIZED)
    mass_component, positions =
        AdaptiveMass(mass_distribution, GaussianKernel{1}())
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )
    energy_component = StandardEnergy(
        positions,
        get_mass_densities(mass_component),
        energy_distribution,
        eos,
    )
    pressures = @show Pressures(mass_component, energy_component, eos)
    updatepressures!(pressures, mass_component, energy_component, eos)
    @show pressures
    true
end

function test_standard_diffusion_component()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    eos = AdiabaticIdealGasEOS(2.5, NORMALIZED)
    mass_component, positions =
        AdaptiveMass(mass_distribution, GaussianKernel{1}())
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )
    energy_component = StandardEnergy(
        positions,
        get_mass_densities(mass_component),
        energy_distribution,
        eos,
    )
    diffusion_component = @show StandardDiffusion(α = 1.0, β = 2.0, ε = 1e-3)
    updatediffusion!(diffusion_component, energy_component, eos)
    @show diffusion_component
    true
end

function test_euler_stepper()
    exterior_boundary_positions = Bounds(0, 3)
    interior_boundary_positions = [1, 2]
    mass_densities = [1, 2, 1]
    specific_energies = [0.1, 0.2, 0.1]
    mass_distribution = PiecewiseUniformMassDistribution(
        exterior_boundary_positions,
        interior_boundary_positions,
        mass_densities,
        10,
    )
    eos = AdiabaticIdealGasEOS(2.5, NORMALIZED)
    W = GaussianKernel{1}()
    mass_component, positions = AdaptiveMass(mass_distribution, W)
    mass_densities = get_mass_densities(mass_component)
    velocity_distribution = StaticVelocityDistribution{1}()
    velocities = Velocities(positions, mass_densities, velocity_distribution)
    energy_distribution = PiecewiseUniformEnergyDistribution(
        interior_boundary_positions,
        specific_energies,
    )
    energy_component =
        StandardEnergy(positions, mass_densities, energy_distribution, eos)
    diffusion_component = StandardDiffusion()
    stepper = @show EulerStepper(
        positions,
        velocities,
        mass_component,
        energy_component,
        diffusion_component,
        W,
        eos,
    )
    step!(stepper, 0.1)
    @show stepper
    true
end

function test_simulation_builder()
    builder = @show SimulationBuilder(
        UniformMassDistribution(Bounds(SA_F64[0, 0], SA_F64[1, 1]), 1.0, 3^2),
        UniformEnergyDistribution{2}(1),
    ) + NoDiffusion
    @show build(builder)
    true
end

function runtests()

    @testset "geometry.jl" begin
        @test test_geometry()
    end

    @testset "eos.jl" begin
        @test test_ideal_gas_eos()
    end

    @testset "kernels.jl" begin
        @test test_gaussian_kernel()
    end

    @testset "initial/mass/piecewise_uniform.jl" begin
        @test test_initial_piecewise_uniform_mass_distr()
    end

    @testset "initial/velocity/static.jl" begin
        @test test_initial_static_velocity_distr()
    end

    @testset "initial/energy/uniform.jl" begin
        @test test_initial_uniform_energy_distr()
    end

    @testset "initial/energy/piecewise_uniform.jl" begin
        @test test_initial_piecewise_uniform_energy_distr()
    end

    @testset "initial/energy/isothermal.jl" begin
        @test test_initial_isothermal_energy_distr()
    end

    @testset "velocity.jl" begin
        @test test_velocity()
    end

    @testset "mass/adaptive.jl" begin
        @test test_adaptive_mass_component()
    end

    @testset "energy/standard.jl" begin
        @test test_standard_energy_component()
    end

    @testset "energy/isothermal.jl" begin
        @test test_isothermal_energy_component()
    end

    @testset "pressure.jl" begin
        @test test_pressure()
    end

    @testset "diffusion/standard.jl" begin
        @test test_standard_diffusion_component()
    end

    @testset "stepping/euler.jl" begin
        @test test_euler_stepper()
    end

    @testset "simulation.jl" begin
        @test test_simulation_builder()
    end
end
export runtests
