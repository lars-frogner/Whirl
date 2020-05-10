export get_positions,
    get_velocities,
    get_mass_densities,
    get_kernel_widths,
    get_specific_energies,
    get_pressures,
    step!

abstract type Stepper{N} end

"""
    get_positions(stepper::Stepper{N}) -> Vector{SVectorF{N}}

Return the positions of the particles.
"""
function get_positions end

"""
    get_velocities(stepper::Stepper{N}) -> Vector{SVectorF{N}}

Return the velocities of the particles.
"""
function get_velocities end

"""
    get_mass_densities(stepper::Stepper) -> Vector{Float}

Return the mass densities of the particles.
"""
function get_mass_densities end

"""
    get_kernel_widths(stepper::Stepper) -> Vector{Float}

Return the widths of the particle smoothing kernels.
"""
function get_kernel_widths end

"""
    get_specific_energies(stepper::Stepper) -> Vector{Float}

Return the specific energies of the particles.
"""
function get_specific_energies end

"""
    get_pressures(stepper::Stepper) -> Vector{Float}

Return the pressures of the particles.
"""
function get_pressures end

"""
    step!(stepper::Stepper, Δt::Number)

Advance the state of the fluid with the given time step.
"""
function step! end

include("stepping/euler.jl")
