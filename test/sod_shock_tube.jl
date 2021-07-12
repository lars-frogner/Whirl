using Whirl
using Plots

gr()

Δt = 1e-4
n_steps = 10000
n_skips = 50

builder = @show sod_shock_tube_builder(number_of_particles = 1000)
stepper = build(builder)

anim = @animate for _ = 1:n_steps
    sphplot(get_positions(stepper), get_mass_densities(stepper))
    step!(stepper, Δt)
end every n_skips
gif(anim, "sod_shock_tube.gif", fps = 15)
