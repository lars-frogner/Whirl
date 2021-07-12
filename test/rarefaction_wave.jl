using Whirl
using Plots

gr()

Δt = 1e-4
n_steps = 10000
n_skips = 50

stepper = build(rarefaction_wave_builder(number_of_particles = 1000))

anim = @animate for _ = 1:n_steps
    sphplot(get_positions(stepper), get_velocities(stepper))
    step!(stepper, Δt)
end every n_skips
gif(anim, "rarefaction_wave.gif", fps = 15)
