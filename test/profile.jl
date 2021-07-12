using Whirl
using Profile
using ProfileView

Profile.init(delay = 1e-2)

stepper = build(sod_shock_tube_builder(number_of_particles = 1000))

function profiling(n)
   for _ = 1:n
      step!(stepper, 1e-4)
   end
end

@profview profiling(1)
@profview profiling(10000)
