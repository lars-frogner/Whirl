module Whirl
using LinearAlgebra
using StaticArrays
using Statistics
using Plots
import Base

include("num.jl")
include("error.jl")
include("debug.jl")
include("geometry.jl")
include("units.jl")
include("eos.jl")
include("kernel.jl")
include("initial.jl")
include("velocity.jl")
include("mass.jl")
include("energy.jl")
include("pressure.jl")
include("diffusion.jl")
include("fluid.jl")
include("stepping.jl")
include("simulation.jl")
include("experiments.jl")
include("plotting.jl")

include("test.jl")

end # module
