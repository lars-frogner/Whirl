"SPH smoothing kernel."
abstract type Kernel{N} end

"Mathematical definition of a smoothing kernel."
abstract type KernelFormula{N} <: Kernel{N} end

"""
    (::Kernel)(q::Number, h::Number) -> Float

Evaluate the kernel with width `h` at relative distance `q = r/h`.
"""
function (W::Kernel{N})(q::Number, h::Number) where {N}
    @dbgassert(h > 0) && W(q) / h^N
end

"""
    ddr(W::Kernel, q::Number, h::Number) -> Float

Evaluate the derivative of the kernel with width `h` at relative distance
`q = r/h`. The derivative is with respect to `r`.
"""
function ddr(W::Kernel{N}, q::Number, h::Number) where {N}
    @dbgassert(h > 0) && ddr(W, q) / h^(N + 1)
end

"""
    ddh(W::Kernel, q::Number, h::Number) -> Float

Evaluate the derivative of the kernel with width `h` at relative distance
`q = r/h`. The derivative is with respect to `h`.
"""
function ddh(W::Kernel{N}, q::Number, h::Number) where {N}
    @dbgassert(h > 0) && ddh(W, q) / h^(N + 1)
end

"""
    get_support(W::Kernel) -> Float

Returns the relative distance `q = r/h` at which the kernel is truncated.
"""
function get_support end

"""
    nonzeroat(W::Kernel, q::Number) -> Bool

Whether the kernel is non-zero at relative distance `q = r/h`.
"""
nonzeroat(W::Kernel, q::Number) = q < get_support(W)

"""
    ddr(W::Kernel, (qᵢ, qⱼ)::Tuple{Number,Number}, (hᵢ, hⱼ)::Tuple{Number,Number})
        -> (Float, Float)

Evaluate the derivative of the kernel at the relative distances `qᵢ = rᵢ/hᵢ`
and `qⱼ = rⱼ/hⱼ`. The derivative is with respect to `r`.
"""
ddr(W::Kernel, (qᵢ, qⱼ)::Tuple{Number,Number}, (hᵢ, hⱼ)::Tuple{Number,Number}) =
    ddr(W, qᵢ, hᵢ), ddr(W, qⱼ, hⱼ)

"""
    nonzeroat(W::Kernel, q::Tuple{Number,Number}) -> Bool

Whether the kernel is non-zero at both `qᵢ` and `qⱼ`.
"""
nonzeroat(W::Kernel, (qᵢ, qⱼ)::Tuple{Number,Number}) =
    nonzeroat(W, qᵢ) && nonzeroat(W, qⱼ)

include("kernel/gaussian.jl")
include("kernel/cubic_spline.jl")
include("kernel/tabulated.jl")
include("kernel/constant_width.jl")
