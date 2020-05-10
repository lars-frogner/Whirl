export GaussianKernel

"Gaussian smoothing kernel."
struct GaussianKernel{N} <: KernelFormula{N}
    support::Float
    σ::Float

    GaussianKernel{1}(; support::Number = 3.0) =
        @check(support > 0) && new(support, 1 / √π)
    GaussianKernel{2}(; support::Number = 3.0) =
        @check(support > 0) && new(support, 1 / π)
    GaussianKernel{3}(; support::Number = 3.0) =
        @check(support > 0) && new(support, 1 / (π * √π))
end

Base.show(io::IO, ::Type{GaussianKernel{N}}) where {N} =
    print(io, "GaussianKernel")

"""
    (W::GaussianKernel)(q::Number) -> Float

Evaluate the kernel with unit width at relative distance `q`.
"""
(W::GaussianKernel)(q::Number) = @dbgassert(q ≥ 0) && W.σ * exp(-q^2)

"""
    ddr(W::GaussianKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `r`.
"""
ddr(W::GaussianKernel, q::Number) = -2q * W(q)

"""
    ddh(W::GaussianKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `h`.
"""
ddh(W::GaussianKernel{N}, q::Number) where {N} = (2 * q^2 - N) * W(q)

get_support(W::GaussianKernel) = W.support
