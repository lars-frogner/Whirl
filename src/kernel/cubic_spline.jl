export CubicSplineKernel, TabulatedCubicSplineKernel

"Cubic spline smoothing kernel."
struct CubicSplineKernel{N} <: KernelFormula{N}
    σ::Float

    CubicSplineKernel{1}() = new(2 / 3)
    CubicSplineKernel{2}() = new(10 / (7π))
    CubicSplineKernel{3}() = new(1 / π)
end

Base.show(io::IO, ::Type{CubicSplineKernel{N}}) where {N} =
    print(io, "CubicSplineKernel")

"""
    (W::CubicSplineKernel)(q::Number) -> Float

Evaluate the kernel with unit width at relative distance `q`.
"""
(W::CubicSplineKernel)(q::Number) =
    @dbgassert(0 ≤ q ≤ 2) &&
    W.σ * (q ≤ 1 ? (1 - 1.5 * q^2 * (1 - 0.5 * q)) : (0.25 * (2 - q)^3))

"""
    ddr(W::CubicSplineKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `r`.
"""
ddr(W::CubicSplineKernel, q::Number) =
    @dbgassert(0 ≤ q ≤ 2) &&
    W.σ * (q ≤ 1 ? q * (2.25 * q - 3) : -0.75 * (2 - q)^2)

"""
    ddh(W::CubicSplineKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `h`.
"""
ddh(W::CubicSplineKernel{N}, q::Number) where {N} =
    @dbgassert(0 ≤ q ≤ 2) &&
    W.σ * (
        q ≤ 1 ? ((3 + 1.5N) * q^2 - (2.25 + 0.75N) * q^3 - N) :
        (-0.5 * (2 - q)^2 * (N - (1.5 + 0.5 * N) * q))
    )

get_support(W::CubicSplineKernel) = 2.0
