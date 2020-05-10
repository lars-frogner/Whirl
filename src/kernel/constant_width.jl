export ConstantWidthKernel

"SPH smoothing kernel with constant width."
struct ConstantWidthKernel{N,K} <: Kernel{N}
    W::K
    h::Float
    h⁻¹::Float
    h⁻ᴺ::Float
    h⁻ᴺ⁻¹::Float

    ConstantWidthKernel(W::K, h::Number) where {N,K<:Kernel{N}} =
        @check(h > 0) && new{N,K}(W, h, 1 / h, 1 / h^N, 1 / h^(N + 1))
end

"""
    (W::ConstantWidthKernel)(q::Number) -> Float

Evaluate the kernel with unit width at relative distance `q`.
"""
(wrapper::ConstantWidthKernel)(q::Number) = wrapper.W(q)

"""
    ddr(W::ConstantWidthKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `r`.
"""
ddr(wrapper::ConstantWidthKernel, q::Number) = ddr(wrapper.W, q)

"""
    ddh(W::ConstantWidthKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `h`.
"""
ddh(wrapper::ConstantWidthKernel, q::Number) = ddh(wrapper.W, q)

"""
    (W::ConstantWidthKernel)(q::Number, ::Number) -> Float

Evaluate the kernel with constant width `h` at relative distance `q = r/h`.
"""
(W::ConstantWidthKernel)(q::Number, ::Number) = W(q) * W.h⁻ᴺ

"""
    ddr(W::ConstantWidthKernel, q::Number, ::Number) -> Float

Evaluate the derivative of the kernel with constant width `h` at relative
distance `q = r/h`. The derivative is with respect to `r`.
"""
ddr(W::ConstantWidthKernel, q::Number, ::Number) = ddr(W, q) * W.h⁻ᴺ⁻¹

"""
    ddh(W::ConstantWidthKernel, q::Number, ::Number) -> Float

Evaluate the derivative of the kernel with constant width `h` at relative
distance `q = r/h`. The derivative is with respect to `h`.
"""
ddh(W::ConstantWidthKernel, q::Number, ::Number) = ddh(W, q) * W.h⁻ᴺ⁻¹

get_support(wrapper::ConstantWidthKernel) = get_support(wrapper.W)

get_kernel_width(wrapper::ConstantWidthKernel) = wrapper.h

get_inverse_kernel_width(wrapper::ConstantWidthKernel) = wrapper.h⁻¹
