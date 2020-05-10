export TabulatedKernel

"Tabulated SPH smoothing kernel."
struct TabulatedKernel{N} <: Kernel{N}
    support::Float
    q_to_idx_scale::Float
    W_lookup::Vector{Float}
    ddr_W_lookup::Vector{Float}
    ddh_W_lookup::Vector{Float}

    function TabulatedKernel(
        W::KernelFormula{N};
        number_of_values::Integer = 10000,
    ) where {N}
        @check number_of_values > 1

        support = get_support(W)

        q_to_idx_scale = (number_of_values - 1) / support
        idx_to_q = 1 / q_to_idx_scale

        W_lookup = [W(i * idx_to_q, 1.0) for i = 0:number_of_values-1]
        ddr_W_lookup = [ddr(W, i * idx_to_q, 1.0) for i = 0:number_of_values-1]
        ddh_W_lookup = [ddh(W, i * idx_to_q, 1.0) for i = 0:number_of_values-1]

        new{N}(support, q_to_idx_scale, W_lookup, ddr_W_lookup, ddh_W_lookup)
    end
end

Base.show(io::IO, ::Type{TabulatedKernel{N}}) where {N} =
    print(io, "TabulatedKernel")

"""
    (W::TabulatedKernel)(q::Number) -> Float

Evaluate the kernel with unit width at relative distance `q`.
"""
(W::TabulatedKernel)(q::Number) = W.W_lookup[lookupidx(W, q)]

"""
    ddr(W::TabulatedKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `r`.
"""
ddr(W::TabulatedKernel, q::Number) = W.ddr_W_lookup[lookupidx(W, q)]

"""
    ddh(W::TabulatedKernel, q::Number) -> Float

Evaluate the derivative of the kernel with unit width at relative distance `q`.
The derivative is with respect to `h`.
"""
ddh(W::TabulatedKernel{N}, q::Number) where {N} =
    W.ddh_W_lookup[lookupidx(W, q)]

get_support(W::TabulatedKernel) = W.support

lookupidx(W::TabulatedKernel, q::Number)::UInt =
    @dbgassert(0 ≤ q ≤ get_support(W)) && 1 + trunc(UInt, q * W.q_to_idx_scale)
