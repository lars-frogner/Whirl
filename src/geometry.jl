export DIM1, DIM2, DIM3, Bounds

SVectorF{N} = SVector{N,Float}

struct Dim{N}
    Dim{1}() = new()
    Dim{2}() = new()
    Dim{3}() = new()
end

const DIM1 = Dim{1}()
const DIM2 = Dim{2}()
const DIM3 = Dim{3}()

struct Bounds{N,T}
    lower::SVector{N,T}
    upper::SVector{N,T}

    function Bounds(
        lower::SVector{N,T},
        upper::SVector{N,T},
    ) where {N,T<:Number}
        @check all(upper > lower)
        new{N,T}(lower, upper)
    end

    function Bounds(lower::T, upper::T) where {T<:Number}
        @check upper > lower
        new{1,T}(SVector(lower), SVector(upper))
    end
end

lower(bounds::Bounds) = bounds.lower

upper(bounds::Bounds) = bounds.upper

"""
    extents(bounds::Bounds{N,T}) => SVector{N,T}

Compute the extents of the region within the given bounds.
"""
extents(bounds::Bounds) = bounds.upper - bounds.lower

"""
    extended(bounds::Bounds{N,T}, lengths::SVector{N,T}) => Bounds{N,T}

Compute new bounds extended symmetrically by the given lengths.
"""
function extended(bounds::Bounds{N,T}, lengths::SVector{N,T}) where {N,T}
    half_lengths = 0.5 * lengths
    Bounds(bounds.lower - half_lengths, bounds.upper + half_lengths)
end
