import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import Normalization: NormUnion, AbstractNormalization, nansafe
import InverseFunctions: square
using Peaks

export UnitPower, findpeaks, maskpeaks!, maskpeaks, upsample

# function 𝑝(x::RegularTimeseries)
#     dur = duration(x)
#     if ~isnothing(unit(dur))
#         return sum(x.^2)/dur
#     else
#         @warn "No time units found for unit power normalization. Assuming seconds."
#         return sum(x.^2)/(dur*u"s")
#     end
# end

"""
    UnitPower <: AbstractNormalization

A normalization that sets the total power of a signal to unity.
"""
mutable struct UnitPower{T} <: AbstractNormalization{T}
    dims::Any
    p::NTuple{1, AbstractArray{T}}
end;
rootpower(x::RegularTimeseries) = sqrt(sum(map(square, x)) / duration(x))
unitpower(r𝑃) = Base.Fix2(/, r𝑃)
Normalization.estimators(::Type{N}) where {N <: UnitPower} = (rootpower,);
Normalization.forward(::Type{N}) where {N <: UnitPower} = unitpower

function findpeaks(x::DimensionalData.AbstractDimVector, w = 1; minprom = nothing,
                   maxprom = nothing,
                   strict = true, N = nothing)
    minprom isa Function && (minprom = minprom(x))
    maxprom isa Function && (maxprom = maxprom(x))
    _pks, vals = findmaxima(x, w)
    pks, proms = peakproms(_pks, x; minprom, maxprom, strict)
    if !isempty(pks)
        pks, widths, leftedge, rightedge = peakwidths(pks, x, proms)
        leftedge = [only(lookup(x))[ceil(Int, l)] for l in leftedge]
        rightedge = [only(lookup(x))[floor(Int, r)] for r in rightedge]
    else
        leftedge = []
        rightedge = []
    end
    idxs = indexin(pks, _pks) .|> Int
    vals = vals[idxs]
    proms = set(vals, proms)
    widths = set(vals, [l .. r for (l, r) in zip(leftedge, rightedge)])
    if !isnothing(N)
        ps = sortperm(proms; rev = true)
        vals = vals[ps[1:N]]
        widths = widths[ps[1:N]]
    end
    return vals, proms, widths
end

function findpeaks(x::DimensionalData.AbstractDimArray, args...; dims = 1, kwargs...)
    @assert length(dims) == 1
    _dims = DimensionalData.dims(x)[DimensionalData.dims(x) .!= [DimensionalData.dims(x,
                                                                                      dims)]]
    P = findpeaks.(eachslice(x; dims = _dims); kwargs...)
    return [getindex.(P, i) for i in 1:3] # vals, proms, widths
end

function maskpeaks!(y, x::DimensionalData.AbstractDimVector, args...; kwargs...)
    vals, proms, widths = findpeaks(x, args...; kwargs...)
    y .= 0
    for (i, I) in enumerate(widths)
        y[I] .= i
    end
    return y
end
function maskpeaks(x::DimensionalData.AbstractDimVector, args...; kwargs...)
    y = set(x, similar(x, Int))
    maskpeaks!(y, x, args...; kwargs...)
    return y
end

function maskpeaks(x::DimensionalData.AbstractDimArray, args...; dims = 1, kwargs...)
    @assert length(dims) == 1
    _dims = DimensionalData.dims(x)[DimensionalData.dims(x) .!= [DimensionalData.dims(x,
                                                                                      dims)]]
    y = similar(x, Int)
    maskpeaks!.(eachslice(y; dims = _dims), eachslice(x; dims = _dims), args...;
                kwargs...)
    return y
end

function upsample(d::DimensionalData.Dimension{<:RegularIndex}, factor::Number)
    rebuild(d, range(start = minimum(d), stop = maximum(d), step = step(d) / factor))
end
function upsample(d::DimensionalData.Dimension, factor)
    rebuild(d,
            range(start = minimum(d), stop = maximum(d),
                  step = mean(diff(lookup(d))) / factor))
end
