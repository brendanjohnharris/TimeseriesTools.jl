import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import Normalization: NormUnion, AbstractNormalization, nansafe
import InverseFunctions: square
using Peaks
using LinearAlgebra

export findpeaks, maskpeaks!, maskpeaks, upsample, madev

function findpeaks(x::DimensionalData.AbstractDimVector, w = 1;
                   minprom = nothing,
                   maxprom = nothing,
                   strict = true, N = nothing)
    minprom isa Function && (minprom = minprom(x))
    maxprom isa Function && (maxprom = maxprom(x))
    _pks, vals = findmaxima(x, w)
    pks, proms = peakproms(_pks, x; min = minprom, max = maxprom, strict)
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

_default_lags(x::AbstractVector) = range(1, length(x) - Int(length(x) ÷ 2), step = 1)
_default_lags(x::AbstractMatrix) = range(1, size(x, 1) - Int(size(x, 1) ÷ 2), step = 1)

function madev(x::AbstractVector, lags = _default_lags(x); p = 1)
    if !issorted(lags)
        throw(ArgumentError("Lags must be sorted"))
    end
    l = length(x)
    result = similar(x, length(lags))
    @inbounds for (i, k) in enumerate(lags)
        if k >= l
            result[i] = 0.0
        else
            n_pairs = l - k
            x1 = @view x[1:n_pairs]           # x[1:n-k]
            x2 = @view x[(k + 1):(k + n_pairs)]   # x[k+1:n]

            result[i] = norm(x1 .- x2, p) / n_pairs
        end
    end
    return result
end

function madev(x::UnivariateRegular, _lags = _default_lags(x); kwargs...)
    lags = round.(Int, _lags ./ samplingperiod(x))
    return Timeseries(madev(parent(x), lags; kwargs...), _lags)
end
function madev(x::MultivariateRegular, _lags = _default_lags(x); kwargs...)
    lags = round.(Int, _lags .* samplingperiod(x))
    d = dims(x)[2:end]
    m = mapslices(x -> madev(x, lags; kwargs...), parent(x); dims = 1)
    return Timeseries(m, _lags, d...)
end
