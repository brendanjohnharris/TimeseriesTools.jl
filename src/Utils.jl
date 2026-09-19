import DimensionalData.Dimensions: At, Near, Dimension
import Normalization: NormUnion, AbstractNormalization, nansafe
using Peaks
using LinearAlgebra

export findpeaks, maskpeaks!, maskpeaks, madev

"""
    findpeaks(x, w = 1; minprom = nothing, maxprom = nothing, strict = true, N = nothing)

Find local maxima of `x` that remain maxima over a window of `w` points either side, returning
`(values, prominences, widths)`. Each is indexed by the position of its peak along the lookup of
`x`; `widths` holds the interval spanned by each peak at half its prominence.

`minprom`/`maxprom` bound the prominences kept, and may be given as a function of `x` (e.g.
`minprom = x -> 3std(x)`). `N` keeps only the `N` most prominent peaks. For an array of more than
one dimension, pass `dims` to select the axis to search; the other dimensions are carried through.

Wraps [Peaks.jl](https://github.com/halleysfifthinc/Peaks.jl). See [`maskpeaks`](@ref) to label
the samples each peak occupies.
"""
function findpeaks(
        x::DimensionalData.AbstractDimVector, w = 1;
        minprom = nothing,
        maxprom = nothing,
        strict = true, N = nothing
    )
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
    _dims = DimensionalData.dims(x)[
        DimensionalData.dims(x) .!= [
            DimensionalData.dims(
                x,
                dims
            ),
        ],
    ]
    P = findpeaks.(eachslice(x; dims = _dims); kwargs...)
    return [getindex.(P, i) for i in 1:3] # vals, proms, widths
end

"""
    maskpeaks!(y, x, args...; kwargs...)

In-place [`maskpeaks`](@ref), writing the labels into `y`.
"""
function maskpeaks!(y, x::DimensionalData.AbstractDimVector, args...; kwargs...)
    vals, proms, widths = findpeaks(x, args...; kwargs...)
    y .= 0
    for (i, I) in enumerate(widths)
        y[I] .= i
    end
    return y
end
"""
    maskpeaks(x, args...; dims = 1, kwargs...)

Label every sample of `x` with the index of the peak it belongs to, returning an integer array of
the same shape: samples within the width of the `i`th peak take the value `i`, and samples in no
peak take `0`. Arguments are passed to [`findpeaks`](@ref), so peak selection is controlled the
same way. See [`maskpeaks!`](@ref) for the in-place form.
"""
function maskpeaks(x::DimensionalData.AbstractDimVector, args...; kwargs...)
    y = set(x, similar(x, Int))
    maskpeaks!(y, x, args...; kwargs...)
    return y
end

function maskpeaks(x::DimensionalData.AbstractDimArray, args...; dims = 1, kwargs...)
    @assert length(dims) == 1
    _dims = DimensionalData.dims(x)[
        DimensionalData.dims(x) .!= [
            DimensionalData.dims(
                x,
                dims
            ),
        ],
    ]
    y = similar(x, Int)
    maskpeaks!.(
        eachslice(y; dims = _dims), eachslice(x; dims = _dims), args...;
        kwargs...
    )
    return y
end


_default_lags(x::AbstractVector) = range(1, length(x) - Int(length(x) ÷ 2), step = 1)
_default_lags(x::AbstractMatrix) = range(1, size(x, 1) - Int(size(x, 1) ÷ 2), step = 1)

"""
    madev(x, lags = 1:(length(x) ÷ 2); p = 1)

Mean absolute deviation of `x` at each of `lags`: the `p`-norm of `x[k+1:n] - x[1:n-k]`, divided by
the number of pairs. For a `RegularTimeseries` the lags are given in time units and the result is
returned as a `Timeseries` indexed by lag.

Unlike the mean-squared displacement this is finite for heavy-tailed processes, since it does not
rely on a variance: for a Lévy process with `α < 2` the MSD diverges while `madev` does not. Lags
must be sorted; a lag at or beyond the length of `x` contributes `0`.
"""
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
    lags = round.(Int, _lags ./ samplingperiod(x))
    d = dims(x)[2:end]
    m = mapslices(x -> madev(x, lags; kwargs...), parent(x); dims = 1)
    return Timeseries(m, _lags, d...)
end
