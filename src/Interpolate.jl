export interpolate, upsample, downsample, resample

# * `interpolate`, `resample` and `downsample` are expanded in `ext/DataInterpolationsExt.jl`
# * and `ext/DSPExt.jl`.

function interpolate end
function resample end
function downsample end


function upsample(d::DimensionalData.Dimension{<:RegularIndex}, factor::Number)
    return rebuild(d, range(start = minimum(d), stop = maximum(d), step = step(d) / factor))
end
function upsample(d::DimensionalData.Dimension, factor)
    return rebuild(
        d,
        range(
            start = minimum(d), stop = maximum(d),
            step = mean(diff(lookup(d))) / factor
        )
    )
end

"""
    upsample(x::AbstractDimArray, factor, args...; dims = 1, kwargs...)

Increase the sampling density of `x` by `factor` along `dims` by interpolation. A thin
wrapper over [`resample`](@ref): for each dimension in `dims` it builds a denser target
grid (`step / factor`, or `mean(diff) / factor` for irregular lookups) and resamples that
axis onto it. `args...`/`kwargs...` (e.g. the interpolation type) are forwarded to
`resample`.

To *reduce* the sampling rate it is suggested
to use [`downsample`](@ref), which anti-alias filters before decimating.
"""
function upsample(
        x::DimensionalData.AbstractDimArray, factor::Number, args...; dims = 1, kwargs...
    )
    y = x
    for d in Tuple(dims)
        target = upsample(DimensionalData.dims(y, d), factor)
        y = resample(y, target, args...; dims = d, kwargs...)
    end
    return y
end

# A regular grid spanning `t` with step `dt`, from `first(t)` up to (not exceeding) `last(t)`.
# `stop` is a soft bound that `range` snaps to the last whole step. Shared by the
# interpolation extensions for `Number`-period resampling; `dt` must share the lookup's units.
_regular_grid(t, dt) = range(start = first(t), step = dt, stop = last(t))

# Map `f` over the slices of `x` orthogonal to `dims`, restacking into a same-named array.
# Each slice is the 1-D view along `dims`; `f` returns the transformed slice. Used by the
# interpolation extensions to apply a per-axis operation across a multidimensional array.
function _mapslices_dims(f, x::DimensionalData.AbstractDimArray; dims = 1)
    negdims = setdiff(1:ndims(x), DimensionalData.dimnum(x, dims)) |> Tuple
    isempty(negdims) && return f(x)
    y = stack(map(f, eachslice(x; dims = negdims)))
    return permutedims(y, DimensionalData.dimnum.([y], DimensionalData.dims(x)))
end
