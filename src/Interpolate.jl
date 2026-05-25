export interpolate, upsample, downsample, resample

# * `interpolate`, `resample` and `downsample` are expanded in `ext/DataInterpolationsExt.jl`
# * and `ext/DSPExt.jl`.

"""
    interpolate(x, args...; kwargs...)

Fit an interpolant to a time series. Method-only function: the implementation is
provided by extensions, and the available signatures depend on which interpolation
backend is loaded.

- `using DataInterpolations` enables a per-axis 1-D method on `AbstractDimVector`,
  parameterised by a `DataInterpolations.AbstractInterpolation` type.
- `using DataInterpolationsND` enables a joint N-D method on `AbstractDimArray`,
  parameterised by `DataInterpolationsND.AbstractInterpolationDimension` types (or
  instances), one per axis.

Without either loaded, calling this function throws `MethodError`. See
[`resample`](@ref), [`upsample`](@ref), [`downsample`](@ref), and [`impute`](@ref) for
the higher-level operations built on top.
"""
function interpolate end

"""
    resample(x, target, args...; dims = 1, kwargs...)

Evaluate an interpolant of `x` on `target`, returning a new series with the same
non-interpolated dimensions. Method-only function: provided by extensions.

- `using DataInterpolations`: per-axis resampling along `dims = 1`. `target` may be an
  `AbstractVector` of new sample points, a `DimensionalData.Dimension`, or a `Number`
  (interpreted as a sampling period). Pass an
  `AbstractInterpolation` subtype as `args[1]` to select the interpolation method
  (default `AkimaInterpolation`).
- `using DataInterpolationsND`: joint multi-axis resampling. `target` is a tuple of new
  lookups (one per dimension of `x`) or a `Number` (common period along every axis), and
  `args[1]` is a `DataInterpolationsND.AbstractInterpolationDimension` type or tuple.

For *reducing* a regular sampling rate prefer [`downsample`](@ref) (filter-then-decimate);
plain `resample` onto a coarser grid does not anti-alias.
"""
function resample end

"""
    downsample(x::RegularTimeseries, factor::Integer; antialias = true)

Reduce the sampling rate of `x` by an integer `factor`. Method-only function: provided
by `DSPExt` (loaded with `using DSP`).

- `antialias = true` (default): apply an anti-aliasing FIR filter before decimating, so
  content above the new Nyquist is removed rather than aliased back into the kept band.
- `antialias = false`: plain `x[1:factor:end]`, no filter. Use to simulate having
  physically sampled the process at the lower rate.

See also [`upsample`](@ref), [`resample`](@ref).
"""
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
