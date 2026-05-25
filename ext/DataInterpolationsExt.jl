module DataInterpolationsExt
using DataInterpolations
import DataInterpolations: AbstractInterpolation
using DimensionalData
import DimensionalData: AbstractDimArray, Dimension
using TimeseriesTools
using Unitful
import TimeseriesTools: interpolate, resample, impute

function interpolate(
        x::DimensionalData.AbstractDimVector,
        interp::Type{<:AbstractInterpolation} = AkimaInterpolation,
        args...; extrapolation = ExtrapolationType.Extension, kwargs...
    )
    DimensionalData.isforward(x) ||
        throw(ArgumentError("`interpolate` only supports forward-ordered data"))
    return interp(
        parent(x) .|> ustrip, x |> lookup |> only .|> ustrip, args...; extrapolation,
        kwargs...
    )
end
function (itp::AbstractInterpolation)(x::DimensionalData.Dimension; kwargs...)
    return DimensionalData.dimconstructor((x,))(
        (x |> lookup .|> ustrip |> itp),
        (x,);
        kwargs...
    )
end

_matches_missing(v, spec::Type) = v isa spec
_matches_missing(v, spec) = ismissing(spec) ? ismissing(v) :
    (v isa Number && spec isa Number && isnan(spec)) ? (v isa Number && isnan(v)) :
    isequal(v, spec)

function _impute_vec(
        x::DimensionalData.AbstractDimVector, interp, args...; replace, kwargs...
    )
    masked = map(parent(x)) do v
        any(spec -> _matches_missing(v, spec), replace) ? missing : v
    end
    kept = skipmissing(masked)
    u = isempty(kept) ? NoUnits : unit(first(kept))
    stripped = map(v -> ismissing(v) ? missing : ustrip(v), masked)
    itp = interpolate(rebuild(x; data = stripped), interp, args...; kwargs...)
    return rebuild(x; data = itp.(ustrip.(lookup(x, 1))) * u)
end

"""
    impute(x, interp = AkimaInterpolation, args...; dims = 1, replace = [NaN, Nothing, Missing], kwargs...)

Fill the flagged entries of a time series by interpolation.

Entries of `x` matching any element of `replace` are set to `missing`. An interpolant
of type `interp` is fit to the remaining samples (`DataInterpolations` drops the
`missing` pairs) and evaluated at every original time point. The returned series has
the same times as `x`. For arrays of more than one dimension, each slice along `dims`
is imputed independently.

`replace` accepts sentinel *values* (matched by `isequal`, with `NaN` matched by
`isnan`) and *types* (matched by `isa`). It defaults to `[NaN, Nothing, Missing]`; e.g.
`replace = [NaN, Nothing, Missing, Complex]` additionally treats every complex-valued entry
as missing. `args...`/`kwargs...` are forwarded to `interpolate`.
"""
function impute(
        x::DimensionalData.AbstractDimArray,
        interp::Type{<:AbstractInterpolation} = AkimaInterpolation,
        args...; dims = 1, replace = [NaN, Nothing, Missing], kwargs...
    )
    return TimeseriesTools._mapslices_dims(x; dims) do s
        _impute_vec(s, interp, args...; replace, kwargs...)
    end
end

"""
    resample(x, t_new, interp = AkimaInterpolation, args...; dims = 1, kwargs...)

Fit an interpolant to `x` along `dims` and evaluate it at `t_new`.

`x` may be regularly or irregularly sampled and uni- or multivariate. `t_new` can be:
- an `AbstractVector` or `DimensionalData.Dimension` of target sample points, or
- a `Number` (sampling period); the new grid runs from the first sample with step
  `t_new`, up to but not exceeding the last sample.

`interp` is any `DataInterpolations.AbstractInterpolation` subtype; `args...` and
`kwargs...` are forwarded to its constructor (e.g. `extrapolation`,
`cache_parameters`, BSpline degree, etc.).
"""
function resample(
        x::DimensionalData.AbstractDimArray, t_new, args...; dims = 1, kwargs...
    )
    d_old = DimensionalData.dims(x, dims)
    target = if t_new isa DimensionalData.Dimension
        t_new
    elseif t_new isa Number
        rebuild(d_old, TimeseriesTools._regular_grid(lookup(d_old), t_new))
    else
        rebuild(d_old, t_new)
    end

    u = unit(eltype(x))
    y = TimeseriesTools._mapslices_dims(x; dims) do s
        interpolate(s, args...; kwargs...)(target)
    end
    return y * u
end

end # module
