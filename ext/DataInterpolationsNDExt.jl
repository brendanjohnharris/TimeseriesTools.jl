module DataInterpolationsNDExt
using DataInterpolationsND
import DataInterpolationsND: NDInterpolation, AbstractInterpolationDimension,
    LinearInterpolationDimension, BSplineInterpolationDimension,
    ConstantInterpolationDimension
using DimensionalData
using TimeseriesTools
using Unitful
import TimeseriesTools: interpolate, upsample, resample

# A type or instance describing how to interpolate one axis. Methods here dispatch on
# this so they never collide with the 1-D `DataInterpolations` (`AbstractInterpolation`)
# methods.
const NDDimSpec = Union{
    AbstractInterpolationDimension,
    Type{<:AbstractInterpolationDimension},
}

# Build an interpolation-dimension object for a single axis from its (unitless) lookup.
_nd_dim(::Type{LinearInterpolationDimension}, t) = LinearInterpolationDimension(collect(t))
function _nd_dim(::Type{ConstantInterpolationDimension}, t)
    return ConstantInterpolationDimension(collect(t))
end
function _nd_dim(::Type{BSplineInterpolationDimension}, t; degree = 2)
    # The data array supplies B-spline *control points*, one per lookup entry. With a default
    # clamped knot vector the basis-function count is `nknots + degree - 1`, so to match the
    # `length(t)` control points we use a reduced knot vector of `length(t) - degree + 1`
    # evenly spaced knots over the lookup span. For degree 1 this is the lookup itself (true
    # piecewise-linear interpolation); for higher degree the curve smooths the data rather
    # than passing through it (see the note in `interpolate`).
    nknots = length(t) - degree + 1
    nknots ≥ degree + 1 ||
        throw(ArgumentError("BSpline of degree $degree needs at least $(2degree) samples along a dimension"))
    knots = collect(range(first(t), last(t); length = nknots))
    return BSplineInterpolationDimension(knots, degree)
end
# Per-dimension spec: a single type/instance is broadcast across every axis; a tuple gives
# one spec per axis.
_nd_specs(spec, n) = ntuple(_ -> spec, n)
_nd_specs(spec::Tuple, n) = (length(spec) == n || throw(ArgumentError("expected $n dimension specs, got $(length(spec))")); spec)

_build_dim(spec::Type{<:AbstractInterpolationDimension}, t; kwargs...) = _nd_dim(spec, t; kwargs...)
_build_dim(spec::AbstractInterpolationDimension, _) = spec  # already constructed

"""
    interpolate(x, dimspec, args...; kwargs...)

Fit a joint N-dimensional interpolant to `x` using `DataInterpolationsND`.

`dimspec` is a `DataInterpolationsND.AbstractInterpolationDimension` type (or instance),
or a tuple of one per dimension of `x` — e.g. `LinearInterpolationDimension` to use linear
interpolation on every axis, or `(LinearInterpolationDimension, ConstantInterpolationDimension)`
for a 2-D array. Unlike the 1-D `DataInterpolations` method this fits all axes jointly, so
`x` must be a complete (gap-free) array.

`LinearInterpolationDimension` and `ConstantInterpolationDimension` pass through the samples
exactly. `BSplineInterpolationDimension` is **not** a true interpolation for `degree > 1`:
`DataInterpolationsND` treats the data as B-spline control points, so the values along that
axis are used as control points over a reduced knot vector — the result smooths the data
rather than passing through it (roughly, a smoothing followed by spline evaluation). For
`degree == 1` it reduces to ordinary piecewise-linear interpolation. The `degree` keyword
(default 2) is forwarded per BSpline axis.

Returns an `NDInterpolation`, which is callable on a tuple of target `DimensionalData`
dimensions to produce a rewrapped array on the new grid.
"""
function interpolate(
        x::DimensionalData.AbstractDimArray,
        dimspec::Union{NDDimSpec, Tuple},
        args...; kwargs...
    )
    specs = _nd_specs(dimspec, ndims(x))
    idims = ntuple(ndims(x)) do i
        _build_dim(specs[i], ustrip.(lookup(x, i)), args...; kwargs...)
    end
    return NDInterpolation(ustrip.(parent(x)), idims)
end

# Evaluate an NDInterpolation on a Cartesian grid given by target dimensions, returning a
# ToolsArray with those dimensions (units reattached by the caller).
function (itp::NDInterpolation)(targets::Tuple{Vararg{DimensionalData.Dimension}})
    grids = map(d -> ustrip.(lookup(d)), targets)
    out = Array{eltype(itp.u)}(undef, length.(grids)...)
    for I in CartesianIndices(out)
        pt = ntuple(k -> grids[k][I[k]], length(grids))
        out[I] = itp(pt)
    end
    return ToolsArray(out, targets)
end

"""
    upsample(x, factor, dimspec; dims = 1:ndims(x), kwargs...)

Increase the sampling density of `x` by `factor` along `dims`, using a joint
`DataInterpolationsND` fit (`dimspec` selects the per-axis interpolation method, as in
[`interpolate`](@ref)). Unlike the 1-D `upsample`, all `dims` are interpolated jointly in a
single pass rather than axis-by-axis. `x` must be gap-free.
"""
function upsample(
        x::DimensionalData.AbstractDimArray, factor::Number,
        dimspec::Union{NDDimSpec, Tuple}, args...;
        dims = Tuple(1:ndims(x)), kwargs...
    )
    u = unit(eltype(x))
    itp = interpolate(x, dimspec, args...; kwargs...)
    dimnums = dimnum(x, dims)
    targets = ntuple(ndims(x)) do i
        d = DimensionalData.dims(x, i)
        i in dimnums ? upsample(d, factor) : d
    end
    return (itp(targets)) * u
end

"""
    resample(x, targets, dimspec; kwargs...)

Resample `x` onto a new grid using a joint `DataInterpolationsND` fit. `targets` is a tuple
of new lookups/dimensions (one per dimension of `x`; pass the existing lookup to leave an
axis unchanged), or a single `Number` giving a common sampling period along every axis (the
grid runs from each axis' first sample up by that step, not exceeding the last). `dimspec`
selects the per-axis interpolation method, as in [`interpolate`](@ref). `x` must be gap-free.
"""
function resample(
        x::DimensionalData.AbstractDimArray, targets, dimspec::Union{NDDimSpec, Tuple},
        args...; kwargs...
    )
    u = unit(eltype(x))
    itp = interpolate(x, dimspec, args...; kwargs...)
    tdims = _nd_targets(x, targets)
    return (itp(tdims)) * u
end

# Build the tuple of target dimensions from a user `targets` argument.
function _nd_targets(x, targets::Number)
    return ntuple(ndims(x)) do i
        d = DimensionalData.dims(x, i)
        rebuild(d, TimeseriesTools._regular_grid(lookup(d), targets))
    end
end
function _nd_targets(x, targets::Tuple)
    length(targets) == ndims(x) ||
        throw(ArgumentError("expected $(ndims(x)) target lookups, got $(length(targets))"))
    return ntuple(ndims(x)) do i
        t = targets[i]
        t isa DimensionalData.Dimension ? t : rebuild(DimensionalData.dims(x, i), t)
    end
end

end # module
