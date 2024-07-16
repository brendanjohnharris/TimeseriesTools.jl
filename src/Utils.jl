import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import DimensionalData: print_array, _print_array_ctx, _print_matrix
import Normalization: NormUnion, AbstractNormalization, nansafe
using Peaks

export times, samplingrate, duration, samplingperiod, UnitPower, dimname, dimnames,
       describedim, describedims, describename, interlace, _buffer, buffer, window,
       delayembed, circularmean, circularstd, circularvar, resultant,
       resultantlength,
       centraldiff!, centraldiff, centralderiv!, centralderiv,
       rightdiff!, rightdiff, rightderiv!, rightderiv,
       rectify, phasegrad, addrefdim, addmetadata,
       findpeaks, maskpeaks, align, upsample, matchdim, nansafe, coarsegrain, meansquaredisp

import LinearAlgebra.mul!
function mul!(a::AbstractVector, b::AbstractTimeSeries, args...; kwargs...)
    mul!(a, b.data, args...; kwargs...)
end

Selectors = [:At, :Between, :Touches, :Near, :Where, :Contains]
# Allow dims to be passed directly to selectors
[:($(S)(D::Dimension) = $(S)(D.val.data)) for S in Selectors] .|> eval

description(x) = "$(size(x)) $(typeof(x).name.name)"
function print_array(io::IO, mime, A::AbstractDimArray{T, 0}) where {T <: AbstractArray}
    print(_print_array_ctx(io, T), "\n", description.(A[]))
end
function print_array(io::IO, mime, A::AbstractToolsArray{T, 1}) where {T <: AbstractArray}
    Base.print_matrix(_print_array_ctx(io, T), description.(A))
end
function print_array(io::IO, mime, A::AbstractToolsArray{T, 2}) where {T <: AbstractArray}
    Base.print_matrix(_print_array_ctx(io, T), description.(A))
end
function print_array(io::IO, mime, A::AbstractToolsArray{T, 3}) where {T <: AbstractArray}
    i3 = firstindex(A, 3)
    frame = view(parent(A), :, :, i3)
    println(io, "[:, :, $i3]")
    _print_matrix(_print_array_ctx(io, T), description.(frame), lookup(A, (1, 2)))
    nremaining = size(A, 3) - 1
    nremaining > 0 &&
        printstyled(io, "\n[and $nremaining more slices...]"; color = :light_black)
end
function print_array(io::IO, mime,
                     A::AbstractToolsArray{T, N}) where {T <: AbstractArray, N}
    o = ntuple(x -> firstindex(A, x + 2), N - 2)
    frame = view(A, :, :, o...)
    onestring = join(o, ", ")
    println(io, "[:, :, $(onestring)]")
    _print_matrix(_print_array_ctx(io, T), description.(frame), lookup(A, (1, 2)))
    nremaining = prod(size(A, d) for d in 3:N) - 1
    nremaining > 0 &&
        printstyled(io, "\n[and $nremaining more slices...]"; color = :light_black)
end
function print_array(io::IO, mime, A::SpikeTrain{Bool, 1})
    _print_array_ctx(io, Bool)
end

function Base.stack(D::DimensionalData.Dimension,
                    args::AbstractVector{<:AbstractToolsArray};
                    dims = nothing, kwargs...)
    x = first(args) # Is this allocating?

    isnothing(dims) && (dims = ndims(x) + 1)
    if !all([size(x)] .== size.(args))
        error("Input arrays must have the same dimensionality and size")
    end
    unidims = DimensionalData.dims(first(args))
    if !all(DimensionalData.dims.(args) .== [unidims])
        error("Input arrays must have the same dimensions")
    end
    if dims isa Val && typeof(dims).parameters[1] isa Integer
        dims = typeof(dims).parameters[1]
    end
    if !(dims isa Integer) && !(dims isa Val)
        idx = DimensionalData.dims(x) isa dims ? 1 :
              findfirst(isa.(DimensionalData.dims(x), [dims]))
        isnothing(idx) && error("Dimension $dims not found in input array")
        dims = first(idx)
    end
    outsize = size(x) |> collect
    insert!(outsize, dims, length(D))

    X = Array{eltype(x), length(outsize)}(undef, Tuple(outsize))

    for (i, _x) in enumerate(eachslice(X; dims, drop = false))
        _x .= args[i]
    end

    ds = Vector{Any}([DimensionalData.dims(x)...])
    insert!(ds, dims, D)
    y = ToolsArray(X, (ds...,); refdims = refdims(x), name = name(x),
                   metadata = metadata(x))
    return y
end

"""
    Base.cat(D::DimensionalData.Dimension, args...; kwargs...)
Concatenate the arrays given in `args...`, and give the resulting extra axis dimensions `D`.
Note that unlike `Base.cat` without the first `Dim` argument, this increments all existing dimensions greater than `dims` by one (so N n×n arrays concatenated at `dims=1` will result in an N×n×n array).
`args...` can be a splatted collection of `ToolsArray`s, but this will give the same behaviour as if `args...` is a single vector of `ToolsArray`s; the latter is much more performant.
"""
function Base.cat(D::DimensionalData.Dimension, x::AbstractToolsArray,
                  y::AbstractToolsArray,
                  args...;
                  dims = nothing,
                  kwargs...) # TODO Refactor this method, to call the method above.
    isnothing(dims) && (dims = ndims(x) + 1)
    if !all([size(x)] .== size.(args))
        error("Input arrays must have the same dimensionality and size")
    end
    if dims isa Val && typeof(dims).parameters[1] isa Integer
        dims = typeof(dims).parameters[1]
    end
    if !(dims isa Integer) && !(dims isa Val)
        idx = DimensionalData.dims(x) isa dims ? 1 :
              findfirst(isa.(DimensionalData.dims(x), [dims]))
        isnothing(idx) && error("Dimension $dims not found in input array")
        dims = first(idx)
    end
    function rf(x)
        r = size(x) |> collect
        insert!(r, dims, 1)
        reshape(x, r...)
    end
    _x = rf(x.data)
    _y = rf(y.data)
    _args = rf.(getfield.(args, [:data]))
    x′ = cat(_x, _y, _args...; dims, kwargs...)
    ds = Vector{Any}([DimensionalData.dims(x)...])
    insert!(ds, dims, D)
    y = ToolsArray(x′, (ds...,); refdims = refdims(x), name = name(x),
                   metadata = metadata(x))
    # if hasdim(y, Ti)
    #     ts = times(y)
    #     y = set(y, Ti => ts .- minimum(ts))
    # end
    return y
end

"""
    times(x::AbstractTimeSeries)

Returns the time indices of the [`AbstractTimeSeries`](@ref) `x`.

## Examples
```@example 1
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x);
julia> times(ts) == t
```
"""
times(x::AbstractTimeSeries) = dims(x, Ti).val.data

"""
    step(x::RegularTimeSeries)

Returns the step size (time increment) of a regularly sampled [`RegularTimeSeries`](@ref).

## Examples
```@example 1
julia> t = 1:100;
julia> x = rand(100);
julia> rts = TimeSeries(t, x);
julia> step(rts) == 1
```
"""
Base.step(x::RegularTimeSeries; dims = Ti) = DimensionalData.dims(x, dims).val.data |> step

"""
    samplingrate(x::RegularTimeSeries)

Returns the sampling rate (inverse of the step size) of a regularly sampled [`RegularTimeSeries`](@ref).

## Examples
```@example 1
julia> t = 1:100;
julia> x = rand(100);
julia> rts = TimeSeries(t, x);
julia> samplingrate(rts) == 1
```
"""
samplingrate(x::RegularTimeSeries; kwargs...) = 1 / step(x; kwargs...)

"""
    samplingperiod(x::RegularTimeSeries)

Returns the sampling period (step size) of a regularly sampled [`RegularTimeSeries`](@ref).

## Examples
```@example 1
julia> t = 1:100;
julia> x = rand(100);
julia> rts = TimeSeries(t, x);
julia> samplingperiod(rts) == 1
```
"""
samplingperiod(x::RegularTimeSeries; kwargs...) = step(x; kwargs...)

"""
    duration(x::AbstractTimeSeries)

Returns the duration of the [`AbstractTimeSeries`](@ref) `x`.

## Examples
```@example 1
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x);
julia> TimeseriesTools.duration(ts) == 99
```
"""
duration(x::AbstractTimeSeries) = (last ∘ times)(x) - (first ∘ times)(x)

"""
    IntervalSets.Interval(x::AbstractTimeSeries)

Returns an interval representing the range of the [`AbstractTimeSeries`](@ref) `x`.

## Examples
```@example 1
julia> using IntervalSets;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x);
julia> IntervalSets.Interval(ts) == (1..100)
```
"""
IntervalSets.Interval(x::AbstractTimeSeries) = (first ∘ times)(x) .. (last ∘ times)(x)

# function 𝑝(x::RegularTimeSeries)
#     dur = duration(x)
#     if ~isnothing(unit(dur))
#         return sum(x.^2)/dur
#     else
#         @warn "No time units found for unit power normalization. Assuming seconds."
#         return sum(x.^2)/(dur*u"s")
#     end
# end
𝑝(x::RegularTimeSeries) = sum(x .^ 2) / duration(x) # * Assume seconds.
# TODO !!!
"""
    UnitPower <: AbstractNormalization

A normalization that sets the total power of a signal to unity.

# Fields
- `dims`: The dimensions to normalize over.
- `p`: Computed normalization parameters.
- `𝑝`: A function that returns the power from a given time series.
- `𝑓`: The normalization method
- `𝑓⁻¹`: The inverse normalization method.

"""
mutable struct UnitPower{T} <: AbstractNormalization{T}
    dims::Any
    p::NTuple{1, AbstractArray{T}}
    𝑝::NTuple{1, Function}
    𝑓::Function
    𝑓⁻¹::Function
end;

function UnitPower{T}(; dims = nothing,
                      p = (Vector{T}(),),
                      𝑝 = (𝑝,),
                      𝑓 = (x, 𝑃) -> x .= x ./ sqrt.(𝑃),
                      𝑓⁻¹ = (y, 𝑃) -> y .= y .* sqrt.(𝑃)) where {T}
    UnitPower(((isnothing(dims) || length(dims) < 2) ? dims : sort(dims)), p, 𝑝, 𝑓, 𝑓⁻¹)
end

UnitPower(; kwargs...) = UnitPower{Nothing}(; kwargs...);

dimname(d::DimensionalData.Dimension) = name(d) |> string
dimname(x::AbstractDimArray, dim) = dims(x, dim) |> dimname
dimname(x::AbstractDimArray) = map(dimname, dims(x))
dimnames = dimname

function describedim(d::DimensionalData.Dimension)
    if d isa DimensionalData.TimeDim
        n = "Time"
    elseif d isa FrequencyDim
        n = "Frequency"
    elseif d isa VariableDim
        n = "Variable"
    else
        n = name(d)
    end
    units = d |> eltype |> unit
    isnothing(units) && (n = n * " ($units)")
    n
end

function describedim(x::Tuple)
    @assert eltype(x) <: DimensionalData.Dimension # Move this into function args
    ns = describedim.(x)
    return ns
end
describedim(x::AbstractDimArray, i) = dims(x, i) |> describedim
describedim(x::AbstractDimArray) = x |> dims |> describedim
describedims = describedim

function describename(x::AbstractDimArray)
    n = name(x)
    n isa DimensionalData.NoName && (n = "")
    if ~isnothing(n)
        units = x |> eltype |> unit
        isnothing(units) && (n = n * " ($units)")
    end
    n
end

function interlace(x::AbstractTimeSeries, y::AbstractTimeSeries)
    ts = vcat(times(x), times(y))
    idxs = sortperm(ts)
    ts = ts[idxs]
    data = vcat(x.data, y.data)
    data = data[idxs]
    return TimeSeries(ts, data)
end

function _buffer(x, n::Integer, p::Integer = 0; discard::Bool = true)
    y = [@views x[i:min(i + n - 1, end)] for i in 1:(n - p):length(x)]
    while discard && length(y[end]) < n
        pop!(y)
    end
    y
end
function _buffer(x::AbstractMatrix, n::Integer, p::Integer = 0; discard::Bool = true)
    y = [@views x[i:min(i + n - 1, end), :] for i in 1:(n - p):size(x, 1)]
    while discard && size(y[end], 1) < n
        pop!(y)
    end
    y
end
buffer(x::AbstractVector, args...; kwargs...) = _buffer(x, args...; kwargs...)

"""
    buffer(x::RegularTimeSeries, n::Integer, p::Integer; kwargs...)

Buffer a time series `x` with a given window length and overlap between successive buffers.

## Arguments
- `x`: The regular time series to be buffered.
- `n`: The number of samples in each buffer.
- `p`: The number of samples of overlap betweeen the buffers.
    - `0` indicates no overlap
    - +`2` indicates `2` samples of overlap between successive buffers
    - -`2` indicates `2` samples of gap between buffers

See also: [`window`](@ref), [`delayembed`](@ref), [`coarsegrain`](@ref)
"""
function buffer(x::RegularTimeSeries, args...; kwargs...)
    y = _buffer(x, args...; kwargs...)
    t = _buffer(times(x), args...; kwargs...) .|> mean
    # For a regular time series, the buffer centres are regular
    ts = range(first(t), last(t), length(y))
    y = TimeSeries(ts, y)
end

"""
    window(x::RegularTimeSeries, n::Integer, p::Integer; kwargs...)

Window a time series `x` with a given window length and step between successive windows.

## Arguments
- `x`: The regular time series to be windows.
- `n`: The number of samples in each window.
- `p`: The number of samples to slide each successive window.

See also: [`buffer`](@ref), [`delayembed`](@ref), [`coarsegrain`](@ref)
"""
window(x, n, p = n, args...; kwargs...) = buffer(x, n, n - p, args...; kwargs...)

function _delayembed(x::AbstractVector, n, τ, p = 1; kwargs...) # A delay embedding with dimension `n`, delay `τ`, and skip length of `p`
    y = window(x, n * τ, p; kwargs...)
    y = map(y) do _y
        @view _y[1:τ:end]
    end
end
delayembed(x::AbstractVector, args...; kwargs...) = _delayembed(x, args...; kwargs...)

"""
    delayembed(x::UnivariateRegular, n::Integer, τ::Integer, p::Integer=1; kwargs...)

Delay embed a univariate time series `x` with a given dimension `n`, delay `τ`, and skip length of `p`

## Arguments
- `x`: The regular time series to be delay embedded.
- `n`: The embedding dimension, i.e., the number of samples in each embedded vector.
- `τ`: The number of original sampling periods between each sample in the embedded vectors.
- `p`: The number of samples to skip between each successive embedded vector.

See also: [`buffer`](@ref), [`window`](@ref)
"""
function delayembed(x::UnivariateRegular, n, τ, p = 1, args...; kwargs...)
    y = _delayembed(x, n, τ, p, args...; kwargs...)
    ts = last.(times.(y))  # Time of the head of the vector
    dt = step(x) * p
    ts = ts[1]:dt:(ts[1] + dt * (length(y) - 1))
    δt = τ * p * step(x)
    delays = (-(δt * (n - 1))):δt:0
    y = set.(y, [Ti => Dim{:delay}(delays)])
    y = cat(Ti(ts), y..., dims = Dim{:delay})
end

function rectify(ts::DimensionalData.Dimension; tol = 4, zero = false, extend = false,
                 atol = nothing)
    u = unit(eltype(ts))
    ts = collect(ts)
    origts = ts
    stp = ts |> diff |> mean
    err = ts |> diff |> std
    tol = Int(tol - round(log10(stp |> ustripall)))

    if isnothing(atol) && ustripall(err) > exp10(-tol - 1)
        @warn "Step $stp is not approximately constant (err=$err, tol=$(exp10(-tol-1))), skipping rectification"
    else
        if !isnothing(atol)
            tol = atol
        end
        stp = u == NoUnits ? round(stp; digits = tol) : round(u, stp; digits = tol)
        t0, t1 = u == NoUnits ? round.(extrema(ts); digits = tol) :
                 round.(u, extrema(ts); digits = tol)
        if zero
            origts = t0:stp:(t1 + (10000 * stp))
            t1 = t1 - t0
            t0 = 0
        end
        if extend
            ts = t0:stp:(t1 + (10000 * stp))
        else
            ts = range(start = t0, step = stp, length = length(ts))
        end
    end
    return ts, origts
end
rectifytime(ts::Ti; kwargs...) = rectify(ts; kwargs...)

function rectify(X::AbstractDimArray; dims, tol = 4, zero = false, kwargs...) # tol gives significant figures for rounding
    if !(dims isa Tuple || dims isa AbstractVector)
        dims = [dims]
    end
    for dim in dims
        ts, origts = rectify(DimensionalData.dims(X, dim); tol, zero, extend = true,
                             kwargs...)
        ts = ts[1:size(X, dim)] # Should be ok?
        @assert length(ts) == size(X, dim)
        X = set(X, dim => ts)
        if zero
            X = rebuild(X; metadata = (Symbol(dim) => origts, pairs(metadata(X))...))
        end
    end
    return X
end

"""
    rectifytime(X::IrregularTimeSeries; tol = 6, zero = false)

Rectifies the time values of an [`IrregularTimeSeries`](@ref). Checks if the time step of
the input time series is approximately constant. If it is, the function rounds the time step
and constructs a [`RegularTimeSeries`](@ref) with range time indices. If the time step is
not approximately constant, a warning is issued and the rectification is skipped.

# Arguments
- `X::IrregularTimeSeries`: The input time series.
- `tol::Int`: The number of significant figures for rounding the time step. Default is 6.
- `zero::Bool`: If `true`, the rectified time values will start from zero. Default is
  `false`.
"""
rectifytime(X::IrregularTimeSeries; kwargs...) = rectify(X; dims = Ti, kwargs...)

function rectifytime(X::AbstractVector; tol = 6, zero = false) # ! Legacy
    # Generate some common time indices as close as possible to the rectified times of each element of the input vector
    ts = times.(X)
    mint = maximum(minimum.(ts)) - exp10(-tol) .. minimum(maximum.(ts)) + exp10(-tol)
    X = [x[Ti(mint)] for x in X]
    X = [x[1:minimum(size.(X, 1))] for x in X]
    ts = mean(times.(X))
    ts, origts = rectifytime(Ti(ts); tol, zero)
    ts = ts[1:size(X[1], Ti)] # Should be ok?
    if any([any(ts .- times(x) .> std(ts) / 10.0^(-tol)) for x in X])
        @error "Cannot find common times within tolerance"
    end

    X = [set(x, Ti => ts) for x in X]
    return X
end

function matchdim(X::AbstractVector{<:AbstractDimArray}; dims = 1, tol = 4, zero = false,
                  kwargs...)
    # Generate some common time indices as close as possible to the rectified times of each element of the input vector. At most this will change each time index by a maximum of 1 sampling period. We could do better--maximum of a half-- but leave that for now.
    u = lookup(X |> first, dims) |> eltype |> unit
    ts = lookup.(X, [dims])
    mint = (maximum(minimum.(ts)) - exp10(-tol) * u) ..
           (minimum(maximum.(ts)) + exp10(-tol) * u)
    X = map(X) do x
        d = rebuild(DimensionalData.dims(x, dims), mint)
        x = getindex(x, d)
    end
    L = minimum(size.(X, dims))
    X = map(X) do x
        d = rebuild(DimensionalData.dims(x, dims), 1:L)
        x = getindex(x, d) # Should now have same length for all inputs
    end

    ts = mean(lookup.(X, [dims]))
    ts, origts = rectify(rebuild(DimensionalData.dims(X[1], dims), ts); tol, zero,
                         kwargs...)
    if any([any(ts .- lookup(x, dims) .> std(ts) / exp10(-tol)) for x in X])
        @error "Cannot find common dimension indices within tolerance"
    end
    X = [set(x, rebuild(DimensionalData.dims(x, dims), ts)) for x in X]
    return X
end

phasegrad(x::Real, y::Real) = mod(x - y + π, 2π) - π # +pi - pi because we want the difference mapped from -pi to +pi, so we can represent negative changes.
phasegrad(x, y) = phasegrad.(x, y)
phasegrad(x::Complex, y::Complex) = phasegrad(angle(x), angle(y))

function _centraldiff!(x; grad = -, dims = nothing) # Dims unused
    # a = x[2] # Save here, otherwise they get mutated before we use them
    # b = x[end - 1]
    if grad == -
        x[2:(end - 1)] .= grad(x[3:end], x[1:(end - 2)]) / 2
    else # For a non-euclidean metric, we need to calculate both sides individually
        x[2:(end - 1)] .= (grad(x[3:end], x[2:(end - 1)]) +
                           grad(x[2:(end - 1)], x[1:(end - 2)])) / 2
    end
    # x[[1, end]] .= [grad(a, x[1]), grad(x[end], b)]
    x[[1, end]] .= [copy(x[2]), copy(x[end - 1])]
    return nothing
end

_diff!(x::UnivariateRegular, f!; kwargs...) = f!(x; kwargs...)
function _diff!(x::AbstractDimArray, f!; dims = 1, kwargs...)
    if !(DimensionalData.lookup(x, dims).data isa AbstractRange)
        error("Differencing dimension must be regularly sampled")
    end
    f!(eachslice(x; dims); kwargs...)
end

"""
    centraldiff!(x::RegularTimeSeries; dims=Ti, grad=-)

Compute the central difference of a regular time series `x`, in-place.
The first and last elements are set to the forward and backward difference, respectively.
The dimension to perform differencing over can be specified as `dims`, and the differencing function can be specified as `grad` (defaulting to the euclidean distance, `-`)
"""
centraldiff!(args...; kwargs...) = _diff!(args..., _centraldiff!; kwargs...)

function _diff(x::RegularTimeSeries, f!; kwargs...)
    y = deepcopy(x)
    f!(y; kwargs...)
    return y
end
"""
    centraldiff(x::RegularTimeSeries; dims=Ti, grad=-)

Compute the central difference of a regular time series `x`.
The first and last elements are set to the forward and backward difference, respectively.
The dimension to perform differencing over can be specified as `dims`, and the differencing function can be specified as `grad` (defaulting to the euclidean distance, `-`)
See [`centraldiff!`](@ref).
"""
centraldiff(args...; kwargs...) = _diff(args..., centraldiff!; kwargs...)

function checkderivdims(dims)
    if dims isa Tuple || dims isa AbstractVector
        error("Only one dimension can be specified for derivatives.")
    end
end

function _deriv!(x::RegularTimeSeries, f!; dims = Ti, kwargs...)
    checkderivdims(dims)
    f!(x; dims, kwargs...)
    x ./= step(x; dims)
    nothing
end

"""
    centralderiv!(x::RegularTimeSeries; kwargs...)

Compute the central derivative of a regular time series `x`, in-place.
See [`centraldiff!`](@ref) for available keyword arguments.
"""
centralderiv!(args...; kwargs...) = _deriv!(args..., centraldiff!; kwargs...)

function _deriv(x::RegularTimeSeries, f!; dims = Ti, kwargs...)
    y = deepcopy(x)
    if unit(step(x; dims)) == NoUnits # Can safely mutate
        f!(y; kwargs...)
    else
        y = ustripall(y)
        f!(y; dims, kwargs...)
        newu = unit(eltype(x)) / unit(step(x; dims))
        y = set(x, y .* newu)
    end
    return y
end
"""
    centralderiv(x::AbstractTimeSeries)

Compute the central derivative of a time series `x`.
See [`centraldiff`](@ref) for available keyword arguments.
Also c.f. [`centralderiv!`](@ref).
"""
centralderiv(args...; kwargs...) = _deriv(args..., centralderiv!; kwargs...)

function _rightdiff!(x; grad = -, dims = nothing) # Dims unused
    x[1:(end - 1)] .= grad(x[2:end], x[1:(end - 1)])
    # x[[1, end]] .= [grad(a, x[1]), grad(x[end], b)]
    x[[end]] .= [copy(x[end - 1])]
    return nothing
end
rightdiff!(args...; kwargs...) = _diff!(args..., _rightdiff!; kwargs...)
rightdiff(args...; kwargs...) = _diff(args..., rightdiff!; kwargs...)
rightderiv!(args...; kwargs...) = _deriv!(args..., rightdiff!; kwargs...)
rightderiv(args...; kwargs...) = _deriv(args..., rightderiv!; kwargs...)

function _leftdiff!(x; grad = -, dims = nothing) # Dims unused
    x[2:end] .= grad(x[2:end], x[1:(end - 1)])
    # x[[1, end]] .= [grad(a, x[1]), grad(x[end], b)]
    x[[1]] .= [copy(x[2])]
    return nothing
end
leftdiff!(args...; kwargs...) = _diff!(args..., _leftdiff!; kwargs...)
leftdiff(args...; kwargs...) = _diff(args..., leftdiff!; kwargs...)
leftderiv!(args...; kwargs...) = _deriv!(args..., leftdiff!; kwargs...)
leftderiv(args...; kwargs...) = _deriv(args..., leftderiv!; kwargs...)

Base.abs(x::AbstractTimeSeries) = Base.abs.(x)
Base.angle(x::AbstractTimeSeries) = Base.angle.(x)

# * See https://en.wikipedia.org/wiki/Directional_statistics
resultant(θ; kwargs...) = mean(exp.(im .* θ); kwargs...)
resultantlength(θ; kwargs...) = abs.(resultant(θ; kwargs...))
circularmean(θ; kwargs...) = angle.(resultant(θ; kwargs...))
circularvar(θ; kwargs...) = 1 - resultantlength(θ; kwargs...)
circularstd(θ; kwargs...) = sqrt.(-2 * log.(resultantlength(θ; kwargs...)))

## Add refdims to a DimArray
function addrefdim(X::AbstractDimArray, dim::DimensionalData.Dimension)
    rebuild(X; dims = dims(X),
            metadata = DimensionalData.metadata(X),
            name = DimensionalData.name(X),
            refdims = (DimensionalData.refdims(X)..., dim))
end

function addmetadata(X::AbstractDimArray; kwargs...)
    p = DimensionalData.metadata(X)
    p = p isa DimensionalData.Metadata ? pairs(p.val) : pairs(p)
    if any(keys(kwargs) .∈ [keys(p)])
        @warn "Metadata already contains one of the keys, overwriting $(collect(pairs(kwargs))[keys(kwargs) .∈ [keys(p)]])"
    end
    md = DimensionalData.Metadata(p..., kwargs...)
    rebuild(X; dims = dims(X),
            metadata = md,
            name = DimensionalData.name(X),
            refdims = DimensionalData.refdims(X))
end

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

"""
    align(x::AbstractDimArray, ts, dt; dims = 1)

Align a `DimArray` `x` to each of a set of dimension values `ts`, selecting a window given by `dt` centered at each element of `ts`.
`dt` can be a two-element vector/tuple, or an interval.
The `dims` argument specifies the dimension along which the alignment is performed.
Each element of the resulting `DimArray` is an aligned portion of the original `x`.
"""
function align(x::DimensionalData.AbstractDimArray, ts,
               dt::Union{<:Tuple, <:AbstractVector}; dims = 1, zero = true)
    @assert length(dims) == 1
    dims isa Integer &&
        (dims = DimensionalData.dims(x, dims))
    ints = [Interval((t .+ dt)...) for t in ts]
    x = TimeSeries(ts, [view(x, rebuild(dims, i)) for i in ints])
    if zero
        x = set(x, map(enumerate(x)) do (i, _x)
                    set(_x, dims => lookup(_x, dims) .- ts[i])
                end)
    end
    return x
end
align(x, ts, dt::Interval; kwargs...) = align(x, ts, extrema(dt); kwargs...)

function upsample(d::DimensionalData.Dimension{<:RegularIndex}, factor::Number)
    rebuild(d, range(start = minimum(d), stop = maximum(d), step = step(d) / factor))
end
function upsample(d::DimensionalData.Dimension, factor)
    rebuild(d,
            range(start = minimum(d), stop = maximum(d),
                  step = mean(diff(lookup(d))) / factor))
end

"""
    stitch(x, args...)

Stitch multiple time series together by concatenatign along the time dimension generating new contiguous time indices. The time series must be of the same type (`UnivariateRegular`, `MultivariateRegular`, or `AbstractArray`), and the sampling period and dimensions of the data arrays must match. If the arguments are `MultivariateRegular, they must have the same dimensions (except for the time dimension).

# Arguments
- `X`: The first time series.
- `args...`: Additional time series.

# Returns
- A new time series containing the concatenated data.
"""
function stitch(x::UnivariateRegular, y::UnivariateRegular)
    dt = samplingperiod(x)
    @assert dt == samplingperiod(y)
    z = vcat(x.data, y.data)
    z = TimeSeries(dt:dt:(dt * size(z, 1)), z)
end
stitch(x::AbstractArray, y::AbstractArray) = vcat(x, y)
function stitch(x::MultivariateRegular, y::MultivariateRegular)
    dt = samplingperiod(x)
    @assert dt == samplingperiod(y)
    @assert all(dims(x)[2:end] .== dims(y)[2:end])
    z = vcat(x.data, y.data)
    z = TimeSeries(dt:dt:(dt * size(z, 1)), dims(x)[2:end]..., z)
end
stitch(X, Y, args...) = reduce(stitch, (X, Y, args...))

"""
    coarsegrain(X::AbstractArray; dims = nothing, newdim=ndims(X)+1)
Coarse-grain an array by taking every second element over the given dimensions `dims` and concatenating them in the dimension `newdim`. `dims` are coarse-grained in sequence, from last to first. If `dims` is not specified, we iterate over all dimensions that are not `newdim`. If the array has an odd number of slices in any `dims`, the last slice is discarded.
This is more flexibile than the conventional, mean-based definition of coarse graining: it can be used to generate coarse-grained distributions from an array. To recover this conventional mean-based coarse-graining:
```julia
    C = coarsegrain(X)
    mean(C, dims=ndims(C))
```
"""
function coarsegrain(X::AbstractArray; dims = nothing, newdim = ndims(X) + 1)
    if isnothing(dims)
        dims = collect(1:ndims(X))
        dims = setdiff(dims, newdim)
    end
    dims = collect(Tuple(dims))
    if newdim ∈ dims
        error("`dims` cannot contain `newdim`")
    end
    all(size(X)[dims] .> 1) ||
        error("Cannot coarse-grain a dimension with only one element")
    while !isempty(dims)
        dim = pop!(dims)
        𝒳 = eachslice(X; dims = dim)
        N = floor(Int, length(𝒳) / 2)
        X = cat(stack(𝒳[1:2:(N * 2)], dims = dim), stack(𝒳[2:2:(N * 2)], dims = dim),
                dims = newdim)
    end
    return X
end

function coarsegrain(X::AbstractDimArray; dims = nothing,
                     newdim = ndims(X) + 1)
    if isnothing(dims)
        dims = DimensionalData.dims(X)
    end
    _dims = [dimnum(X, dims)...]
    dims = DimensionalData.dims.([X], _dims)
    if hasdim(X, newdim)
        _newdim = dimnum(X, newdim)
        newdim = DimensionalData.dims(X, _newdim)
    else
        _newdim = ndims(X) + 1
    end
    while !isempty(_dims)
        _dim = pop!(_dims)
        _X = coarsegrain(X.data; dims = _dim, newdim = _newdim)
        N = floor(Int, size(X, _dim) / 2)
        if hasdim(X, newdim)
            newdim = rebuild(DimensionalData.dims(X, newdim),
                             vcat(DimensionalData.dims(X, _newdim).val,
                                  DimensionalData.dims(X, _newdim).val))
            newdims = collect(Any, DimensionalData.dims(X))
            newdims[_dim] = rebuild(newdims[_dim],
                                    (newdims[_dim][1:2:(N * 2)] .+
                                     newdims[_dim][2:2:(N * 2)]) ./ 2)
            newdims[_newdim] = newdim
        else
            newdims = collect(Any, DimensionalData.dims(X))
            newdims[_dim] = rebuild(newdims[_dim],
                                    (newdims[_dim][1:2:(N * 2)] .+
                                     newdims[_dim][2:2:(N * 2)]) ./ 2)
            newdims = [newdims..., DimensionalData.AnonDim(1:size(_X, _newdim))]
            newdim = newdims[newdim]
        end
        X = ToolsArray(_X, Tuple(newdims); refdims = refdims(X), name = name(X),
                       metadata = metadata(X))
    end

    return X
end
