import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import DimensionalData: print_array, _print_array_ctx
import Normalization: NormUnion, AbstractNormalization

export times, samplingrate, duration, samplingperiod, UnitPower, dimname, dimnames,
       describedim, describedims, describename, interlace, _buffer, buffer, window,
       delayembed, circmean, circvar, circstd, circresultant, circlength,
       centraldiff!, centraldiff, centralderiv!, centralderiv,
       rightdiff!, rightdiff, rightderiv!, rightderiv,
       rectify, phasegrad, addrefdim, addmetadata

import LinearAlgebra.mul!
function mul!(a::AbstractVector, b::AbstractTimeSeries, args...; kwargs...)
    mul!(a, b.data, args...; kwargs...)
end

Selectors = [:At, :Between, :Touches, :Near, :Where, :Contains]
# Allow dims to be passed directly to selectors
[:($(S)(D::Dimension) = $(S)(D.val.data)) for S in Selectors] .|> eval

describate(x) = "$(size(x)) $(typeof(x).name.name)"
function print_array(io::IO, mime, A::AbstractDimArray{T, 0}) where {T <: AbstractArray}
    print(_print_array_ctx(io, T), "\n", describate.(A[]))
end
function print_array(io::IO, mime, A::AbstractDimArray{T, 1}) where {T <: AbstractArray}
    Base.print_matrix(_print_array_ctx(io, T), describate.(A))
end
function print_array(io::IO, mime, A::AbstractDimArray{T, 2}) where {T <: AbstractArray}
    Base.print_matrix(_print_array_ctx(io, T), describate.(A))
end
function print_array(io::IO, mime, A::AbstractDimArray{T, 3}) where {T <: AbstractArray}
    i3 = firstindex(A, 3)
    frame = view(parent(A), :, :, i3)
    println(io, "[:, :, $i3]")
    _print_matrix(_print_array_ctx(io, T), describate.(frame), lookup(A, (1, 2)))
    nremaining = size(A, 3) - 1
    nremaining > 0 &&
        printstyled(io, "\n[and $nremaining more slices...]"; color = :light_black)
end
function print_array(io::IO, mime, A::AbstractDimArray{T, N}) where {T <: AbstractArray, N}
    o = ntuple(x -> firstindex(A, x + 2), N - 2)
    frame = view(A, :, :, o...)
    onestring = join(o, ", ")
    println(io, "[:, :, $(onestring)]")
    _print_matrix(_print_array_ctx(io, T), describate.(frame), lookup(A, (1, 2)))
    nremaining = prod(size(A, d) for d in 3:N) - 1
    nremaining > 0 &&
        printstyled(io, "\n[and $nremaining more slices...]"; color = :light_black)
end
function print_array(io::IO, mime, A::SpikeTrain{Bool, 1})
    _print_array_ctx(io, Bool)
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
𝑝(x::RegularTimeSeries) = sum(x .^ 2) / duration(x)
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
function buffer(x::RegularTimeSeries, args...; kwargs...)
    y = _buffer(x, args...; kwargs...)
    t = _buffer(times(x), args...; kwargs...) .|> mean
    # For a regular time series, the buffer centres are regular
    ts = range(first(t), last(t), length(y))
    y = TimeSeries(ts, y)
end
window(x, n, p = n, args...; kwargs...) = buffer(x, n, n - p, args...; kwargs...)

function _delayembed(x::AbstractVector, n, τ, p = 1; kwargs...) # A delay embedding with dimension `n`, delay `τ`, and skip length of `p`
    y = window(x, n * τ, p; kwargs...)
    y = map(y) do _y
        @view _y[1:τ:end]
    end
end
delayembed(x::AbstractVector, args...; kwargs...) = _delayembed(x, args...; kwargs...)
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

function rectify(ts::DimensionalData.Dimension; tol = 6, zero = false)
    ts = collect(ts)
    origts = ts
    stp = ts |> diff |> mean
    err = ts |> diff |> std
    if err > stp / 10.0^(-tol)
        @warn "Step is not approximately constant, skipping rectification"
    else
        stp = round(stp; digits = tol)
        t0, t1 = round.(extrema(ts); digits = tol)
        if zero
            origts = t0:stp:(t1 + (10000 * stp))
            t1 = t1 - t0
            t0 = 0
        end
        ts = t0:stp:(t1 + (10000 * stp))
    end
    return ts, origts
end
rectifytime(ts::Ti; kwargs...) = rectify(ts; kwargs...)

function rectify(X::AbstractDimArray; dims, tol = 6, zero = false) # tol gives significant figures for rounding
    if !(dims isa Tuple || dims isa AbstractVector)
        dims = [dims]
    end
    for dim in dims
        ts, origts = rectify(DimensionalData.dims(X, dim); tol, zero)
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

function rectifytime(X::AbstractVector; tol = 6, zero = false)
    # Generate some common time indices as close as possible to the rectified times of each element of the input vector
    ts = times.(X)
    mint = maximum(minimum.(ts)) - exp10(tol) .. minimum(maximum.(ts)) + exp10(tol)
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

"""
    centraldiff!(x::RegularTimeSeries; dims=Ti, grad=-)

Compute the central difference of a regular time series `x`, in-place.
The first and last elements are set to the forward and backward difference, respectively.
The dimension to perform differencing over can be specified as `dims`, and the differencing function can be specified as `grad` (defaulting to the euclidean distance, `-`)
"""
centraldiff!(x::UnivariateRegular; kwargs...) = _centraldiff!(x; kwargs...)
function centraldiff!(x::MultidimensionalTimeSeries; dims = Ti, kwargs...)
    _centraldiff!(eachslice(x; dims); kwargs...)
end

"""
    centraldiff(x::RegularTimeSeries; dims=Ti, grad=-)

Compute the central difference of a regular time series `x`.
The first and last elements are set to the forward and backward difference, respectively.
The dimension to perform differencing over can be specified as `dims`, and the differencing function can be specified as `grad` (defaulting to the euclidean distance, `-`)
See [`centraldiff!`](@ref).
"""
function centraldiff(x::AbstractTimeSeries; kwargs...)
    y = deepcopy(x)
    centraldiff!(y; kwargs...)
    return y
end

"""
    centralderiv!(x::RegularTimeSeries; kwargs...)

Compute the central derivative of a regular time series `x`, in-place.
See [`centraldiff!`](@ref) for available keyword arguments.
"""
function centralderiv!(x::RegularTimeSeries; dims = Ti, kwargs...)
    if dims isa Tuple || dims isa AbstractVector
        error("Only one dimension can be specified for central derivatives.")
    end
    centraldiff!(x; dims, kwargs...)
    x ./= samplingperiod(x; dims)
    nothing
end

"""
    centralderiv(x::RegularTimeSeries)

Compute the central derivative of a regular time series `x`.
See [`centraldiff`](@ref) for available keyword arguments.
Also c.f. [`centralderiv!`](@ref).
"""
function centralderiv(x::RegularTimeSeries; dims = Ti, kwargs...)
    y = deepcopy(x)
    centralderiv!(y; dims, kwargs...)
    return y
end

function _rightdiff!(x; grad = -, dims = nothing) # Dims unused
    x[1:(end - 1)] .= grad(x[2:end], x[1:(end - 1)])
    # x[[1, end]] .= [grad(a, x[1]), grad(x[end], b)]
    x[[end]] .= [copy(x[end - 1])]
    return nothing
end

rightdiff!(x::UnivariateRegular; kwargs...) = _rightdiff!(x; kwargs...)
function rightdiff!(x::MultidimensionalTimeSeries; dims = Ti, kwargs...)
    _rightdiff!(eachslice(x; dims); kwargs...)
end
function rightdiff(x::AbstractTimeSeries; kwargs...)
    y = deepcopy(x)
    rightdiff!(y; kwargs...)
    return y
end
function rightderiv!(x::RegularTimeSeries; dims = Ti, kwargs...)
    if dims isa Tuple || dims isa AbstractVector
        error("Only one dimension can be specified for right derivatives.")
    end
    rightdiff!(x; dims, kwargs...)
    x ./= samplingperiod(x; dims)
    nothing
end

function rightderiv(x::RegularTimeSeries; dims = Ti, kwargs...)
    y = deepcopy(x)
    rightderiv!(y; dims, kwargs...)
    return y
end

function _leftdiff!(x; grad = -, dims = nothing) # Dims unused
    x[2:end] .= grad(x[2:end], x[1:(end - 1)])
    # x[[1, end]] .= [grad(a, x[1]), grad(x[end], b)]
    x[[1]] .= [copy(x[2])]
    return nothing
end

leftdiff!(x::UnivariateRegular; kwargs...) = _leftdiff!(x; kwargs...)
function leftdiff!(x::MultidimensionalTimeSeries; dims = Ti, kwargs...)
    _leftdiff!(eachslice(x; dims); kwargs...)
end
function leftdiff(x::AbstractTimeSeries; kwargs...)
    y = deepcopy(x)
    leftdiff!(y; kwargs...)
    return y
end
function leftderiv!(x::RegularTimeSeries; dims = Ti, kwargs...)
    if dims isa Tuple || dims isa AbstractVector
        error("Only one dimension can be specified for left derivatives.")
    end
    leftdiff!(x; dims, kwargs...)
    x ./= samplingperiod(x; dims)
    nothing
end

function leftderiv(x::RegularTimeSeries; dims = Ti, kwargs...)
    y = deepcopy(x)
    leftderiv!(y; dims, kwargs...)
    return y
end

Base.abs(x::AbstractTimeSeries) = Base.abs.(x)
Base.angle(x::AbstractTimeSeries) = Base.angle.(x)

# * See https://en.wikipedia.org/wiki/Directional_statistics
circresultant(θ; kwargs...) = mean(exp.(im .* θ); kwargs...)
circlength(θ; kwargs...) = abs.(circresultant(θ; kwargs...))
circmean(θ; kwargs...) = angle.(circresultant(θ; kwargs...))
circvar(θ; kwargs...) = 1 - circlength(θ; kwargs...)
circstd(θ; kwargs...) = sqrt.(-2 * log.(circlength(θ; kwargs...)))

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
