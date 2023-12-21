import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import DimensionalData: print_array, _print_array_ctx
import Normalization: NormUnion, AbstractNormalization

export times, samplingrate, duration, samplingperiod, UnitPower, dimname, dimnames,
       describedim, describedims, describename, interlace, _buffer, buffer, window,
       delayembed, centraldiff!, centraldiff, centralderiv!, centralderiv

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
Base.step(x::RegularTimeSeries) = x |> times |> step

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
samplingrate(x::RegularTimeSeries) = 1 / step(x)

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
samplingperiod(x::RegularTimeSeries) = step(x)

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

function rectifytime(X::IrregularTimeSeries; tol = 6, zero = false) # tol gives significant figures for rounding
    ts = times(X)
    stp = ts |> diff |> mean
    err = ts |> diff |> std
    if err > stp / 10.0^(-tol)
        @warn "Time step is not approximately constant, skipping rectification"
    else
        stp = round(stp; digits = tol)
        t0, t1 = round.(extrema(ts); digits = tol)
        if zero
            origts = t0:stp:(t1 + (10000 * stp))
            t1 = t1 - t0
            t0 = 0
        end
        ts = t0:stp:(t1 + (10000 * stp))
        ts = ts[1:size(X, Ti)] # Should be ok?
    end
    @assert length(ts) == size(X, Ti)
    X = set(X, Ti => ts)
    if zero
        X = rebuild(X; metadata = (:time => origts, pairs(metadata(X))...))
    end
    return X
end

function _centraldiff!(x)
    a = x[2] # Save here, otherwise they get mutated before we use them
    b = x[end - 1]
    x[2:(end - 1)] .= (x[3:end] - x[1:(end - 2)]) / 2
    x[[1, end]] .= [a - x[1], x[end] - b]
    return nothing
end

"""
    centraldiff!(x::RegularTimeSeries)

Compute the central difference of a regular time series `x`, in-place.
The first and last elements are set to the forward and backward difference, respectively.
"""
centraldiff!(x::UnivariateRegular) = _centraldiff!(x)
function centraldiff!(x::typeintersect(MultivariateTimeSeries, RegularTimeSeries))
    _centraldiff!(eachslice(x; dims = Ti))
end

"""
    centraldiff(x::RegularTimeSeries)

Compute the central difference of a regular time series `x`.
The first and last elements are set to the forward and backward difference, respectively.
See [`centraldiff!`](@ref).
"""
function centraldiff(x::AbstractTimeSeries)
    y = deepcopy(x)
    centraldiff!(y)
    return y
end

"""
    centralderiv!(x::RegularTimeSeries)

Compute the central derivative of a regular time series `x`, in-place.
See [`centraldiff!`](@ref).
"""
function centralderiv!(x::AbstractTimeSeries)
    centraldiff!(x)
    x ./= samplingperiod(x)
    nothing
end

"""
    centralderiv(x::RegularTimeSeries)

Compute the central derivative of a regular time series `x`.
See [`centraldiff`](@ref) and  [`centralderiv!`](@ref).
"""
function centralderiv(x::AbstractTimeSeries)
    y = deepcopy(x)
    centralderiv!(y)
    return y
end
