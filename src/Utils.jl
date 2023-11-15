import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension
import Normalization: NormUnion, AbstractNormalization

export times, samplingrate, duration, samplingperiod, UnitPower, dimname, dimnames,
       describedim, describedims, describename, interlace, _buffer, buffer

import LinearAlgebra.mul!
function mul!(a::AbstractVector, b::AbstractTimeSeries, args...; kwargs...)
    mul!(a, b.data, args...; kwargs...)
end

Selectors = [:At, :Between, :Touches, :Near, :Where, :Contains]
# Allow dims to be passed directly to selectors
[:($(S)(D::Dimension) = $(S)(D.val.data)) for S in Selectors] .|> eval

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

function buffer(x::UnivariateTimeSeries, n, p = 0, discard = true)
    y = [@view x[i:min(i + n - 1, end)] for i in 1:(n - p):length(x)]
    if discard && length(y[end]) < n
        pop!(y)
    end
    y
end
function buffer(x::MultivariateTimeSeries, n, p = 0, discard = true)
    y = [@view x[i:min(i + n - 1, end), :] for i in 1:(n - p):size(x, 1)]
    if discard && size(y[end], 1) < n
        pop!(y)
    end
    y
end
