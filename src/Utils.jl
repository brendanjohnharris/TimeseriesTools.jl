import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension

export times, samplingrate, duration, samplingperiod

Selectors = [:At, :Between, :Touches, :Near, :Where, :Contains]
# Allow dims to be passed directly to selectors
[:($(S)(D::Dimension) = $(S)(D.val.data)) for S in Selectors] .|> eval

"""
    times(x::AbstractTimeSeries)

Returns the time indices of the [`AbstractTimeSeries`](@ref) `x`.

## Examples
```jldoctest
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
```jldoctest
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
```jldoctest
julia> t = 1:100;
julia> x = rand(100);
julia> rts = TimeSeries(t, x);
julia> samplingrate(rts) == 1
```
"""
samplingrate(x::RegularTimeSeries) = 1/step(x)

"""
    samplingperiod(x::RegularTimeSeries)

Returns the sampling period (step size) of a regularly sampled [`RegularTimeSeries`](@ref).

## Examples
```jldoctest
julia> t = 1:100;
julia> x = rand(100);
julia> rts = TimeSeries(t, x);
julia> samplingperiod(rts) == 1
```
"""
samplingperiod(x::RegularTimeSeries) = step(x)

"""
    duration(x::AbstractTimeSeries)

Returns the duration of the [@ref](AbstractTimeSeries) `x`.

## Examples
```jldoctest
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x);
julia> duration(ts) == 99
```
"""
duration(x::AbstractTimeSeries) = (last∘times)(x) - (first∘times)(x)

"""
    IntervalSets.Interval(x::AbstractTimeSeries)

Returns an interval representing the range of the [`AbstractTimeSeries`](@ref) `x`.

## Examples
```jldoctest
julia> using IntervalSets;
julia> t = 1:100;
julia> x = rand(100);
julia> ts = TimeSeries(t, x);
julia> IntervalSets.Interval(ts) == (1..100)
```
"""
IntervalSets.Interval(x::AbstractTimeSeries) = (first∘times)(x)..(last∘times)(x)
