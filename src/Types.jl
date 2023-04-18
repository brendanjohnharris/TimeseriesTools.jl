export  AbstractTimeSeries, AbstractTS,
        UnivariateTimeSeries, UnivariateTS,
        MultivariateTimeSeries, MultivariateTS,
        RegularTimeSeries, RegularTS,
        IrregularTimeSeries, IrregularTS,
        TimeIndex, RegularIndex, RegularTimeIndex,
        TimeSeries, Timeseries, TS, Var

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
TimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
AbstractTimeSeries = AbstractTS = AbstractDimArray{T, N, <:TimeIndex, B} where {T, N, B}

"""
    UnivariateTimeSeries{T}

A type alias for a time series with one variable (a vector with only a `Ti` dimension).
"""
UnivariateTimeSeries = UnivariateTS = AbstractTimeSeries{T, 1} where T

"""
    MultivariateTimeSeries{T}

A type alias for a multivariate time series (A matrix, with a first `Ti` dimension and an arbitrary second dimension).
"""
MultivariateTimeSeries = MultivariateTS = AbstractTimeSeries{T, 2} where T

abstract type VariableDim{T} <: DimensionalData.IndependentDim{T} end
DimensionalData.@dim Var VariableDim "Var"

"""
    Var

A DimensionalData.jl dimension representing the variables of a multivariate time series.
"""
Var

"""
    RegularIndex

A type alias for a regularly sampled dimension, wrapping an `AbstractRange`.
"""
RegularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R} where {T, R<:AbstractRange}

"""
    RegularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
RegularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
RegularTimeSeries = RegularTS = AbstractDimArray{T, N, <:RegularTimeIndex, B} where {T, N, B}

"""
    IrregularTimeSeries

A type alias for a potentially irregularly sampled time series.
"""
IrregularTimeSeries = AbstractTimeSeries


"""
    TimeSeries(t, x)

Constructs a univariate time series with time `t` and data `x`. Alteratively, use `TS(t, x)`

## Examples
```@example 1
julia> using TimeseriesTools, Unitful;
julia> t = 1:100
julia> x = rand(100)
julia> ts = TimeSeries(t, x)
julia> ts isa Union{UnivariateTimeSeries, RegularTimeSeries}
```
"""
TimeSeries(t, x) = DimArray(x, (Ti(t),))

"""
    TimeSeries(t, v, x)

Constructs a multivariate time series with time t, variable v, and data x.

## Examples
```@example 1
julia> t = 1:100;
julia> v = [:a, :b, :c];
julia> x = rand(100, 3);
julia> mts = TimeSeries(t, v, x)
julia> mts isa Union{MultivariateTimeSeries, RegularTimeSeries}
```
"""
TimeSeries(t, v, x) = DimArray(x, (Ti(t), Var(v)))

TS = Timeseries = TimeSeries

convertconst(a, _) = a
