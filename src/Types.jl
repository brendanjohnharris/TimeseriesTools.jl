import DimensionalData: Dimension, TimeDim

export AbstractTimeSeries, AbstractTS,
       UnivariateTimeSeries, UnivariateTS,
       MultivariateTimeSeries, MultivariateTS,
       RegularTimeSeries, RegularTS,
       UnivariateRegular, MultivariateRegular,
       IrregularTimeSeries, IrregularTS,
       TimeIndex, RegularIndex, RegularTimeIndex,
       IrregularIndex, IrregularTimeIndex,
       TimeSeries, Timeseries, TS, Var,
       stitch,
       IrregularBinaryTimeSeries, SpikeTrain, spiketrain,
       MultidimensionalIndex, MultidimensionalTimeSeries, MultidimensionalTS

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
const TimeIndex = Tuple{A, Vararg{Dimension}} where {A <: TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
const AbstractTimeSeries = AbstractTS = AbstractDimArray{T, N, <:TimeIndex,
                                                         B} where {T, N, B}

"""
    UnivariateTimeSeries{T}

A type alias for a time series with one variable (a vector with only a `Ti` dimension).
"""
const UnivariateTimeSeries = UnivariateTS = AbstractTimeSeries{T, 1} where {T}

"""
    MultivariateTimeSeries{T}

A type alias for a multivariate time series (A matrix, with a first `Ti` dimension and an arbitrary second dimension).
"""
const MultivariateTimeSeries = MultivariateTS = AbstractTimeSeries{T, 2} where {T}

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
const RegularIndex = Dimensions.LookupArrays.Sampled{T, R} where {T, R <: AbstractRange}

"""
    RegularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const RegularTimeIndex = Tuple{A, Vararg{Dimension}} where {A <: TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
const RegularTimeSeries = RegularTS = AbstractDimArray{T, N, <:RegularTimeIndex,
                                                       B} where {T, N, B}

const MultidimensionalIndex = Tuple{A,
                                    Vararg{Dimension{B}}} where {
                                                                 A <:
                                                                 TimeDim{<:RegularIndex},
                                                                 B <:
                                                                 RegularIndex
                                                                 }

"""
A multidimensional time series has a regular sampling over a dimension other than time; a one-dimensional time series can be thought of as a field over an even grid in 1 dimension that fluctuates over time.
"""
const MultidimensionalTimeSeries = AbstractDimArray{T, N, <:MultidimensionalIndex,
                                                    B} where {T, N, B}
const MultidimensionalTS = MultidimensionalTimeSeries

"""
    IrregularIndex

A type alias for an irregularly sampled dimension, wrapping an `AbstractVector`.
"""
const IrregularIndex = Dimensions.LookupArrays.Sampled{T,
                                                       R} where {T,
                                                                 R <:
                                                                 AbstractVector
                                                                 }

"""
    IrregularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const IrregularTimeIndex = Tuple{A,
                                 Vararg{Dimension}} where {A <:
                                                           TimeDim{<:IrregularIndex}}

"""
    IrregularTimeSeries

A type alias for a potentially irregularly sampled time series.
"""
const IrregularTimeSeries = IrregularTS = AbstractDimArray{T, N, <:IrregularTimeIndex,
                                                           B} where {T, N, B}

"""
    BinaryTimeSeries

A type alias for a time series of bits.
"""
const BinaryTimeSeries = SpikeTrain = BinaryTS = AbstractDimArray{T, N, <:TimeIndex,
                                                                  B} where {T <: Bool, N, B}

"""
    SpikeTrain

A type alias for a spike-train time series, which contains spike times in the time dimension and `true` for all values corresponding to a spike. The spike times can be retrieved with `times(x)`.
"""
SpikeTrain

function spiketrain(x; kwargs...)
    TimeSeries(sort(x), trues(length(x)); kwargs...)
end

"""
    TimeSeries(t, x)

Constructs a univariate time series with time `t` and data `x`. Alteratively, use `TS(t, x)`

## Examples
```@example 1
julia> using TimeseriesTools, Unitful;
julia> t = 1:100
julia> x = rand(100)
julia> ts = TimeSeries(t, x)
julia> ts isa typeintersect(UnivariateTimeSeries, RegularTimeSeries)
```
"""
TimeSeries(t, x; kwargs...) = DimArray(x, (Ti(t),); kwargs...)
TimeSeries(t::TimeDim, x; kwargs...) = DimArray(x, (t,); kwargs...)

"""
    TimeSeries(t, v, x)

Constructs a multivariate time series with time t, variable v, and data x.

## Examples
```@example 1
julia> t = 1:100;
julia> v = [:a, :b, :c];
julia> x = rand(100, 3);
julia> mts = TimeSeries(t, v, x)
julia> mts isa typeintersect(MultivariateTimeSeries, RegularTimeSeries)
```
"""
function TimeSeries(t::TimeDim, v::Dimension, x; kwargs...)
    DimArray(x, (t, v); kwargs...)
end
function TimeSeries(t::TimeDim, v, x; kwargs...)
    DimArray(x, (t, Var(v)); kwargs...)
end
function TimeSeries(t, v::Dimension, x; kwargs...)
    DimArray(x, (Ti(t), v); kwargs...)
end
TimeSeries(t, v, x; kwargs...) = DimArray(x, (Ti(t), Var(v)); kwargs...)

function TimeSeries(t::TimeDim, a::Dimension,
                    b::Dimension, x; kwargs...)
    DimArray(x, (t, a, b); kwargs...)
end
function TimeSeries(t, a::Dimension,
                    b::Dimension, x; kwargs...)
    DimArray(x, (Ti(t), a, b); kwargs...)
end

import DimensionalData.data
function data(x::AbstractTimeSeries)
    _x = x.data
    while _x isa AbstractTimeSeries
        _x = _x.data
    end
    return _x
end
function TimeSeries(t, x::AbstractDimArray; kwargs...)
    TimeSeries(t, DimensionalData.data(x); kwargs...)
end
function TimeSeries(t, v, x::AbstractDimArray; kwargs...)
    TimeSeries(t, v, DimensionalData.data(x); kwargs...)
end

"""
    TimeSeries(t, f::Function)

Construct a time series by mapping a function `f` over the time points `t`.
"""
TimeSeries(t, f::Function; kwargs...) = TimeSeries(t, f.(t), kwargs...)

const TS = Timeseries = TimeSeries

convertconst(a, _) = a

const UnivariateRegular = typeintersect(UnivariateTimeSeries, RegularTimeSeries)
const MultivariateRegular = typeintersect(MultivariateTimeSeries, RegularTimeSeries)
