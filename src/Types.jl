export AbstractTimeSeries, AbstractTS,
       UnivariateTimeSeries, UnivariateTS,
       MultivariateTimeSeries, MultivariateTS,
       RegularTimeSeries, RegularTS,
       IrregularTimeSeries, IrregularTS,
       TimeIndex, RegularIndex, RegularTimeIndex,
       IrregularIndex, IrregularTimeIndex,
       TimeSeries, Timeseries, TS, Var,
       stitch,
       IrregularBinaryTimeSeries, SpikeTrain

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
TimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A <: DimensionalData.TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
AbstractTimeSeries = AbstractTS = AbstractDimArray{T, N, <:TimeIndex, B} where {T, N, B}

"""
    UnivariateTimeSeries{T}

A type alias for a time series with one variable (a vector with only a `Ti` dimension).
"""
UnivariateTimeSeries = UnivariateTS = AbstractTimeSeries{T, 1} where {T}

"""
    MultivariateTimeSeries{T}

A type alias for a multivariate time series (A matrix, with a first `Ti` dimension and an arbitrary second dimension).
"""
MultivariateTimeSeries = MultivariateTS = AbstractTimeSeries{T, 2} where {T}

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
RegularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R
                                                               } where {T,
                                                                        R <: AbstractRange}

"""
    RegularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
RegularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}
                         } where {A <: DimensionalData.TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
RegularTimeSeries = RegularTS = AbstractDimArray{T, N, <:RegularTimeIndex, B
                                                 } where {T, N, B}

"""
    IrregularIndex

A type alias for an irregularly sampled dimension, wrapping an `AbstractVector`.
"""
IrregularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R
                                                                 } where {T,
                                                                          R <:
                                                                          AbstractVector}

"""
    IrregularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
IrregularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}
                           } where {A <: DimensionalData.TimeDim{<:IrregularIndex}}

"""
    IrregularTimeSeries

A type alias for a potentially irregularly sampled time series.
"""
IrregularTimeSeries = IrregularTS = AbstractDimArray{T, N, <:IrregularTimeIndex, B
                                                     } where {T, N, B}

"""
    BinaryTimeSeries

A type alias for a time series of bits.
"""
BinaryTimeSeries = BinaryTS = AbstractDimArray{<:Bool, N, <:TimeIndex, B
                                               } where {N, B}

IrregularBinaryTimeSeries = SpikeTrain = typeintersect(BinaryTimeSeries,
                                                       IrregularTimeSeries)

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
TimeSeries(t, v, x; kwargs...) = DimArray(x, (Ti(t), Var(v)); kwargs...)
function TimeSeries(t, v::DimensionalData.Dimension, x; kwargs...)
    DimArray(x, (Ti(t), v); kwargs...)
end

TS = Timeseries = TimeSeries

convertconst(a, _) = a

"""
    Base.cat(D::DimensionalData.Dimension, args...; kwargs...)
Concatenate the arrays given in `args...`, and give the resulting extra axis dimensions `D`.
"""
function Base.cat(D::DimensionalData.Dimension, x::AbstractDimArray, args...; kwargs...)
    x′ = cat(x.data, getfield.(args, [:data])...; kwargs...)
    y = DimArray(x′, (dims(x)..., D); refdims = refdims(x), name = name(x),
                 metadata = metadata(x))
    ts = times(y)
    set(y, Ti => ts .- minimum(ts))
end

UnivariateRegular = typeintersect(UnivariateTimeSeries, RegularTimeSeries)
MultivariateRegular = typeintersect(MultivariateTimeSeries, RegularTimeSeries)

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
    @assert(dims(x)[2:end].==dims(y)[2:end])
    z = vcat(x.data, y.data)
    z = TimeSeries(dt:dt:(dt * size(z, 1)), dims(x)[2:end]..., z)
end
stitch(X, Y, args...) = reduce(stitch, (X, Y, args...))
