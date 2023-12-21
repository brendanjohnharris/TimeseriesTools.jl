export AbstractTimeSeries, AbstractTS,
       UnivariateTimeSeries, UnivariateTS,
       MultivariateTimeSeries, MultivariateTS,
       RegularTimeSeries, RegularTS,
       IrregularTimeSeries, IrregularTS,
       TimeIndex, RegularIndex, RegularTimeIndex,
       IrregularIndex, IrregularTimeIndex,
       TimeSeries, Timeseries, TS, Var,
       stitch,
       IrregularBinaryTimeSeries, SpikeTrain, spiketrain

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
const TimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}
                        } where {A <: DimensionalData.TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
const AbstractTimeSeries = AbstractTS = AbstractDimArray{T, N, <:TimeIndex, B
                                                         } where {T, N, B}

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
const RegularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R
                                                                     } where {T,
                                                                              R <:
                                                                              AbstractRange}

"""
    RegularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const RegularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}
                               } where {A <: DimensionalData.TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
const RegularTimeSeries = RegularTS = AbstractDimArray{T, N, <:RegularTimeIndex, B
                                                       } where {T, N, B}

"""
    IrregularIndex

A type alias for an irregularly sampled dimension, wrapping an `AbstractVector`.
"""
const IrregularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R
                                                                       } where {T,
                                                                                R <:
                                                                                AbstractVector
                                                                                }

"""
    IrregularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const IrregularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}
                                 } where {A <: DimensionalData.TimeDim{<:IrregularIndex}}

"""
    IrregularTimeSeries

A type alias for a potentially irregularly sampled time series.
"""
const IrregularTimeSeries = IrregularTS = AbstractDimArray{T, N, <:IrregularTimeIndex, B
                                                           } where {T, N, B}

"""
    BinaryTimeSeries

A type alias for a time series of bits.
"""
const BinaryTimeSeries = SpikeTrain = BinaryTS = AbstractDimArray{T, N, <:TimeIndex, B
                                                                  } where {T <: Bool, N, B}

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

"""
    TimeSeries(t, f::Function)

Construct a time series by mapping a function `f` over the time points `t`.
"""
TimeSeries(t, f::Function; kwargs...) = TimeSeries(t, f.(t), kwargs...)

const TS = Timeseries = TimeSeries

convertconst(a, _) = a

"""
    Base.cat(D::DimensionalData.Dimension, args...; kwargs...)
Concatenate the arrays given in `args...`, and give the resulting extra axis dimensions `D`.
Note that unlike `Base.cat` without the first `Dim` argument, this increments all existing dimensions greater than `dims` by one (so N n×n arrays concatenated at `dims=1` will result in an N×n×n array).
"""
function Base.cat(D::DimensionalData.Dimension, x::AbstractDimArray, args...;
                  dims = nothing,
                  kwargs...)
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
    _args = rf.(getfield.(args, [:data]))
    x′ = cat(_x, _args...; dims, kwargs...)
    ds = Vector{Any}([DimensionalData.dims(x)...])
    insert!(ds, dims, D)
    y = DimArray(x′, (ds...,); refdims = refdims(x), name = name(x),
                 metadata = metadata(x))
    # if hasdim(y, Ti)
    #     ts = times(y)
    #     y = set(y, Ti => ts .- minimum(ts))
    # end
    return y
end

const UnivariateRegular = typeintersect(UnivariateTimeSeries, RegularTimeSeries)
const MultivariateRegular = typeintersect(MultivariateTimeSeries, RegularTimeSeries)

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
