import DimensionalData: Dimension, TimeDim, NoName, NoMetadata, format

export AbstractToolsArray, ToolsArray,
       AbstractTimeSeries, AbstractTS,
       UnivariateTimeSeries, UnivariateTS,
       MultivariateTimeSeries, MultivariateTS,
       RegularTimeSeries, RegularTS,
       UnivariateRegular, MultivariateRegular,
       IrregularTimeSeries, IrregularTS,
       TimeIndex, RegularIndex, RegularTimeIndex,
       IrregularIndex, IrregularTimeIndex,
       TimeSeries, Timeseries, TS,
       stitch,
       IrregularBinaryTimeSeries, SpikeTrain, spiketrain,
       MultidimensionalIndex, MultidimensionalTimeSeries, MultidimensionalTS,
       ToolsDimension, ToolsDim,
       𝑡, 𝑥, 𝑦, 𝑧, 𝑓, Var

"""
A local type to avoid overloading and piracy issues with DimensionalData.jl
"""
abstract type AbstractToolsArray{T, N, D, A} <: DimensionalData.AbstractDimArray{T, N, D, A} end

AbstractDimVector = AbstractToolsArray{T, 1} where {T}
AbstractDimMatrix = AbstractToolsArray{T, 2} where {T}

struct ToolsArray{T, N, D <: Tuple, R <: Tuple, A <: AbstractArray{T, N}, Na, Me} <:
       AbstractToolsArray{T, N, D, A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
end

function ToolsArray(data::A, dims::Tuple{D, Vararg};
                    refdims::R = (), name::Na = NoName(),
                    metadata::M = NoMetadata()) where {D <: DimensionalData.Dimension,
                                                       R, A, Na, M}
    ToolsArray(data, format(dims, data), refdims, name, metadata)
end

ToolsArray(x::AbstractArray, D::DimensionalData.Dimension) = ToolsArray(x, (D,))

function ToolsArray(D::DimensionalData.DimArray)
    ToolsArray(D.data, D.dims, D.refdims, D.name, D.metadata)
end

@inline function DimensionalData.rebuild(A::ToolsArray, data::AbstractArray, dims::Tuple,
                                         refdims::Tuple, name, metadata)
    ToolsArray(data, dims, refdims, name, metadata)
end

# * Custom dimensions. Arrays with these dimensions will default to being reconstructed as
#   ToolsArrays
import DimensionalData: TimeDim, XDim, YDim, ZDim
DimensionalData.@dim 𝑡 TimeDim "Time"
DimensionalData.@dim 𝑥 XDim "x"
DimensionalData.@dim 𝑧 YDim "y"
DimensionalData.@dim 𝑦 ZDim "z"

abstract type VariableDim{T} <: Dimension{T} end
DimensionalData.@dim Var VariableDim "Var"

abstract type FrequencyDim{T} <: Dimension{T} end
DimensionalData.@dim 𝑓 FrequencyDim "Frequency"

abstract type ToolsDim{T} <: DimensionalData.Dimension{T} end
ToolsDimension = Union{𝑡, 𝑥, 𝑧, 𝑦, 𝑓, Var, ToolsDim}

function DimensionalData.dimconstructor(::Tuple{ToolsDimension,
                                                Vararg{DimensionalData.Dimension}})
    ToolsArray
end
DimensionalData.dimconstructor(::Tuple{<:ToolsDimension, Vararg}) = ToolsArray
DimensionalData.dimconstructor(dims::ToolsDimension) = ToolsArray

# struct Dim{S, T} <: ToolsDimension{T}
#     val::T
#     function Dim{S}(val; kw...) where {S}
#         if length(kw) > 0
#             val = DimensionalData.AutoVal(val, values(kw))
#         end
#         new{S, typeof(val)}(val)
#     end
#     function Dim{S}(val::AbstractArray; kw...) where {S}
#         if length(kw) > 0
#             val = DimensionalData.AutoLookup(val, values(kw))
#         end
#         Dim{S, typeof(val)}(val)
#     end
#     function Dim{S, T}(val::T) where {S, T}
#         new{S, T}(val)
#     end
# end
# Dim{S}() where {S} = Dim{S}(:)

# DimensionalData.name(::Type{<:Dim{S}}) where {S} = S
# DimensionalData.basetypeof(::Type{<:Dim{S}}) where {S} = Dim{S}
# const TDim = ToolsDim

# abstract type TimeDim{T} <: ToolsDimension{T} end
# DimensionalData.@dim 𝑡 TimeDim "𝑡"
# DimensionalData.@dim 𝑥 ToolsDimension "𝑥"
# DimensionalData.@dim 𝑦 ToolsDimension "𝑦"
# DimensionalData.@dim 𝑧 ToolsDimension "𝑧"

# function DimensionalData.rebuildsliced(f::Function, A::AbstractToolsArray,
#                                        data::AbstractArray, I::Tuple, name = name(A))
#     DimensionalData.rebuildsliced(f, DimArray(A), data, I, name) |> ToolsArray
# end

# function DimensionalData._similar(A::AbstractToolsArray, T::Type, shape::Tuple)
#     data = similar(parent(A), T, map(DimensionalData._parent_range, shape))
#     shape isa Tuple{Vararg{DimensionalData.Dimensions.DimUnitRange}} || return data
#     return ToolsArray(data, dims(shape))
# end
# function DimensionalData._similar(::Type{T}, shape::Tuple) where {T <: AbstractToolsArray}
#     data = similar(T, map(DimensionalData._parent_range, shape))
#     shape isa Tuple{Vararg{DimensionalData.Dimensions.DimUnitRange}} || return data
#     return ToolsArray(data, dims(shape))
# end
# function Base.similar(A::DimensionalData.AbstractDimArrayGenerator, ::Type{T},
#                       D::DimensionalData.DimTuple) where {T <: AbstractToolsArray}
#     ToolsArray(D)(A; data = similar(Array{T}, size(D)), dims = D, refdims = (),
#                   metadata = NoMetadata())
# end
# function Base.similar(A::DimensionalData.AbstractDimArrayGenerator, ::Type{T},
#                       D::Tuple{}) where {T <: AbstractToolsArray}
#     ToolsArray(D)(A; data = similar(Array{T}, ()), dims = (), refdims = (),
#                   metadata = NoMetadata())
# end

TimeSeries(x::DimArray) = ToolsArray(x)

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
const TimeIndex = Tuple{A, Vararg{Dimension}} where {A <: TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
const AbstractTimeSeries = AbstractTS = AbstractToolsArray{T, N, <:TimeIndex,
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
const RegularTimeIndex = Tuple{A,
                               Vararg{Dimension}} where {A <: TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
const RegularTimeSeries = RegularTS = AbstractToolsArray{T, N, <:RegularTimeIndex,
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
const MultidimensionalTimeSeries = AbstractToolsArray{T, N, <:MultidimensionalIndex,
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
const IrregularTimeSeries = IrregularTS = AbstractToolsArray{T, N, <:IrregularTimeIndex,
                                                             B} where {T, N, B}

"""
    BinaryTimeSeries

A type alias for a time series of bits.
"""
const BinaryTimeSeries = SpikeTrain = BinaryTS = AbstractToolsArray{T, N, <:TimeIndex,
                                                                    B} where {T <: Bool, N,
                                                                              B}

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
TimeSeries(t, x; kwargs...) = ToolsArray(x, (𝑡(t),); kwargs...)
TimeSeries(t::TimeDim, x; kwargs...) = ToolsArray(x, (t,); kwargs...)

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
    ToolsArray(x, (t, v); kwargs...)
end
function TimeSeries(t::TimeDim, v, x; kwargs...)
    ToolsArray(x, (t, Var(v)); kwargs...)
end
function TimeSeries(t, v::Dimension, x; kwargs...)
    ToolsArray(x, (𝑡(t), v); kwargs...)
end
TimeSeries(t, v, x; kwargs...) = ToolsArray(x, (𝑡(t), Var(v)); kwargs...)

function TimeSeries(t::TimeDim, a::Dimension,
                    b::Dimension, x; kwargs...)
    ToolsArray(x, (t, a, b); kwargs...)
end
function TimeSeries(t, a::Dimension,
                    b::Dimension, x; kwargs...)
    ToolsArray(x, (𝑡(t), a, b); kwargs...)
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
