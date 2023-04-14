export  AbstractTimeSeries, AbstractTS,
        UnivariateTimeSeries, UnivariateTS,
        MultivariateTimeSeries, MultivariateTS,
        RegularTimeSeries, RegularTS,
        IrregularTimeSeries, IrregularTS,
        TimeIndex, RegularIndex, RegularTimeIndex,
        TimeSeries

## Time-series types
TimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.TimeDim}
AbstractTimeSeries = AbstractTS = AbstractDimArray{T, N, <:TimeIndex, B} where {T, N, B}

UnivariateTimeSeries = UnivariateTS = AbstractTimeSeries{T, 1} where T
MultivariateTimeSeries = MultivariateTS = AbstractTimeSeries{T, 2} where T
abstract type VariableDim{T} <: DimensionalData.IndependentDim{T} end
DimensionalData.@dim Var VariableDim "Var"

RegularIndex = DimensionalData.Dimensions.LookupArrays.Sampled{T, R} where {T, R<:AbstractRange}
RegularTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.TimeDim{<:RegularIndex}}
RegularTimeSeries = RegularTS = AbstractDimArray{T, N, <:RegularTimeIndex, B} where {T, N, B}
IrregularTimeSeries = AbstractTimeSeries



## Methods
TimeSeries(t, x) = DimArray(x, (Ti(t),))
TimeSeries(t, v, x) = DimArray(x, (Ti(t), Var(v)))

convertconst(a, _) = a
