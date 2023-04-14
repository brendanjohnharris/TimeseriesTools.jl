import DimensionalData.Dimensions.LookupArrays: At, Near
import DimensionalData.Dimensions.Dimension

export times, samplingrate, duration, samplingperiod

# Allow dims to be passed directly to selectors
Selectors = [:At, :Between, :Touches, :Near, :Where, :Contains]
[:($(S)(D::Dimension) = $(S)(D.val.data)) for S in Selectors] .|> eval

times(x::AbstractTimeSeries) = dims(x, Ti).val.data
Base.step(x::RegularTimeSeries) = x |> times |> step
samplingrate(x::RegularTimeSeries) = 1/step(x)
samplingperiod(x::RegularTimeSeries) = step(x)
duration(x::AbstractTimeSeries) = (last∘times)(x) - (first∘times)(x)
IntervalSets.Interval(x::AbstractTimeSeries) = (first∘times)(x)..(last∘times)(x)
