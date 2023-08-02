using Dates

export DateIndex, DateTimeIndex, DateTimeSeries

if false
x = 1:100
t = DateTime(1901):Year(1):DateTime(2000);
y = TimeSeries(t, x)

samplingperiod(y)
times(y)
duration(y)
end

DateIndex = DateTIndex = Union{AbstractArray{<:Dates.AbstractTime}, AbstractRange{<:Dates.AbstractTime}, Tuple{<:Dates.AbstractTime}}

DateTimeIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:DimensionalData.Dimension{<:DateIndex}}

DateTimeSeries = AbstractDimArray{T, N, <:DateTimeIndex, B} where {T, N, B}

unit(x::DateTimeSeries) = NoUnits
