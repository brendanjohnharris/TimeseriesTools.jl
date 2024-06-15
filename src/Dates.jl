using Dates

export DateIndex, DateTimeIndex, DateTimeSeries

DateIndex = DateTIndex = Union{AbstractArray{<:Dates.AbstractTime},
                               AbstractRange{<:Dates.AbstractTime},
                               Tuple{<:Dates.AbstractTime}}

DateTimeIndex = Tuple{A,
                      Vararg{DimensionalData.Dimension}} where {A <:
                                                                DimensionalData.Dimension{<:DateIndex}}

DateTimeSeries = AbstractDimArray{T, N, <:DateTimeIndex, B} where {T, N, B}

unit(::DateTimeSeries) = NoUnits
