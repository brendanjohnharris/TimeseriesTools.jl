# import DimensionalData: Dimension, TimeDim
# export AbstractSpectrogram, MultivariateSpectrogram, RegularSpectrogram

# const TimeFreqIndex = Tuple{T, F, Vararg{Dimension}} where {T <: TimeDim, F <: 𝑓}
# const RegularTimeFreqIndex = Tuple{T, F,
#                                    Vararg{Dimension}} where {
#                                                              T <:
#                                                              TimeDim{<:RegularIndex},
#                                                              F <: 𝑓}

# const AbstractSpectrogram = AbstractToolsArray{T, N, <:TimeFreqIndex, B} where {T, N, B}
# times(x::AbstractSpectrogram) = dims(x, 𝑡).val.data
# freqs(x::AbstractSpectrogram) = dims(x, 𝑓).val.data

# const MultivariateSpectrogram = AbstractSpectrogram{T, 3} where {T}

# const RegularSpectrogram = AbstractToolsArray{T, N, <:RegularTimeFreqIndex,
#                                               B} where {T, N, B}
