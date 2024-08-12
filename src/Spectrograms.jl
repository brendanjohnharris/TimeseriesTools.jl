import DimensionalData: Dimension, TimeDim
export AbstractSpectrogram, MultivariateSpectrogram, RegularSpectrogram

const TimeFreqIndex = Tuple{T, F, Vararg{Dimension}} where {T <: ToolsTimeDim, F <: Freq}
const RegularTimeFreqIndex = Tuple{T, F,
                                   Vararg{Dimension}} where {
                                                             T <:
                                                             ToolsTimeDim{<:RegularIndex},
                                                             F <: Freq}

const AbstractSpectrogram = AbstractToolsArray{T, N, <:TimeFreqIndex, B} where {T, N, B}
times(x::AbstractSpectrogram) = dims(x, 𝑡).val.data
freqs(x::AbstractSpectrogram) = dims(x, Freq).val.data

const MultivariateSpectrogram = AbstractSpectrogram{T, 3} where {T}

const RegularSpectrogram = AbstractToolsArray{T, N, <:RegularTimeFreqIndex,
                                              B} where {T, N, B}
