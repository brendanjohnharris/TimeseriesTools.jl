import DimensionalData: Dimension, TimeDim
export AbstractSpectrogram, MultivariateSpectrogram, RegularSpectrogram

const TimeFreqIndex = Tuple{T, F, Vararg{Dimension}} where {T <: TimeDim, F <: Freq}
const RegularTimeFreqIndex = Tuple{T, F,
                                   Vararg{Dimension}} where {T <: TimeDim{<:RegularIndex},
                                                             F <: Freq}

const AbstractSpectrogram = AbstractDimArray{T, N, <:TimeFreqIndex, B} where {T, N, B}
times(x::AbstractSpectrogram) = dims(x, Ti).val.data
freqs(x::AbstractSpectrogram) = dims(x, Freq).val.data

const MultivariateSpectrogram = AbstractSpectrogram{T, 3} where {T}

const RegularSpectrogram = AbstractDimArray{T, N, <:RegularTimeFreqIndex, B} where {T, N, B}
