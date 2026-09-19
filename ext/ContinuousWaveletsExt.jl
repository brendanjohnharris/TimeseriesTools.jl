module ContinuousWaveletsExt
using ContinuousWavelets
using IntervalSets
using DimensionalData
using TimeseriesTools
import TimeseriesTools: _waveletfreqs, _waveletspectrogram, waveletspectrogram
import TimeseriesTools.TimeseriesBase.Spectra: RegularSpectrogram

function _waveletfreqs(t; moth = Morlet(2π), β = 1, Q = 32)
    n = length(t)
    fs = 1.0 ./ step(t) # Assume rectified time dim
    W = ContinuousWavelets.computeWavelets(n, wavelet(moth; β, Q))[1]
    freqs = getMeanFreq(W, fs)
    freqs[1] = 0
    return freqs
end
function waveletfreqs(t; pass = nothing, kwargs...)
    freqs = _waveletfreqs(t; kwargs...)
    isnothing(pass) && return freqs
    return freqs[freqs .∈ [ClosedInterval(0, maximum(pass))]]
end

# `cwt` reads the compute device off the outermost array wrapper, so it rejects views, reshapes
# and adjoints against a dense daughter matrix. `copy` materialises them without leaving the
# device: a slice of a `CuArray` copies to a `CuArray`.
_dense(x::DenseArray) = x
_dense(x::AbstractArray) = copy(x)

# `cwt` transforms along the first axis and batches over the rest, returning time x frequency x
# (remaining axes).
function _waveletspectrogram(x::AbstractArray; moth, β, Q) # β = 1 means linear in log space
    c = wavelet(moth; β, Q)
    return ContinuousWavelets.cwt(_dense(x), c)
end

function _waveletspectrogram(t, x::AbstractArray; pass = nothing, moth, β, Q)
    if isnothing(pass)
        return _waveletspectrogram(x; moth, β, Q)
    end
    n = size(x, 1)
    @assert length(t) == n
    c = wavelet(moth; β, Q)
    W = ContinuousWavelets.computeWavelets(n, c)[1]
    freqs = getMeanFreq(W, 1.0 ./ step(t))
    keep = freqs .∈ [ClosedInterval(0, maximum(pass))]
    res = ContinuousWavelets.cwt(_dense(x), c, W[:, keep])
    return res[:, keep, ntuple(_ -> Colon(), ndims(res) - 2)...]
end

function _waveletspectrogram(
        x::RegularTimeseries; moth = Morlet(2π), β = 1,
        Q = 32,
        pass = nothing
    )::RegularSpectrogram
    t = times(x)
    res = _waveletspectrogram(t, x.data; moth, β, Q, pass)
    freqs = waveletfreqs(t; moth, β, Q, pass)
    return Timeseries(
        res, t, 𝑓(freqs), dims(x)[2:end]...;
        metadata = DimensionalData.metadata(x),
        refdims = refdims(x)
    )
end

_waveletspectrogram(x, s::Symbol; kwargs...) = _waveletspectrogram(x, Val(s); kwargs...)

function waveletspectrogram(x::RegularTimeseries, args...; kwargs...)::RegularSpectrogram
    return _waveletspectrogram(x, args...; kwargs...)
end

end # module
