module ContinuousWaveletsExt
using ContinuousWavelets
using Mmap
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

function _waveletspectrogram(x::AbstractVector; moth, β, Q) # β = 1 means linear in log space
    c = wavelet(moth; β, Q)
    return res = ContinuousWavelets.cwt(x, c)
end

function _waveletspectrogram(t, x::AbstractVector; pass = nothing, moth, β, Q)
    if isnothing(pass)
        return _waveletspectrogram(x; moth, β, Q)
    end
    n = size(x, 1)
    @assert length(t) == n
    c = wavelet(moth; β, Q)
    W = ContinuousWavelets.computeWavelets(n, c)[1]
    freqs = getMeanFreq(W, 1.0 ./ step(t))
    pass = ClosedInterval(0, maximum(pass))
    W = W[:, freqs .∈ [pass]]
    return res = ContinuousWavelets.cwt(x, c, W)[:, freqs .∈ [pass]]
end

function _waveletspectrogram(
        x::RegularTimeseries; moth = Morlet(2π), β = 1,
        Q = 32,
        pass = nothing
    )::RegularSpectrogram
    t = times(x)
    res = _waveletspectrogram(t, x.data; moth, β, Q, pass)
    freqs = waveletfreqs(t; moth, β, Q, pass)
    return res = Timeseries(
        res, t, 𝑓(freqs); metadata = DimensionalData.metadata(x),
        refdims = refdims(x)
    )
end

# function _waveletspectrogram(x::RegularTimeseries, ::Val{:mmap}; window = 50000,
#                            kwargs...)::RegularSpectrogram
#     md = DimensionalData.metadata(x)
#     rd = DimensionalData.refdims(x)
#     𝓍 = _slidingwindow(x, window; tail = :overlap)
#     t = dims(x, 𝑡)
#     e = step(t) / 2
#     freqs = waveletfreqs(dims(𝓍[1], 𝑡); kwargs...)
#     sz = (length(t), length(freqs))
#     fname = tempname()
#     s = open(fname, "w+")
#     write.((s,), sz)
#     W = mmap(s, Matrix{ComplexF32}, sz)
#     res = DimArray(W, (t, 𝑓(freqs)); metadata = (; md..., file = fname),
#                    refdims = rd)
#     threadlog, threadmax = (0, length(𝓍))
#     @withprogress name="Wavelet transform" begin
#         for _x in 𝓍
#             subres = _waveletspectrogram(_x; kwargs...)
#             tx = extrema(dims(subres, 1))
#             fx = extrema(dims(subres, 2))
#             tilims = Interval{:closed, :closed}(tx[1] - e, tx[2] + e)
#             flims = Interval{:closed, :closed}(fx[1] - e, fx[2] + e)
#             res[ 𝑡(tilims), 𝑓(flims)] .= subres
#             if threadmax > 1
#                 Threads.threadid() == 1 && (threadlog += 1) % 1 == 0 &&
#                     @logprogress threadlog / threadmax
#             end
#         end
#     end
#     close(s)
#     return res
# end

_waveletspectrogram(x, s::Symbol; kwargs...) = _waveletspectrogram(x, Val(s); kwargs...)

function waveletspectrogram(x::RegularTimeseries, args...; kwargs...)::RegularSpectrogram
    return _waveletspectrogram(x, args...; kwargs...)
end
function waveletspectrogram(x::MultivariateTimeseries, args...; kwargs...)
    f = x -> waveletspectrogram(x)
    return ToolsArray(map(f, eachslice(x; dims = 2)), dims(x, 2)) |> stack
end

# ! Add normalized energy spectra and power spectra for wavelet transform

end # module
