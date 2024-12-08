module AutocorrelationsExt
using TimeseriesTools
using Autocorrelations
using DimensionalData
using LinearAlgebra

import Autocorrelations: fftacf, fftacf!, dotacf, dotacf!, default_lags
import TimeseriesTools.msdist

function fftacf!(r::AbstractVector, x::RegularTimeSeries, args...; kwargs...)
    fftacf!(r, parent(x), args...; kwargs...)
end
function fftacf!(r::AbstractVector, x::IrregularTimeSeries, args...; kwargs...)
    throw(MethodError(fftacf!, (r, x, args...)))
end
function dotacf!(r::AbstractVector, x::RegularTimeSeries, args...; kwargs...)
    dotacf!(r, parent(x), args...; kwargs...)
end
function dotacf!(r::AbstractVector, x::IrregularTimeSeries, args...; kwargs...)
    throw(MethodError(dotacf!, (r, x, args...)))
end
function fftacf(x::RegularTimeSeries, lags = default_lags(x); kwargs...)
    TimeSeries(lags .* samplingperiod(x), fftacf(parent(x), lags; kwargs...))
end
function fftacf(x::IrregularTimeSeries, args...; kwargs...)
    throw(MethodError(fftacf, (x, args...)))
end
function dotacf(x::RegularTimeSeries, lags = default_lags(x); kwargs...)
    TimeSeries(lags .* samplingperiod(x), dotacf(parent(x), lags; kwargs...))
end
function dotacf(x::IrregularTimeSeries, args...; kwargs...)
    throw(MethodError(dotacf, (x, args...)))
end

"""
Inspired by https://github.com/mastrof/MeanSquaredDisplacement.jl/
"""
function msdist(x::AbstractVector,
                lags = range(0, length(x) - 1, step = 1))
    if !issorted(lags)
        throw(ArgumentError("Lags must be sorted"))
    end
    l = length(x)
    S2 = acf(x, lags)
    D = [0.0; dot.(x, x); 0.0]
    Q = 2 * sum(D)
    S1 = similar(S2)
    for k in 0:lags[end]
        Q -= (D[k - 1 + 2] + D[l - k + 2]) # Maps indices from -1 to +1
        if k ∈ lags
            i = searchsortedfirst(lags, k)
            S1[i] = Q / (l - k) - 2S2[i]
        end
    end
    return S1
end
function msdist(x::UnivariateRegular, lags = range(0, length(x) - 1, step = 1))
    return TimeSeries(lags .* samplingperiod(x), msdist(parent(x), lags))
end
function msdist(x::MultivariateRegular,
                lags = range(0, size(x, 1) - 1, step = 1))
    d = dims(x)[2:end]
    m = mapslices(x -> msdist(x, lags), parent(x); dims = 1)
    return TimeSeries(lags .* samplingperiod(x), d..., m)
end

end
