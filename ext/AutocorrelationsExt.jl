module AutocorrelationsExt
using TimeseriesTools
using Autocorrelations
using DimensionalData
using LinearAlgebra

import Autocorrelations: fftacf, fftacf!, dotacf, dotacf!, default_lags
import TimeseriesTools.msdist

function fftacf!(r::AbstractVector, x::RegularTimeseries, args...; kwargs...)
    fftacf!(r, parent(x), args...; kwargs...)
end
function fftacf!(r::AbstractVector, x::IrregularTimeseries, args...; kwargs...)
    throw(MethodError(fftacf!, (r, x, args...)))
end
function dotacf!(r::AbstractVector, x::RegularTimeseries, args...; kwargs...)
    dotacf!(r, parent(x), args...; kwargs...)
end
function dotacf!(r::AbstractVector, x::IrregularTimeseries, args...; kwargs...)
    throw(MethodError(dotacf!, (r, x, args...)))
end
function fftacf(x::RegularTimeseries, lags = default_lags(x); kwargs...)
    Timeseries(fftacf(parent(x), lags; kwargs...), lags .* samplingperiod(x))
end
function fftacf(x::IrregularTimeseries, args...; kwargs...)
    throw(MethodError(fftacf, (x, args...)))
end
function dotacf(x::RegularTimeseries, lags = default_lags(x); kwargs...)
    Timeseries(dotacf(parent(x), lags; kwargs...), lags .* samplingperiod(x))
end
function dotacf(x::IrregularTimeseries, args...; kwargs...)
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
    return Timeseries(msdist(parent(x), lags), lags .* samplingperiod(x))
end
function msdist(x::MultivariateRegular,
                lags = range(0, size(x, 1) - 1, step = 1))
    d = dims(x)[2:end]
    m = mapslices(x -> msdist(x, lags), parent(x); dims = 1)
    return Timeseries(m, lags .* samplingperiod(x), d...)
end

end
