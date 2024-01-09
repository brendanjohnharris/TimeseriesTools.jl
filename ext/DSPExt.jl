# module DSPExt

using TimeseriesTools
using ..Unitful
import TimeseriesTools: bandpass, phasestitch
import ..DSP
import ..DSP: hilbert, Bandpass, digitalfilter, filtfilt, unwrap!

hilbert(X::AbstractTimeSeries) = set(X, hilbert(X.data))

analyticphase(x) = x |> hilbert .|> angle
analyticamplitude(x) = x |> hilbert .|> abs
function instantaneousfreq(x)
    y = analyticphase(x)
    unwrap!(y)
    centralderiv!(y)
    return y ./ (2π)
end
function bandpass(x::AbstractArray, fs::A,
                  pass::AbstractVector{B};
                  designmethod = DSP.Butterworth(4)) where {A <: Real, B <: Real}
    DSP.filtfilt(digitalfilter(DSP.Bandpass(pass...; fs), designmethod), x)
end
function bandpass(x::AbstractArray, fs::A,
                  pass::AbstractVector{B};
                  designmethod = DSP.Butterworth(4)) where {A <: Quantity, B <: Quantity}
    DSP.filtfilt(digitalfilter(DSP.Bandpass(ustrip.(pass)...; fs = ustrip(fs)),
                               designmethod), ustrip.(x)) * unit(eltype(x))
end

function bandpass(x::AbstractTimeSeries, fs::A,
                  pass::AbstractVector{B};
                  kwargs...) where {A <: Quantity, B <: Quantity}
    set(x, bandpass(x.data, fs, pass; kwargs...))
end
function bandpass(x::AbstractTimeSeries, fs::A,
                  pass::AbstractVector{B}; kwargs...) where {A <: Real, B <: Real}
    set(x, bandpass(x.data, fs, pass; kwargs...))
end

function bandpass(x::RegularTimeSeries, pass::AbstractVector; kwargs...)
    bandpass(x, samplingrate(x), pass; kwargs...)
end

bandpass(x, pass::NTuple{2}) = bandpass(x, collect(pass))
bandpass(x, fs, pass::NTuple{2}) = bandpass(x, fs, collect(pass))

TimeseriesTools.isoamplitude(x::AbstractVector) = sin.(hilbert(x) .|> angle)
function TimeseriesTools.isoamplitude(x::AbstractArray; dims = 1)
    mapslices(TimeseriesTools.isoamplitude, x; dims)
end

phasewrap(ϕ::Number) = mod(ϕ + π, 2π) - π

function _phasestitch(a::Tuple, b::Tuple; tol = 0.05) # a = (LFP1, PHI1)
    x, xp = a
    y, yp = b

    # ! Remove the half a period at the interface to account for hilbert edge effects
    # c = findlast(xp .< phasewrap(xp[end] - π))
    # x = x[1:c]
    # xp = xp[1:c]
    # c = findfirst(yp .> phasewrap(yp[end] - π))
    # y = y[c:end]
    # yp = yp[c:end]

    # ! Remove one tenth of the samples at the interface to account for hilbert edge effects. Rough, not great
    c = floor(Int, length(xp) / 10)
    x = x[1:(end - c)]
    xp = xp[1:(end - c)]
    c = floor(Int, length(yp) / 10)
    y = y[c:end]
    yp = yp[c:end]

    idxs = 0 .< (yp .- xp[end]) .< tol # Phase is close in value
    i = findfirst(idxs)
    if isnothing(i)
        _, ix = findmin(abs.(yp .- xp[end]))
        ss = (yp .- xp[end])
        @warn "No matching phases found. b is of length $(length(b[1])). The final phase of a is $(xp[end]). The extrema of phases in b is $(extrema(yp)). The smallest difference in phase is $(ss[ix]), at index $ix."
        return a
    else
        x = stitch(x, y[i:end])
        xp = stitch(xp, yp[i:end])
    end
    return x
end

function phasestitch(a::UnivariateTimeSeries, b::UnivariateTimeSeries; kwargs...)
    pha = hilbert(a) .|> angle
    phb = hilbert(b) .|> angle
    return _phasestitch((a, pha), (b, phb); kwargs...)
end

"""
    phasestitch(a::UnivariateTimeSeries, b::UnivariateTimeSeries, [pass]; kwargs...)

Perform phase stitching on two univariate time series `a` and `b` using a specified frequency passband `pass`.
The function applies a bandpass filter to both time series, computes the phase using the Hilbert transform, and then stitches the phases together.

# Arguments
- `a::UnivariateTimeSeries`: The first univariate time series.
- `b::UnivariateTimeSeries`: The second univariate time series.
- `pass`: The frequency passband for the bandpass filter.

# Keyword Arguments
- `kwargs...`: Additional keyword arguments to be passed to the `_phasestitch` function.

# Returns
- A tuple containing the stitched time series and the stitched phase.

"""
function phasestitch(a::UnivariateTimeSeries, b::UnivariateTimeSeries, pass; kwargs...)
    _a = bandpass(a, pass)
    _b = bandpass(b, pass)
    pha = hilbert(_a) .|> angle
    phb = hilbert(_b) .|> angle
    return _phasestitch((a, pha), (b, phb); kwargs...)
end

"""
    phasestitch(X::Union{Tuple, AbstractVector}, [P]; tol = 0.05)

The `phasestitch` function stitches together multiple univariate time series by matching their phases, discarding 1/th of initial and final samples to account for edge effecs.

## Arguments
- `X`: A tuple or vector of univariate time series.
- `P`: Optional. The phase information of the time series. If not provided, it will be calculated using the Hilbert transform.
- `tol`: Optional. The tolerance for matching phases. Default is 0.05.

## Returns
- A single univariate time series obtained by stitching together the input time series.
"""
function TimeseriesTools.phasestitch(X::Union{Tuple{<:UnivariateTimeSeries},
                                              AbstractVector{<:UnivariateTimeSeries}},
                                     P = [hilbert(x) .|> angle for x in X]; tol = 0.05)
    _a = deepcopy(X)
    _ap = deepcopy(P)
    a = []
    ap = []
    aa = []

    # ! Remove one tenth of the samples at the interface to account for hilbert edge effects. Rough, not great
    for i in eachindex(_a)
        c = floor(Int, length(_a[i]) / 10)
        push!(a, _a[i][c:(end - c)])
        push!(ap, _ap[i][c:(end - c)])
    end

    # Now match phases, ready for stitching
    for i in collect(eachindex(a))[2:end]
        x = a[i - 1]
        y = a[i]
        xp = ap[i - 1]
        yp = ap[i]

        idxs = -tol .< (yp .- xp[end]) .< tol
        idx = findfirst(idxs)
        if isnothing(idx)
            _, ix = findmin(abs.(yp .- xp[end]))
            ss = (yp .- xp[end])
            @warn "No matching phases found. y is of length $(length(y)). The final phase of x is $(xp[end]). The extrema of phases in y is $(extrema(yp)). The smallest difference in phase is $(ss[ix]), at index $ix."
        else
            push!(aa, y[idx:end])
        end
    end
    reduce(stitch, aa)
end

function phasestitch(X::Union{Tuple{<:UnivariateTimeSeries},
                              AbstractVector{<:UnivariateTimeSeries}},
                     pass::Union{NTuple{2}, AbstractVector{<:Number}}; kwargs...)
    a = bandpass.(X, [pass])
    P = [hilbert(x) .|> angle for x in a] # Bandpass for phases only
    phasestitch(X, P; kwargs...)
end

# end # module
