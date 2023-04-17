using FFTW
using Statistics

export  FreqDim, Freq,
        AbstractSpectrum, RegularSpectrum, MultivariateSpectrum,
        spectrum, energyspectrum, powerspectrum,
        _energyspectrum, _powerspectrum,
        FreqIndex, RegularFreqIndex

abstract type FreqDim{T} <: DimensionalData.IndependentDim{T} end
DimensionalData.@dim Freq FreqDim "Freq"
FreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim}
AbstractSpectrum = AbstractDimArray{T, N, <:FreqIndex, B} where {T, N, B}
RegularFreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim{<:RegularIndex}}
RegularSpectrum = AbstractDimArray{T, N, <:RegularFreqIndex, B} where {T, N, B}
MultivariateSpectrum = AbstractSpectrum{T, 2} where T

function _energyspectrum(x::RegularTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000))
    fs = samplingrate(x)
    n = length(x)
    validfreqs = rfftfreq(n, fs)
    if f_min == 0
        f_min = validfreqs[2]
        nfft = floor(Int, n/2)*2
    else
        nfft = ceil(Int, ustrip(fs)/ustrip(f_min))
    end
    if ustrip(f_min) < ustrip(validfreqs[2])
        error("Cannot resolve an `f_min` of $f_min")
    end

    isodd(nfft) && (nfft += 1)
    window = nfft÷2
    overlap = window÷2
    hann_window = 0.5 .- 0.5 .* cos.(2 * π * (0:(nfft - 1)) / (nfft - 1))
    n_segments = floor(Int, (n - nfft) / (nfft - overlap) + 1)

    # Get the type of the spectrum
    u = unit(eltype(x)) * unit(eltype(dims(x, Ti)))
    S̄ = zeros(nfft ÷ 2 + 1, n_segments)*u^2
    for i in 1:n_segments
        start_idx = (i - 1) * (nfft - overlap) + 1
        end_idx = start_idx + nfft - 1
        segment = x[start_idx:end_idx] .* hann_window

        y = rfft(segment) / nfft
        y |> eltype |> unit == NoUnits && (y = y * u)
        S̄[:, i] .= (abs.(y).^2) / sum(hann_window.^2)
    end


    # Calculate the frequencies
    freqs = range(convertconst(0, fs), stop=fs/2, length=size(S̄, 1))
    df = step(freqs)

    # Normalize the mean energy spectrum to obey Parseval's theorem
    meanS̄ = mean(S̄, dims=2)
    S̄ = S̄ ./ ustrip((sum(meanS̄) - 0.5.*meanS̄[1]) .* df) # Subtract the zero frequency component twice, so that it doesn't bias when we divide by a half
    S̄ = 0.5 * S̄ .* ustrip(sum(x.^2) ./ fs) # Normalized to have total energy equal to energy of signal. Ala parseval. 0.5 because we only return the positive half of the spectrum.

    return DimArray(S̄, (Freq(freqs), Dim{:window}(1:n_segments)))
end

_energyspectrum(x::typeintersect(MultivariateTS, RegularTS), args...; kwargs...) = cat([_energyspectrum(_x, args...; kwargs...) for _x in eachslice(x, dims=2)]..., dims=dims(x, 2))

energyspectrum(x::RegularTimeSeries, args...; kwargs...) = dropdims(mean(_energyspectrum(x, args...; kwargs...), dims=Dim{:window}); dims=Dim{:window})


function _powerspectrum(x::AbstractTimeSeries, args...; kwargs...)
    S̄ = _energyspectrum(x, args...; kwargs...)
    return S̄ ./ duration(x)
end

powerspectrum(x::AbstractTimeSeries, args...; kwargs...) = dropdims(mean(_powerspectrum(x, args...; kwargs...), dims=Dim{:window}); dims=Dim{:window})

spectrum = powerspectrum
