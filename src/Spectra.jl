using FFTW

export  FreqDim, Freq,
        AbstractSpectrum, RegularSpectrum,
        spectrum,
        FreqIndex, RegularFreqIndex

abstract type FreqDim{T} <: DimensionalData.IndependentDim{T} end
DimensionalData.@dim Freq FreqDim "Freq"
FreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim}
AbstractSpectrum = AbstractDimArray{T, N, <:FreqIndex, B} where {T, N, B}
RegularFreqIndex = Tuple{A, Vararg{DimensionalData.Dimension}} where {A<:FreqDim{<:RegularIndex}}
RegularSpectrum = AbstractDimArray{T, N, <:RegularFreqIndex, B} where {T, N, B}

function spectrum(x::RegularTimeSeries, f_min=samplingrate(x)/min(length(x)÷4, 1000))
    fs = samplingrate(x)
    n = length(x)
    validfreqs = rfftfreq(n, fs)
    if f_min == 0
        f_min = validfreqs[2]
        nfft = floor(Int, n/2)*2
    else
        nfft = ceil(Int, fs/f_min)
    end
    if f_min < validfreqs[2]
        error("Cannot resolve an `f_min` of $f_min")
    end

    isodd(nfft) && (nfft += 1)
    window = nfft÷2
    overlap = window÷2
    hann_window = 0.5 .- 0.5 .* cos.(2 * π * (0:(nfft - 1)) / (nfft - 1))
    n_segments = floor(Int, (n - nfft) / (nfft - overlap) + 1)

    # Get the type of the spectrum
    X = rfft(x[1:nfft]).^2
    Pxx = convertconst.(zeros(nfft ÷ 2 + 1, n_segments), (first(X),))

    for i in 1:n_segments
        start_idx = (i - 1) * (nfft - overlap) + 1
        end_idx = start_idx + nfft - 1
        segment = x[start_idx:end_idx] .* hann_window

        X = rfft(segment) / nfft
        Pxx[:, i] .= (abs.(X).^2) ./ sum(hann_window.^2)
    end

    # Normalize the averaged periodogram
    Pxx[2:end-1, :] .*= 2

    # Calculate the frequencies
    freqs = range(convertconst(0, fs), stop=fs/2, length=size(Pxx, 1))
    df = step(freqs)

    # Normalize the power spectrum to obey Parseval's theorem
    Pxx = Pxx ./ sum(Pxx, dims=1).*(df).*unit(df) # Maybe do proper integration
    Pxx = 0.5 .* Pxx .* sum(x.^2)./(fs)

    return DimArray(Pxx, (Freq(freqs), Dim{:window}(1:n_segments)))
end

function spectrum(X::typeintersect(MultivariateTS, RegularTS))
    mapslices(spectrum, X, dims=Ti)
end
