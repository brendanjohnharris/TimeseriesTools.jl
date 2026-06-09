using Unitful
using FFTW
import Unitful.unit

function FFTW.rfft(x::UnitfulTimeseries{<:Quantity{T}}) where {T} # 🐕
    # In this case the rfft is treated as an approximation to the continuous Fourier transform
    t = samplingperiod(x)
    a = x |> eltype |> unit
    x_re = reinterpret(T, parent(x))
    ℱ = rfft(x_re)
    # CTFT has units of amplitude*time. Normalise the DFT to have bins the width of the sampling period.
    return ℱ * (a * t)
end
