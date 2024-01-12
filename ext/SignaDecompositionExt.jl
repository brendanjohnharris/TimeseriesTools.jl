module SignaDecompositionExt
using TimeseriesTools
using SignalDecomposition
import SignaDecompsition.decompose

function decompose(x::AbstractTimeSeries, args...; kwargs...)
    s, n = decompose(TimeseriesTools.decompose(x)..., args...; kwargs...)
    return TimeSeries(times(x), s), TimeSeries(times(x), n)
end

end
