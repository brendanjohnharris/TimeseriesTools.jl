module SignalDecompositionExt
using TimeseriesTools
using SignalDecomposition
import SignalDecompsition.decompose

function decompose(x::AbstractTimeSeries, args...; kwargs...)
    s, n = decompose(TimeseriesTools.decompose(x)..., args...; kwargs...)
    return TimeSeries(times(x), s), TimeSeries(times(x), n)
end

end
