module SignalDecompositionExt
using TimeseriesTools
using SignalDecomposition
import SignalDecompsition.decompose

function decompose(x::AbstractTimeseries, args...; kwargs...)
    s, n = decompose(TimeseriesTools.decompose(x)..., args...; kwargs...)
    return Timeseries(s, times(x)), Timeseries(n, times(x))
end

end
