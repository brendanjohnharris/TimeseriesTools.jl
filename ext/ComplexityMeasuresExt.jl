module ComplexityMeasuresExt
using TimeseriesTools
using ComplexityMeasures

import ComplexityMeasures: allprobabilities_and_outcomes, SSSet

function allprobabilities_and_outcomes(x::AbstractDimArray, args...; kwargs...)
    return allprobabilities_and_outcomes(TimeseriesTools.decompose(x)..., args...;
                                         kwargs...)
end

end # module
