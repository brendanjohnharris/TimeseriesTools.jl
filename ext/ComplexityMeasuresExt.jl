module ComplexityMeasuresExt
using TimeseriesTools
using ComplexityMeasures

import ComplexityMeasures: allprobabilities_and_outcomes

function allprobabilities_and_outcomes(x::AbstractToolsArray, args...; kwargs...)
    return allprobabilities_and_outcomes(TimeseriesTools.decompose(x)..., args...;
                                         kwargs...)
end

function ComplexityMeasures.StateSpaceSets.StateSpaceSet(x::AbstractToolsArray)
    ComplexityMeasures.StateSpaceSets.StateSpaceSet(parent(x))
end

end # module
