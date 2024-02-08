module GeneralizedPhaseExt
using GeneralizedPhase
using TimeseriesTools

import GeneralizedPhase: _generalized_phase, generalized_phase

function _generalized_phase(X::RegularTimeSeries, lp = 0.0)
    mapslices(x -> _generalized_phase(x, samplingrate(X), lp), X; dims = Ti)
end

function generalized_phase(X::RegularTimeSeries, lp = 0.0)
    mapslices(x -> generalized_phase(x, samplingrate(X), lp), X; dims = Ti)
end

end
