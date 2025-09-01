module GeneralizedPhaseExt
using GeneralizedPhase
using TimeseriesTools

import GeneralizedPhase: _generalized_phase, generalized_phase

function _generalized_phase(X::RegularTimeseries, lp = 0.0)
    Y = deepcopy(X)
    dt = ustripall(samplingrate(X))
    set(Y, mapslices(x -> _generalized_phase(parent(x), dt, lp), X; dims = 𝑡))
end

function generalized_phase(X::RegularTimeseries, lp = 0.0)
    Y = deepcopy(X)
    dt = ustripall(samplingrate(X))
    set(Y, mapslices(x -> generalized_phase(parent(x), dt, lp), X; dims = 𝑡))
end

end
