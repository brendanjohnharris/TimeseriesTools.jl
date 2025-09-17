module OptimExt
using Optim
using DimensionalData
import TimeseriesTools: oneoneff, fit_oneoneff

function fit_oneoneff(logspectrum::AbstractDimVector, initial_params; kwargs...) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    function objective(params)
        pred = oneoneff(lookup(logspectrum, 1), params)
        return sum((logspectrum .- pred) .^ 2) # May want to choose smarter loss
    end
    result = optimize(objective, initial_params, LBFGS(); kwargs...)
    return Optim.minimizer(result)
end

end
