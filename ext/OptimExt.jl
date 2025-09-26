module OptimExt
using Optim
using DimensionalData
using StatsAPI
import TimeseriesTools: apple, fit_apple, APPLE, UnivariateSpectrum, Log10𝑓, frequency_check

function fit_apple(log_f, log_s, initial_params;
                   algorithm = LBFGS(), autodiff = :forward, kwargs...) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)
    function objective(params)
        pred = apple(f, params)
        pred = map(log10, pred)
        return sum((log_s .- pred) .^ 2) # May want to choose smarter loss
    end
    result = optimize(objective, initial_params, algorithm, Optim.Options(; kwargs...);
                      autodiff)
    return Optim.minimizer(result)
end

"""
    fit!(m::APPLE, logspectrum; kwargs...)
Refine the parameters of a APPLE model `m` to fit the provided `spectrum`.
"""
function StatsAPI.fit!(m::APPLE, spectrum::AbstractDimVector{T, D};
                       kwargs...) where {T, d, D <: Tuple{<:d}}
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(log10, parent(spectrum))

    frequency_check(lookup(spectrum, 1), log_f)

    params = fit_apple(log_f, log_s, m.params; kwargs...)
    m.params .= params
    sort!(m)
end
function StatsAPI.fit!(m::Type{APPLE}, f::AbstractVector, s::AbstractVector; kwargs...)
    spectrum = ToolsArray(s, f)
    StatsAPI.fit!(m, spectrum; kwargs...)
end

end
