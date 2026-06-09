module OptimExt
using Optim
using ForwardDiff
using DimensionalData
using StatsAPI
using StatsBase
using ComponentArrays
import TimeseriesTools: mapple, fit_mapple, MAPPLE, UnivariateSpectrum, Log10𝑓,
    frequency_check, mapple_sort, _lower_envelope_mask

function mapple_bounds(log_f, log_s, initial_params)
    lower = deepcopy(initial_params)
    upper = deepcopy(initial_params)

    lower.log_A = -Inf
    upper.log_A = 100 * (maximum(log_s) - minimum(log_s))

    lower.transition_width = minimum(diff(log_f)) / 4
    upper.transition_width = (maximum(log_f) - minimum(log_f)) / 3

    for i in eachindex(lower.peaks)
        lower.peaks[i].log_A = lower.log_A
        upper.peaks[i].log_A = upper.log_A

        lower.peaks[i].log_f = minimum(log_f)
        upper.peaks[i].log_f = maximum(log_f)

        lower.peaks[i].log_σ = lower.transition_width / 2
        upper.peaks[i].log_σ = upper.transition_width
    end

    for i in eachindex(lower.components)
        df = maximum(log_f) - minimum(log_f)
        lower.components[i].log_f_stop = minimum(log_f) - df * 2
        upper.components[i].log_f_stop = maximum(log_f) + df * 2

        lower.components[i].β = -Inf
        upper.components[i].β = Inf
    end

    return lower, upper
end

function mapple_loss(params; f, log_s, log_f = log10.(f))
    pred_log = map(log10, mapple(f, params; log_f))
    return sum(abs2, log_s .- pred_log)
end
mapple_loss(; kwargs...) = params -> mapple_loss(params; kwargs...)

"""
    fit_mapple(log_f, log_s, initial_params; multistart = 0, algorithm = LBFGS(), kwargs...)

Refine MAPPLE `initial_params` against log-10 frequencies `log_f` and log-10 spectral density
`log_s` with Optim. With `multistart > 0`, additionally fit from that many randomly perturbed
restarts and keep the lowest-loss result — an escape hatch for hard/degenerate spectra that is
usually unnecessary once peaks are well initialised. Extra `kwargs` go to `Optim.Options`.
"""
function fit_mapple(
        log_f, log_s, initial_params;
        algorithm = LBFGS(), autodiff = Optim.ADTypes.AutoForwardDiff(),
        multistart = 0, kwargs...
    ) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)
    lower, upper = mapple_bounds(log_f, log_s, initial_params)
    objective = mapple_loss(; f, log_s, log_f)

    refine(x0) = Optim.minimizer(
        optimize(
            objective, lower, upper, deepcopy(x0), Fminbox(algorithm),
            Optim.Options(; kwargs...); autodiff
        )
    )

    # Refine from the supplied initialisation, then from `multistart` perturbed restarts,
    # keeping the lowest-loss result, so restarts never worsen a good fit. (The detrended
    # peak init already keeps single-descent fits out of the degenerate β / log_A basin on
    # the benchmark, so multistart defaults off.)
    best = refine(initial_params)
    best_loss = objective(best)
    for _ in 1:multistart
        cand = refine(_perturb(initial_params, lower, upper))
        cand_loss = objective(cand)
        if cand_loss < best_loss
            best, best_loss = cand, cand_loss
        end
    end
    return best
end

# A randomly perturbed copy of `params`, clamped to `[lower, upper]`, for multistart.
function _perturb(params, lower, upper)
    x = deepcopy(params)
    for i in eachindex(x)
        x[i] = clamp(x[i] + 0.5 * randn(), lower[i], upper[i])
    end
    return x
end

# Full zero-peak background fit used by `TimeseriesTools._background_trend` to detrend the
# spectrum before peak detection (chosen automatically when Optim is loaded): fit the smooth
# broken power law with no peaks, then return its log-space prediction as the trend.
function optim_background_trend(log_f, log_s, components, transition_width; iters = 2)
    nopeaks = map(_ -> ComponentArray(; log_f = 0.0, log_σ = 0.0, log_A = 0.0), 1:0)
    background = ComponentArray(;
        log_A = first(log_s), peaks = nopeaks,
        components = components, transition_width = transition_width
    )
    f = map(exp10, log_f)

    # Robust fit: peaks only push power up, so refit the zero-peak background through the
    # lower envelope (excluding positive excursions) so oscillations don't lift the trend.
    fitted = background
    mask = trues(length(log_f))
    for _ in 1:iters
        fitted = fit_mapple(log_f[mask], log_s[mask], background)
        resid = log_s .- map(log10, mapple(f, fitted))
        newmask = _lower_envelope_mask(resid)
        (count(newmask) ≤ length(fitted) || newmask == mask) && break
        mask = newmask
    end
    return map(log10, mapple(f, fitted))
end

"""
    fit!(m::MAPPLE, logspectrum; kwargs...)
Refine the parameters of a MAPPLE model `m` to fit the provided `spectrum`.
"""
function StatsAPI.fit!(m::MAPPLE, spectrum::AbstractDimVector; kwargs...)
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(log10, parent(spectrum))
    frequency_check(lookup(spectrum, 1), log_f)
    params = fit_mapple(log_f, log_s, m.params; kwargs...)
    m.params .= params
    return sort!(m)
end

end
