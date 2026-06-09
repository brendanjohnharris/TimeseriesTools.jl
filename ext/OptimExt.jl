module OptimExt
using Optim
using ForwardDiff
using DimensionalData
using StatsAPI
using StatsBase
using ComponentArrays
import TimeseriesTools: mapple, fit_mapple, MAPPLE, UnivariateSpectrum, Log10𝑓,
    frequency_check, mapple_sort

function mapple_bounds(log_f, log_s, initial_params)
    lower = deepcopy(initial_params)
    upper = deepcopy(initial_params)

    # `log_A` (and each peak `log_A`) lives on the same scale as `log_s` — a base amplitude is
    # `10^log_A` — so the ceiling must be ABSOLUTE, not range-relative. The rough init
    # `log_A = first(log_s)` is an absolute level; a range-relative `100*(max-min)` collapses to
    # 0 on a flat spectrum, putting the init outside the box so Fminbox rejects the start.
    losS, hiS = extrema(log_s)
    ampcap = hiS + max(hiS - losS, oneunit(hiS))   # strictly above first(log_s) ≤ hiS
    lower.log_A = -Inf
    upper.log_A = ampcap

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
    pred = mapple(f, params; log_f)
    # Fuse log10 + residual + sum-of-squares into one pass over `pred`, avoiding the intermediate
    # `map(log10, …)` array and the `log_s .- pred_log` broadcast temporary on every objective
    # evaluation (each is a full-length allocation, ×hundreds of evals ×ForwardDiff Duals).
    acc = zero(eltype(pred))
    @inbounds for i in eachindex(pred, log_s)
        d = log_s[i] - log10(pred[i])
        acc += d * d
    end
    return acc
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

    # Strictly clamp the initialisation inside the box before refining: a rough-init parameter
    # that lands on a recomputed bound makes Fminbox's log-barrier infinite. This generalises the
    # absolute-log_A fix to every parameter (transition_width, log_σ, log_f_stop, …).
    x0 = deepcopy(initial_params)
    for i in eachindex(x0)
        x0[i] = _strictclamp(x0[i], lower[i], upper[i])
    end

    # Refine from the supplied initialisation, then from `multistart` perturbed restarts,
    # keeping the lowest-loss result, so restarts never worsen a good fit. (The detrended
    # peak init already keeps single-descent fits out of the degenerate β / log_A basin on
    # the benchmark, so multistart defaults off.)
    best = refine(x0)
    best_loss = objective(best)
    # Multistart is best-effort: a perturbed restart can land where Fminbox cannot build a finite
    # barrier (e.g. a vanishing gradient on a near-perfect fit). Skip such restarts rather than
    # failing the whole fit; the unperturbed result is always retained.
    for _ in 1:multistart
        cand = try
            refine(_perturb(x0, lower, upper))
        catch err
            err isa InterruptException && rethrow()
            continue
        end
        cand_loss = objective(cand)
        if cand_loss < best_loss
            best, best_loss = cand, cand_loss
        end
    end
    return best
end

# Clamp `v` strictly INSIDE `(lo, hi)`, never onto a boundary where Fminbox's log-barrier is
# `Inf`. Either bound may be infinite, in which case that side is left open.
function _strictclamp(v, lo, hi)
    if isfinite(lo) && isfinite(hi)
        ϵ = 1.0e-6 * (hi - lo)
        return clamp(v, lo + ϵ, hi - ϵ)
    elseif isfinite(lo)
        return max(v, lo + 1.0e-6 * max(abs(lo), one(lo)))
    elseif isfinite(hi)
        return min(v, hi - 1.0e-6 * max(abs(hi), one(hi)))
    else
        return v
    end
end

# A randomly perturbed copy of `params`, kept strictly inside `[lower, upper]`, for multistart.
function _perturb(params, lower, upper)
    x = deepcopy(params)
    for i in eachindex(x)
        x[i] = _strictclamp(x[i] + 0.5 * randn(), lower[i], upper[i])
    end
    return x
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
