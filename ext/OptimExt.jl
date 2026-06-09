module OptimExt
using Optim
using ForwardDiff
using DimensionalData
using StatsAPI
using StatsBase
using ComponentArrays
import TimeseriesTools: mapple, fit_mapple, MAPPLE, UnivariateSpectrum, Log10𝑓,
    frequency_check, mapple_sort

# `box_peaks` controls the peak position/width box. When `true`, each peak is confined near its
# detected centre and to a sub-decade width; this stops the joint refine sliding several peaks
# together and ballooning one into a background-like blob (the "small peaks vanish next to large
# ones" failure on dense fields). When `false`, peaks are free over the whole band — needed when
# the box would trap the optimizer (e.g. a steep multi-segment background whose slopes can only be
# fixed by letting peaks roam transiently). `fit_mapple` refines under both and keeps the better.
function mapple_bounds(log_f, log_s, initial_params; box_peaks = true)
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

        if box_peaks
            f0 = initial_params.peaks[i].log_f
            σ0 = initial_params.peaks[i].log_σ
            win = max(3 * σ0, oftype(σ0, 0.2))
            lower.peaks[i].log_f = max(minimum(log_f), f0 - win)
            upper.peaks[i].log_f = min(maximum(log_f), f0 + win)
            lower.peaks[i].log_σ = lower.transition_width / 2
            upper.peaks[i].log_σ = min(upper.transition_width, oftype(σ0, 0.5))
        else
            lower.peaks[i].log_f = minimum(log_f)
            upper.peaks[i].log_f = maximum(log_f)
            lower.peaks[i].log_σ = lower.transition_width / 2
            upper.peaks[i].log_σ = upper.transition_width
        end
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
`log_s` with Optim. Each fit is refined twice — once with peaks boxed near their detected
centre/width and once with peaks free over the band — and the lower-loss result is kept; the
boxed refine prevents peaks running away on dense fields while the free refine wins where the box
would trap the optimizer. The supplied (clamped) initialisation is also a candidate, so the result
never worsens it. With `multistart > 0`, additionally fit from that many randomly perturbed restarts.
Extra `kwargs` go to `Optim.Options`.
"""
function fit_mapple(
        log_f, log_s, initial_params;
        algorithm = LBFGS(), autodiff = Optim.ADTypes.AutoForwardDiff(),
        multistart = 0, kwargs...
    ) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)
    objective = mapple_loss(; f, log_s, log_f)

    boxed = mapple_bounds(log_f, log_s, initial_params; box_peaks = true)
    free = mapple_bounds(log_f, log_s, initial_params; box_peaks = false)

    # Strictly clamp the initialisation inside a box before refining: a rough-init parameter that
    # lands on a recomputed bound makes Fminbox's log-barrier infinite. This generalises the
    # absolute-log_A fix to every parameter (transition_width, log_σ, log_f_stop, …).
    function clamp_into((lower, upper))
        x = deepcopy(initial_params)
        for i in eachindex(x)
            x[i] = _strictclamp(x[i], lower[i], upper[i])
        end
        return x
    end

    refine(x0, (lower, upper)) = Optim.minimizer(
        optimize(
            objective, lower, upper, deepcopy(x0), Fminbox(algorithm),
            Optim.Options(; kwargs...); autodiff
        )
    )
    # Best-effort: a (boxed/perturbed) refine can land where Fminbox cannot build a finite barrier
    # (e.g. a vanishing gradient on a near-perfect fit). Skip such attempts rather than failing the
    # whole fit; another candidate (in the worst case the clamped init) is always retained.
    tryrefine(x0, bnds) =
        try
            refine(x0, bnds)
        catch err
            err isa InterruptException && rethrow()
            nothing
        end

    # Baseline candidate: the clamped init, so the result never worsens the supplied fit.
    best = clamp_into(free)
    best_loss = objective(best)
    function consider(cand)
        cand === nothing && return
        l = objective(cand)
        isfinite(l) && l < best_loss && ((best, best_loss) = (cand, l))
        return
    end

    consider(tryrefine(clamp_into(boxed), boxed))   # boxed peaks: no runaway on dense fields
    consider(tryrefine(clamp_into(free), free))     # free peaks: escapes a box that would trap

    x0free = clamp_into(free)
    for _ in 1:multistart
        consider(tryrefine(_perturb(x0free, free...), free))
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
    fit!(m::MAPPLE, spectrum; kwargs...)
Refine the parameters of a MAPPLE model `m` in place to fit the provided `spectrum` (linear
frequency lookup, linear spectral density). `kwargs` are forwarded to [`fit_mapple`](@ref).
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
