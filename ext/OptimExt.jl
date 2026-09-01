module OptimExt
using Optim
using ForwardDiff
using DimensionalData
using StatsAPI
using StatsBase
using ComponentArrays
import TimeseriesTools: mapple, fit_mapple, MAPPLE, UnivariateSpectrum, Log10𝑓,
    frequency_check, mapple_sort, _safelog10, logweights, _held
using LinearAlgebra

# `box_peaks` controls the peak position/width box. When `true`, each peak is confined near its
# detected centre and to a sub-decade width; this stops the joint refine sliding several peaks
# together and ballooning one into a background-like blob (the "small peaks vanish next to large
# ones" failure on dense fields). When `false`, peaks are free over the whole band — needed when
# the box would trap the optimizer (e.g. a steep multi-segment background whose slopes can only be
# fixed by letting peaks roam transiently). `fit_mapple` refines under both and keeps the better.
#
# These bound EVERY parameter; the ones being held (see `_held`) are dropped from the optimisation
# by `fit_mapple` and their bounds simply go unused.
function mapple_bounds(log_f, log_s, initial_params; box_peaks = true)
    lower = deepcopy(initial_params)
    upper = deepcopy(initial_params)

    # `log_A` is the background amplitude AT f = 1 (`10^log_A`). A data-magnitude ceiling assumes
    # f = 1 sits near the data's largest values, which only holds for DECREASING (spectral) laws.
    # For an INCREASING law (e.g. a structure function / MAD) whose band lies below f = 1, the
    # f = 1 intercept is an upward extrapolation far above the data, so a data-magnitude cap clamps
    # `log_A` and forces the slope shallow. Leave the background uncapped (the loss pins it); keep
    # the peak amplitudes capped at `ampcap` below, where a finite ceiling still guards runaway.
    losS, hiS = extrema(log_s)
    ampcap = hiS + max(hiS - losS, oneunit(hiS))   # strictly above first(log_s) ≤ hiS
    lower.log_A = -Inf
    upper.log_A = Inf

    lower.transition_width = minimum(diff(log_f)) / 4
    upper.transition_width = (maximum(log_f) - minimum(log_f)) / 3

    for i in eachindex(lower.peaks)
        lower.peaks[i].log_A = lower.log_A
        upper.peaks[i].log_A = ampcap

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
        # Keep each breakpoint strictly INSIDE the data band. The old ±2·df window let a knot sit up to two
        # decades OUTSIDE the data, where the segment it bounds collapses to a one-point sliver with an
        # arbitrary slope --- the dominant broken-power-law failure (a MAD rise returning a > 1 off a 1 ms
        # sliver). The `df/6` margin above the low edge guarantees a minimum leading-segment width, so the
        # first component is a real power law rather than an edge sliver.
        lower.components[i].log_f_stop = minimum(log_f) + df / 6
        upper.components[i].log_f_stop = maximum(log_f)

        lower.components[i].β = -Inf
        upper.components[i].β = Inf
    end

    return lower, upper
end

function mapple_loss(params; f, log_s, log_f = log10.(f), w = nothing)
    pred = mapple(f, params; log_f)
    # Fuse log10 + residual + sum-of-squares into one pass over `pred`, avoiding the intermediate
    # `map(log10, …)` array and the `log_s .- pred_log` broadcast temporary on every objective
    # evaluation (each is a full-length allocation, ×hundreds of evals ×ForwardDiff Duals).
    # The weighted branch is a separate loop so the (default) unweighted path keeps its tight one.
    acc = zero(eltype(pred))
    if w === nothing
        @inbounds for i in eachindex(pred, log_s)
            d = log_s[i] - _safelog10(pred[i])
            acc += d * d
        end
    else
        @inbounds for i in eachindex(pred, log_s, w)
            d = log_s[i] - _safelog10(pred[i])
            acc += w[i] * d * d
        end
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

`w` weights the squared log-residuals: `nothing` (the default) fits unweighted, `true` computes
[`logweights`](@ref) from `log_f`, and a vector supplies per-sample weights directly. Weighting
makes the objective minimise `∫ r² d(log f)` rather than `Σ r²`, which matters whenever the fitted
axis is not uniform in log space. Only the scale-free ratio of the weights affects the minimiser,
so normalisation is irrelevant; every refine candidate is scored with the same `w`, leaving the
multistart comparison valid.

`fix` holds arbitrary parameters at given values instead of fitting them: an iterable of
`label => value` pairs, where each label addresses one flat parameter exactly as printed by
`ComponentArrays.labels(m.params)` and each value is in that parameter's own internal units
(log10 units for knots and amplitudes). The initialisation is sorted, so `components[i]` always
refers to ascending-`log_f_stop` order. Recipes: hold an inner knot known a priori so every fitted
curve spans the same segments (`"components[1].log_f_stop" => log10(knee)`); pin a slope known on
theoretical grounds (`"components[1].β" => 0.0` for a Fano factor's flat shot-noise shoulder).
Fixing a knot forbids exactly what a multi-segment fit is usually for --- recovering that knot ---
so pin only quantities that are wanted as constraints, not measurements.

Held parameters are ELIMINATED from the optimisation, not boxed into a narrow interval: they are
removed from the vector Optim searches and spliced back in to evaluate the objective. Boxing leaves
them differentiated on every evaluation (ForwardDiff's cost scales with the free-parameter count),
carried through every line search, and starting a hair from a barrier wall whose gradient distorts
both the convergence test and Fminbox's automatic `μ0` for every other parameter. The last
component's `log_f_stop` is held automatically (see [`_held`](@ref)); an explicit `fix` on it wins.

Extra `kwargs` go to `Optim.Options`.
"""
function fit_mapple(
        log_f, log_s, initial_params;
        algorithm = LBFGS(), autodiff = Optim.ADTypes.AutoForwardDiff(),
        multistart = 0, w = nothing, fix = nothing,
        kwargs...
    ) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)
    initial_params = mapple_sort(initial_params) # so `fix` component indices mean ascending knots
    # `w` has to be an explicit keyword rather than left in `kwargs`: the leftovers are splatted
    # into `Optim.Options` below, which would reject it.
    w === true && (w = logweights(log_f))
    objective = mapple_loss(; f, log_s, log_f, w)

    # Held parameters are ELIMINATED from the optimisation rather than boxed into a sliver. A
    # slivered parameter is still differentiated on every objective evaluation (ForwardDiff's cost
    # scales with the free-parameter count), still carried through every line search, and starts a
    # hair from a barrier wall whose enormous gradient both pollutes Fminbox's convergence test and
    # skews its automatic `μ0` for every OTHER parameter. Measured across 2-4 component fits,
    # eliminating rather than boxing is 1.1-2.0x cheaper per objective call and converges to a
    # markedly better optimum (34x lower loss on a 3-component, 2-peak spectrum).
    held = _held(log_f, initial_params, fix)
    freeidx = setdiff(eachindex(initial_params), keys(held))
    ax = getaxes(initial_params)
    base = collect(initial_params)
    for (i, v) in held
        base[i] = v
    end
    isempty(freeidx) && return ComponentArray(base, ax) # everything held: nothing to optimise

    # Rebuild a full parameter vector from the free coordinates. Generic in the element type so
    # ForwardDiff Duals flow through it.
    function expand(x::AbstractVector{T}) where {T}
        y = Vector{T}(undef, length(base))
        copyto!(y, base)
        @inbounds for (k, i) in enumerate(freeidx)
            y[i] = x[k]
        end
        return ComponentArray(y, ax)
    end
    reduced(x) = objective(expand(x))

    # Bounds, restricted to the coordinates actually being optimised.
    reduce_bounds((lower, upper)) = (collect(lower)[freeidx], collect(upper)[freeidx])
    boxed = reduce_bounds(mapple_bounds(log_f, log_s, initial_params; box_peaks = true))
    free = reduce_bounds(mapple_bounds(log_f, log_s, initial_params; box_peaks = false))

    # Strictly clamp the initialisation inside a box before refining: a rough-init parameter that
    # lands on a recomputed bound makes Fminbox's log-barrier infinite. This generalises the
    # absolute-log_A fix to every parameter (transition_width, log_σ, log_f_stop, …).
    clamp_into((lower, upper)) =
        [_strictclamp(base[i], lower[k], upper[k]) for (k, i) in enumerate(freeidx)]

    # Bound the refine by default so a hard / ill-conditioned spectrum (e.g. a low-dynamic-range curve
    # where the two components are near-degenerate and the objective nearly flat) cannot iterate without
    # limit — the dominant cost when fitting many spectra. Caller `kwargs` override (e.g. a high-accuracy
    # `iterations = 1000` refit, or `time_limit`).
    opts = Optim.Options(; merge((; iterations = 500, outer_iterations = 5), (; kwargs...))...)
    refine(x0, (lower, upper)) = Optim.minimizer(
        optimize(reduced, lower, upper, copy(x0), Fminbox(algorithm), opts; autodiff)
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
    best_loss = reduced(best)
    function consider(cand)
        cand === nothing && return
        l = reduced(cand)
        isfinite(l) && l < best_loss && ((best, best_loss) = (cand, l))
        return
    end

    consider(tryrefine(clamp_into(boxed), boxed))   # boxed peaks: no runaway on dense fields
    consider(tryrefine(clamp_into(free), free))     # free peaks: escapes a box that would trap

    x0free = clamp_into(free)
    for _ in 1:multistart
        consider(tryrefine(_perturb(x0free, free...), free))
    end
    return expand(best)
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
frequency lookup, linear spectral density). `kwargs` are forwarded to [`fit_mapple`](@ref); in
particular `w = true` weights the fit by [`logweights`](@ref) of the spectrum's own axis, so a
grid that is not uniform in log space cannot bias the fit toward its log-denser end, and `fix`
holds arbitrary parameters (addressed by their `ComponentArrays.labels` string) at given values
instead of fitting them.
"""
function StatsAPI.fit!(m::MAPPLE, spectrum::AbstractDimVector; kwargs...)
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(_safelog10, parent(spectrum))
    frequency_check(lookup(spectrum, 1), log_f)
    params = fit_mapple(log_f, log_s, m.params; kwargs...)
    m.params .= params
    return sort!(m)
end


# --- Uncertainty ------------------------------------------------------------------------------
# Laplace (observed-information) uncertainties: invert the Hessian of the residual sum of squares
# at the fitted optimum. One Hessian, no sampler. Only the parameters the fit SEARCHED get an
# uncertainty --- held ones (`fix`, and the inert outermost knot) are constants, and including them
# would make the Hessian singular.
function _laplace(m::MAPPLE, spectrum::AbstractDimVector; fix = nothing, w = nothing)
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(_safelog10, parent(spectrum))
    f = map(exp10, log_f)
    w === true && (w = logweights(log_f))
    p̂ = m.params

    # Only WHICH parameters were held matters here; their values are read from the fitted model, so
    # passing a `fix` with different values cannot move the model off its optimum.
    held = _held(log_f, p̂, fix)
    freeidx = setdiff(eachindex(p̂), keys(held))
    isempty(freeidx) && throw(ArgumentError("every parameter is held; nothing to report on"))
    base = collect(p̂)
    ax = getaxes(p̂)
    function expand(x::AbstractVector{T}) where {T}
        y = Vector{T}(undef, length(base))
        copyto!(y, base)
        @inbounds for (k, i) in enumerate(freeidx)
            y[i] = x[k]
        end
        return ComponentArray(y, ax)
    end
    n = length(log_s)
    # Put a weighted fit on the same footing as an unweighted one, exactly as `_bic` does, so σ²
    # below is a residual variance either way.
    scale = w === nothing ? 1 : n
    obj(x) = scale * mapple_loss(expand(x); f, log_s, log_f, w)

    x̂ = base[freeidx]
    p = length(x̂)
    p ≥ n && throw(ArgumentError("$p free parameters for $n samples: no residual degrees of freedom"))
    H = ForwardDiff.hessian(obj, x̂)
    σ² = obj(x̂) / (n - p)
    # Gaussian residuals: -2 log L = RSS/σ², so the observed information is H/(2σ²).
    # `pinv` rather than `inv` so a (near-)singular direction degrades instead of throwing.
    Σ = 2σ² .* pinv(H)
    κ = cond(H)
    κ > 1.0e10 && @warn "MAPPLE Hessian is near-singular (cond ≈ $(round(κ; sigdigits = 3))); \
        some parameters are not identified by this curve and their uncertainties are unreliable."
    return Σ, freeidx, ax, length(base)
end

"""
    vcov(m::MAPPLE, spectrum; fix = nothing, w = nothing)

Covariance matrix of the fitted parameters, from the Laplace (observed-information) approximation:
`2σ̂² H⁻¹`, where `H` is the Hessian of the residual sum of squares at the optimum. Rows and columns
are the parameters the fit searched, ordered as [`freelabels`](@ref); pass the same `fix` and `w`
the fit used, so the held set and the objective match (only which parameters were held matters --- 
their values come from the fitted model).

!!! warning "Independent residuals are assumed"
    These are exact only for independent, equal-variance residuals. On an AVERAGED curve --- a
    median Fano curve or spectrum --- neighbouring samples share most of their underlying data and
    the residuals are strongly autocorrelated: lag-1 correlations of 0.58 and 0.86 have been
    measured on real curves, understating variances by 3.7x and 13x respectively. Check
    `cor(r[1:(end - 1)], r[2:end])` on [`mapple_residuals`](@ref) and inflate accordingly, or
    prefer resampling the independent units (sessions, trials) when you have them.

!!! note "This is not the uncertainty of an aggregate"
    Fitting a curve that is already an average answers "how well does THIS curve determine the
    parameters", not "how much would they move on a repeat experiment". The latter is usually far
    larger and is best estimated by resampling the units that were averaged.
"""
function StatsAPI.vcov(m::MAPPLE, spectrum::AbstractDimVector; kwargs...)
    return first(_laplace(m, spectrum; kwargs...))
end

"""
    stderror(m::MAPPLE, spectrum; fix = nothing, w = nothing)

Asymptotic standard errors of the fitted parameters, shaped exactly like `m.params` so they can be
read with the same accessors (`betas(stderror(...))`, `.components[2].β`, …). Held parameters are
constants and report `0`. Square roots of the [`vcov`](@ref) diagonal; the same independence
warning applies.
"""
function StatsAPI.stderror(m::MAPPLE, spectrum::AbstractDimVector; kwargs...)
    Σ, freeidx, ax, n = _laplace(m, spectrum; kwargs...)
    se = zeros(eltype(m.params), n)
    se[freeidx] .= sqrt.(abs.(diag(Σ)))
    return ComponentArray(se, ax)
end

end
