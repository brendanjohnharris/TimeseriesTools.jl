export MAPPLE, mapple, fit_mapple
"""
    MAPPLE(params::ComponentArray)
An MAPPLE (Adaptive Peaks and Power-Law Exponents) model for fitting power spectra.
`params` consists of:
- `log_A`: Base log-10 amplitude of the spectrum.
- `components`: An array of components, each with:
    - `log_f_stop`: Log-10 frequency where the component crossfades into the next. The LAST
      component has no next and so is never windowed; its `log_f_stop` bounds nothing, and
      fitting holds it at the top of the fitted band so that it reports where the fit ends.
      `n` components therefore have `n` entries but only `n - 1` crossover points.
    - `β`: Power-law exponent for the component.
- `peaks`: An array of Gaussian peaks, each with:
    - `log_f`: Log-10 center frequency of the peak.
    - `log_σ`: Width of the peak in log-frequency space.
    - `log_A`: Log-10 amplitude of the peak.
- `transition_width`: Width of the transition between components in log-frequency space.
"""
struct MAPPLE <: StatsAPI.RegressionModel
    params::ComponentArray
end

function Base.show(io::IO, m::MAPPLE)
    # Check if we're in a terminal that supports colors
    use_color = get(io, :color, false)

    # Global Parameters (Blue)
    printstyled(io, "\nGlobal Parameters:\n", color = :blue, bold = true)
    printstyled(io, "  log_A:            ")
    println(io, round(m.params.log_A, digits = 4))
    printstyled(io, "  transition_width: ")
    println(io, round(m.params.transition_width, digits = 4))

    # Components (Red). The Hz columns reuse the reporting accessors so the conversions live in
    # exactly one place (see `breakfrequencies`, `peakcentres`, `peakbandwidths`, `peakheights`).
    printstyled(
        io, "\nComponents ($(length(m.params.components))):\n", color = :red,
        bold = true
    )
    bf = breakfrequencies(m)
    for (i, comp) in enumerate(m.params.components)
        printstyled(io, "  Component $i:\n")
        printstyled(io, "    log_f_stop: ")
        println(io, "$(round(comp.log_f_stop, digits = 4))  (≈ $(round(bf[i], sigdigits = 4)) Hz)")
        printstyled(io, "    β:          ")
        println(io, round(comp.β, digits = 4))
    end

    # Peaks (Green)
    printstyled(io, "\nPeaks ($(length(m.params.peaks))):\n", color = :green, bold = true)
    pc, fwhm, ph = peakcentres(m), peakbandwidths(m), peakheights(m)
    for (i, peak) in enumerate(m.params.peaks)
        printstyled(io, "  Peak $i:\n")
        printstyled(io, "    log_f:  ")
        println(io, "$(round(peak.log_f, digits = 4))  (centre ≈ $(round(pc[i], sigdigits = 4)) Hz)")
        printstyled(io, "    log_σ:  ")
        println(io, "$(round(peak.log_σ, digits = 4))  (FWHM ≈ $(round(fwhm[i], sigdigits = 4)) Hz)")
        printstyled(io, "    log_A:  ")
        println(io, "$(round(peak.log_A, digits = 4))  (height ≈ $(round(ph[i], sigdigits = 4)))")
    end
    return
end

# An empty block with a *concrete* ComponentArray eltype (e.g. `Vector{PeakT}()`) makes
# ComponentArrays' `make_idx` divide element-size by a zero count, throwing DivideError.
# An `Any`-eltype empty avoids that path; coerce empty blocks to one.
_emptylike(x) = isempty(x) ? collect(Any, x) : x

function MAPPLE(; peaks, components, log_A, transition_width)
    isempty(components) &&
        throw(ArgumentError("MAPPLE requires at least one component"))
    peaks = _emptylike(peaks)
    components = _emptylike(components)
    return ComponentArray(; log_A, peaks, components, transition_width) |> MAPPLE
end

function mapple_sort!(params)
    _sortblock!(params, :peaks, :log_f)
    _sortblock!(params, :components, :log_f_stop)
    return params
end

# Reorder the nested `block` of `params` in place, sorting its sub-components by field
# `key`. Each sub-component is snapshotted to a plain vector before any write, so the
# writes never alias the shared flat storage; an in-place `sort!` on the nested block
# does alias it and silently duplicates one sub-component over the others.
function _sortblock!(params, block::Symbol, key::Symbol)
    b = getproperty(params, block)
    n = length(b)
    n ≤ 1 && return params
    perm = sortperm([getproperty(b[i], key) for i in 1:n])
    issorted(perm) && return params
    snapshot = [collect(b[i]) for i in 1:n]
    for (newpos, oldpos) in enumerate(perm)
        b[newpos] .= snapshot[oldpos]
    end
    return params
end
mapple_sort(params) = (params = deepcopy(params); mapple_sort!(params); params)

Base.sort!(m::MAPPLE) = (mapple_sort!(m.params); m)   # return the model, per `sort!` convention
function Base.sort(m::MAPPLE)
    m = deepcopy(m)
    sort!(m)
    return m
end

function StatsAPI.predict(m::MAPPLE, freqs)
    return mapple(freqs, m.params)
end
function StatsAPI.predict(m::MAPPLE, freqs::AbstractDimVector)
    return set(freqs, mapple(lookup(freqs, 1), m.params))
end

# A `log10` that never returns `-Inf`/`NaN`: a zero (or underflowed) spectral value is clamped to a
# tiny positive floor. A single silent bin (`s == 0`) would otherwise make the fit objective `-Inf`,
# which (a) drives the rough regression init to `NaN` and (b) leaves Optim with a non-finite landscape
# it cannot descend. Applied to spectral DENSITIES only — frequencies/lags are positive by construction.
@inline _safelog10(x) = log10(max(float(x), eps(float(typeof(x)))))

function frequency_check(f, log_f)
    length(log_f) ≥ 2 ||
        throw(ArgumentError("MAPPLE needs at least two frequencies (got $(length(log_f)))"))
    minimum(log_f) < maximum(log_f) ||
        throw(ArgumentError("frequencies must span a nonzero range"))
    if first(f) ≤ 0
        @warn "Frequencies should be positive"
    end
    df = diff(log_f)
    return if maximum(df) > minimum(df) * 10
        @warn "Frequencies should be evenly spaced in log-space"
    end
end

export logweights

"""
    logweights(log_x)

Cell weights for samples at log-10 positions `log_x`: each point owns a cell centred on it,
reaching half the gap to each neighbour, normalised to sum to one. The two boundary points have
only one neighbour, so their single gap is mirrored and they receive a full-width cell.

Weighting a least-squares fit by these makes it represent the fitted range evenly, however the
sample grid happens to discretise it, rather than weighting by sample count. This matters when the
grid is not uniform in log space; a lag grid rounded to whole samples is the common case, where
short lags land on consecutive integers (linearly spaced) while long lags remain log-spaced, so an
unweighted fit is dominated by whichever end is log-denser.

Mirroring the boundary gap (rather than giving the ends half-cells, as the trapezoidal rule for
`∫ r² d(log x)` over `[x₁, xₙ]` would) is what makes the weighting an *exact* no-op on a uniformly
log-spaced grid: every weight is then `1/n`. Half-cells miss that by `O(1/n)`, which reaches 4% at
`n = 12`, and would down-weight the outermost samples purely for being outermost, even though
nothing physical distinguishes them. Pass `w = true` to [`fit!`](@ref) to have these computed from
the fitted axis.
"""
function logweights(log_x)
    lx = collect(float.(log_x))
    n = length(lx)
    n == 0 && throw(ArgumentError("need at least one sample to weight"))
    n == 1 && return ones(eltype(lx), 1)
    w = similar(lx)
    @inbounds for i in 1:n
        left = i == 1 ? (lx[2] - lx[1]) / 2 : (lx[i] - lx[i - 1]) / 2
        right = i == n ? (lx[n] - lx[n - 1]) / 2 : (lx[i + 1] - lx[i]) / 2
        w[i] = left + right
    end
    return w ./ sum(w)
end

"""
Flat-index => value for every parameter a fit HOLDS rather than searches: the caller's `fix`, plus
the last component's `log_f_stop`. That knot is never windowed (see [`mapple!`](@ref)), so it bounds
nothing and has exactly zero loss-gradient; left in the optimisation it is still boxed, and
Fminbox's log-barrier drifts it through its box for no gain --- 1.25 decades on a two-component fit
--- and towards its neighbour, where a crossing reorders the components and the model changes
discontinuously. Holding it at the top of the fitted band keeps the accessors a uniform shape (`n`
components, `n` entries) and makes `breakpoints(m)[end]` report where the fit ends. A caller `fix`
on the same parameter wins, so the automatic pin is a default rather than a constraint.

Lives here rather than in `OptimExt` because which parameters are structurally held is a property of
the MODEL, and [`_bic`](@ref) has to know them to count degrees of freedom.
"""
function _held(log_f, params, fix)
    labs = labels(params)
    held = Dict{Int, Float64}()
    idx(lab) = findfirst(==(String(lab)), labs)
    if !isempty(params.components)
        i = argmax(k -> params.components[k].log_f_stop, eachindex(params.components))
        held[idx("components[$i].log_f_stop")] = maximum(log_f)
    end
    isnothing(fix) && return held
    for (lab, v) in fix
        i = idx(lab)
        isnothing(i) && throw(
            ArgumentError(
                "fix: no parameter labelled \"$lab\"; valid labels: $(join(labs, ", "))"
            )
        )
        held[i] = v
    end
    return held
end

# Effective free parameters of a fit. A component's breakpoint `log_f_stop` is a *free knot*: a
# nonlinear parameter that hunts for structure (including noise), so its effective complexity is
# ~2 dof, not 1 — standard BIC would otherwise under-penalise and over-select components. We
# therefore charge 3 dof per component (β + a 2-dof knot); peaks keep their 3 literal params.
# Parameters that are HELD (`fix`, and the inert outermost knot) are never searched, so they are
# not charged: a held knot returns its 2 dof, anything else 1.
function _dof(params, held)
    k = 2 + 3 * length(params.components) + 3 * length(params.peaks)
    labs = labels(params)
    for i in keys(held)
        k -= endswith(labs[i], ".log_f_stop") ? 2 : 1
    end
    return k
end

# Bayesian information criterion for a fitted model, treating the log-spectral residuals as
# Gaussian: BIC = n·log(RSS/n) + k·log(n), with `k` from [`_dof`](@ref). Lower is better; used to
# pick the component count.
# `f` is the linear frequency grid (`exp10.(log_f)`), precomputed once by the caller so the
# component-count sweep does not rebuild it per candidate.
# `w` must match the weights the candidates were REFINED under, or selection compares models on a
# different criterion than the one they were fitted to. Scaling by `n` keeps the weighted RSS on the
# same footing as the unweighted one (uniform weights are `1/n`, so `n·Σ wr² == Σ r²` exactly).
function _bic(params, f, log_s, w = nothing, held = Dict{Int, Float64}())
    n = length(f)
    resid = log_s .- map(_safelog10, mapple(f, params))
    rss = w === nothing ? sum(abs2, resid) : n * sum(w .* abs2.(resid))
    return n * log(rss / n) + _dof(params, held) * log(n)
end

# Refine `params` with Optim when OptimExt is loaded, else return unchanged. Component-count
# selection relies on this: the rough fit gives every component the same slope, so only after
# refinement do different counts produce different fits (without Optim the sweep collapses to a
# single component).
function _refine_if_possible(log_f, log_s, params; kwargs...)
    Base.get_extension(@__MODULE__, :OptimExt) === nothing && return params
    return Base.invokelatest(fit_mapple, log_f, log_s, params; kwargs...)
end

# Pick the component count over `candidates` by minimum BIC, refining each candidate so the
# slopes (and fit quality) actually differ. `kwargs` go to the rough peak-finding init; `refine`
# (a NamedTuple) goes to the Optim refinement so refine-only options never leak into `findpeaks`.
# The subset of `fix` addressing parameters this candidate actually has, recording which labels
# matched anywhere in `used`. One `fix` has to serve every candidate in a `:auto` sweep, but the
# candidates have different parameters: `"components[3].β"` is meaningless when fitting two
# components, and `"peaks[2].log_f"` when only one peak was detected. Skipping an inapplicable
# label is the only sensible reading of "hold this if it exists"; `_select_mapple` still rejects a
# label that matched NO candidate, which is a typo rather than a count mismatch.
function _applicable(fix, params, used)
    isnothing(fix) && return nothing
    labs = labels(params)
    keep = [p for p in fix if String(first(p)) ∈ labs]
    for p in keep
        push!(used, String(first(p)))
    end
    return keep
end

function _select_mapple(log_f, log_s, candidates; refine = (;), kwargs...)
    f = map(exp10, log_f)
    # Score candidates on the criterion they were FITTED under: the refine's own weights, and its
    # held parameters, which are not searched and so must not be charged degrees of freedom.
    w = get(refine, :w, nothing)
    w === true && (w = logweights(log_f))
    fix = get(refine, :fix, nothing)
    used = Set{String}()
    best, best_bic = nothing, Inf
    for nc in candidates
        rough = fit_mapple(log_f, log_s; components = nc, kwargs...)
        fix_nc = _applicable(fix, rough, used)
        params = _refine_if_possible(log_f, log_s, rough; merge(refine, (; fix = fix_nc))...)
        bic = _bic(params, f, log_s, w, _held(log_f, params, fix_nc))
        if bic < best_bic
            best, best_bic = params, bic
        end
    end
    if !isnothing(fix)
        unmatched = [String(first(p)) for p in fix if String(first(p)) ∉ used]
        isempty(unmatched) || throw(
            ArgumentError(
                "fix: no candidate model has a parameter labelled " *
                    join(map(l -> "\"$l\"", unmatched), ", ") *
                    "; a label inapplicable at SOME component counts is skipped there, but one that " *
                    "matches none is a typo. Valid labels for the largest candidate: " *
                    join(labels(best), ", ")
            )
        )
    end
    return best
end

"""
    fit(::Type{MAPPLE}, spectrum::UnivariateSpectrum; components = :auto, max_components = 3, refine = (;), kwargs...)

Fit an MAPPLE model to linear frequencies and linear spectral density using peak-finding and
linear regression.

With `components = :auto` (default) the number of broken-power-law segments is chosen by
minimum BIC over `1:max_components`; when Optim is loaded each candidate is refined (so the
selection compares fully-fitted models and the returned model is already refined). Without Optim
loaded the candidates cannot be refined and the `:auto` sweep collapses to a single component (the
rough fits give every component the same slope and so are indistinguishable); load `Optim` for
meaningful component selection. Pass an `Integer` `components` to fix the count and get the rough
fit only, then call [`fit!`](@ref) to refine. `kwargs` (e.g. `peaks`, `w`, `peak_threshold`) are forwarded to the peak-finding init;
`refine` is a NamedTuple of Optim options (e.g. `refine = (; multistart = 4, iterations = 200)`)
forwarded to the refinement of each `:auto` candidate, kept separate so refine-only options do
not leak into the peak finder.

A `refine.fix` label that does not exist for a given candidate is skipped for that candidate rather
than raising: one `fix` has to serve the whole sweep, but `"components[3].β"` is meaningless when
fitting two components, as is `"peaks[2].log_f"` when one peak was detected. A label that matches NO
candidate is still an error --- that is a typo, not a count mismatch. Count-INDEPENDENT labels
(`"log_A"`, `"transition_width"`, `"components[1].*"`) therefore apply to every candidate, while
count-dependent ones apply only where they exist, which is usually what you want: pinning the
outermost slope, say, means something different at each count.

!!! tip "Log-sample noisy spectra first"
    The fit minimises an *unweighted* sum of squared log-10 residuals, so every frequency sample
    counts equally. A linear-frequency spectrum packs most of its samples into the high-frequency
    decades, where periodogram noise is also largest; those noisy samples then dominate the loss.
    Passing the spectrum through [`logsample`](@ref) before fitting gives roughly equal samples per
    decade and averages down the high-frequency noise, yielding a fit that is markedly less
    sensitive to it. This is the recommended preprocessing for any measured (non-synthetic) spectrum.
"""
function StatsAPI.fit(
        ::Type{MAPPLE}, spectrum::AbstractDimVector;
        components = :auto, max_components = 3, refine = (;), kwargs...
    )
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(_safelog10, parent(spectrum))

    frequency_check(lookup(spectrum, 1), log_f)

    params = if components isa Integer
        fit_mapple(log_f, log_s; components, kwargs...)
    else
        _select_mapple(log_f, log_s, 1:max_components; refine, kwargs...)
    end
    return sort(MAPPLE(params))
end
function StatsAPI.fit(::Type{MAPPLE}, f::AbstractVector, s::AbstractVector; kwargs...)
    spectrum = ToolsArray(s, 𝑓(f))
    return StatsAPI.fit(MAPPLE, spectrum; kwargs...)
end

StatsAPI.params(m::MAPPLE) = m.params
peaks(m::MAPPLE) = m.params.peaks
components(m::MAPPLE) = m.params.components

export betas, breakpoints, peakfreqs, peaksigmas, peakamplitudes

# Extract a field across the sub-components of a nested block as a `Vector{Float64}`.
# Broadcasting `getproperty` over a nested block otherwise yields a `Vector{Any}`; an
# empty block has no field vector, so return an empty `Float64[]`.
_field(block, key) = isempty(block) ? Float64[] : collect(Float64, getproperty(block, key))

"Power-law exponents `β` of the broken-power-law components, in fitting order."
betas(m::MAPPLE) = _field(m.params.components, :β)

"""
Component breakpoints as log-10 frequencies (`log_f_stop`), in fitting order. See
[`breakfrequencies`](@ref) for Hz. Only the first `n - 1` are crossover points of the fitted curve:
the final component is never windowed, so its entry bounds nothing and a fit holds it at the top of
the fitted band, where it reports the end of the fitted range rather than a transition.
"""
breakpoints(m::MAPPLE) = _field(m.params.components, :log_f_stop)

"Peak centres in log-10 frequency space (`log_f`). See [`peakcentres`](@ref) for Hz."
peakfreqs(m::MAPPLE) = _field(m.params.peaks, :log_f)

"Peak widths in log-10 frequency space (`log_σ`). See [`peakbandwidths`](@ref) for the FWHM in Hz."
peaksigmas(m::MAPPLE) = _field(m.params.peaks, :log_σ)

"Peak amplitudes in log-10 space (`log_A`). See [`peakheights`](@ref) for the linear height."
peakamplitudes(m::MAPPLE) = _field(m.params.peaks, :log_A)

export peakcentres, peakbandwidths, peakheights, breakfrequencies

# Physical-unit reporting layer. Internally MAPPLE is parameterised in log-10 space (where the
# fit is well-conditioned); these accessors convert to intuitive units at the reporting
# boundary. "Hz" is in whatever frequency unit the fitted spectrum carried, since a centre is
# just `10^log_f`. The log-space accessors above remain for direct parameter access.

"Peak centre frequencies in Hz (`10^log_f`)."
peakcentres(m::MAPPLE) = exp10.(peakfreqs(m))

"Component breakpoint frequencies in Hz (`10^log_f_stop`)."
breakfrequencies(m::MAPPLE) = exp10.(breakpoints(m))

"""
    peakbandwidths(m)
Peak full-width-at-half-maximum in Hz. The Gaussian width in linear frequency is
`σ = f·tanh(log_σ)`, so `FWHM = 2√(2 ln 2)·σ`.
"""
peakbandwidths(m::MAPPLE) = (2 * sqrt(2 * log(2))) .* (peakcentres(m) .* tanh.(peaksigmas(m)))

"Peak heights as linear power added on top of the background (`10^log_A`)."
peakheights(m::MAPPLE) = exp10.(peakamplitudes(m))

export peakparams, peakparams_hz

"All peak parameters as a tuple of vectors `(log_f, log_σ, log_A)` (log-10 space)."
peakparams(m::MAPPLE) = (peakfreqs(m), peaksigmas(m), peakamplitudes(m))

"All peak parameters in physical units as `(centre_Hz, FWHM_Hz, height)`."
peakparams_hz(m::MAPPLE) = (peakcentres(m), peakbandwidths(m), peakheights(m))

export mapple_residuals, mapple_loss, rsquared, vcov, stderror, freelabels

"""
    freelabels(m::MAPPLE, spectrum; fix = nothing)

Labels of the parameters a fit actually searched, in the row/column order of [`vcov`](@ref). The
complement of [`_held`](@ref): every other parameter was held, so it carries no uncertainty.
"""
function freelabels(m::MAPPLE, spectrum::AbstractDimVector; fix = nothing)
    log_f = map(log10, lookup(spectrum, 1))
    held = _held(log_f, m.params, fix)
    return labels(m.params)[setdiff(eachindex(m.params), keys(held))]
end

"""
    mapple_residuals(m::MAPPLE, spectrum)
Log-10 residuals between model `m` and a measured `spectrum` (linear frequency lookup,
linear spectral density), matching the space in which the fit is performed.
"""
function mapple_residuals(m::MAPPLE, spectrum::AbstractDimVector)
    f = lookup(spectrum, 1)
    log_s = _safelog10.(parent(spectrum))
    return log_s .- _safelog10.(mapple(f, m.params))
end

"""
    mapple_loss(m::MAPPLE, spectrum)
Sum of squared log-10 residuals; the objective `fit!` minimises.
"""
mapple_loss(m::MAPPLE, spectrum::AbstractDimVector) = sum(abs2, mapple_residuals(m, spectrum))

"""
    rsquared(m::MAPPLE, spectrum)
Coefficient of determination of the fit, computed in log-10 space.
"""
function rsquared(m::MAPPLE, spectrum::AbstractDimVector)
    log_s = _safelog10.(parent(spectrum))
    ss_res = sum(abs2, mapple_residuals(m, spectrum))
    ss_tot = sum(abs2, log_s .- (sum(log_s) / length(log_s)))
    iszero(ss_tot) && return NaN   # flat spectrum: R² is undefined rather than ±Inf
    return 1 - ss_res / ss_tot
end

"""
    mapple(f, model::ComponentArray; log_f = log10.(f))
    mapple(f, component_params::ComponentArray, peaks::ComponentArray; log_f = log10.(f))

Evaluate the MAPPLE model on linear frequencies `f`, returning the linear spectral density
(the broken-power-law background plus Gaussian peaks). The first form takes a whole-model
`ComponentArray` as held by [`MAPPLE`](@ref); the second keeps the background
(`component_params`) and `peaks` blocks separate, which the optimiser uses to bound the two
blocks independently. Pass a precomputed `log_f = log10.(f)` to avoid recomputing it when the
model is evaluated many times (e.g. inside an objective).
"""
# Whole-model evaluation. mapple! reads `.components`/`.transition_width`/`.log_A` from the first
# argument and `.peaks` from the second, so passing the full model as both avoids splitting it into
# `model[[:log_A,…]]` / `model[[:peaks]]` sub-ComponentArrays (which allocate on every objective eval).
function mapple(f::AbstractVector, model::ComponentArray; log_f = log10.(f))
    ElType = promote_type(eltype(f), eltype(model))
    s = similar(f, ElType)   # mapple! zero-fills before accumulating
    mapple!(s, f, model, model; log_f)
    return s
end
# Split component/peak blocks. One parametric method covers every (Float, ForwardDiff.Dual)
# combination via promote_type — the previous four type-specialised overloads were redundant.
function mapple(f, component_params::ComponentArray, peaks::ComponentArray; log_f = log10.(f))
    ElType = promote_type(eltype(f), eltype(component_params), eltype(peaks))
    s = similar(f, ElType)   # mapple! zero-fills before accumulating
    mapple!(s, f, component_params, peaks; log_f)
    return s
end

function mapple!(s::AbstractVector{El}, f, component_params, peaks; log_f = log10.(f)) where {El}
    fill!(s, zero(El))

    components = component_params.components
    peaks = peaks.peaks

    width = component_params.transition_width
    A_base = exp10(component_params.log_A)

    # Use type-stable operations without mutation
    n_components = length(components)

    # Get sorted indices instead of sorting the array
    sorted_indices = sortperm(components; by = c -> c.log_f_stop)

    # Pre-calculate component amplitudes for continuity
    component_amplitudes = Vector{El}(undef, n_components)
    component_amplitudes[sorted_indices[1]] = A_base

    for j in 2:n_components
        idx_prev = sorted_indices[j - 1]
        idx_curr = sorted_indices[j]

        # Use the actual log_f_stop value, treating the last one specially
        log_f_stop_prev = components[idx_prev].log_f_stop

        # Avoid Inf by using a large finite value or conditional logic
        f_transition = exp10(log_f_stop_prev)
        β_prev = components[idx_prev].β
        β_curr = components[idx_curr].β

        component_amplitudes[idx_curr] = component_amplitudes[idx_prev] *
            f_transition^(β_prev - β_curr)
    end

    @inbounds @fastmath for i in eachindex(f, s)
        # * Add contribution from each component
        for j in 1:n_components
            idx = sorted_indices[j]
            seg = components[idx]
            A_seg = component_amplitudes[idx]

            # * Calculate smooth window weight.
            # A shoulder crossfades one segment into the next, so it belongs only where there IS a
            # next segment. The OUTER edges of the model have none: below the first component and
            # above the last there is nothing to hand over to, and a shoulder there would taper the
            # model to zero off the ends of its own domain rather than blend two power laws. The
            # closing shoulder on the last component in particular made `last(β)` a nuisance
            # parameter absorbing the fade instead of a slope --- a fit whose top segment ended
            # inside the data was attenuated to half at that knot, and the optimiser compensated
            # with a runaway β, or collapsed the knot onto its neighbour to switch the segment off.
            # Consequently `components[end].log_f_stop` is unused: n components have n-1 knots.
            start_weight = j == 1 ? one(El) :
                (one(El) + tanh((log_f[i] - components[sorted_indices[j - 1]].log_f_stop) / width)) / 2
            stop_weight = j == n_components ? one(El) :
                (one(El) + tanh((seg.log_f_stop - log_f[i]) / width)) / 2

            # * Add weighted contribution
            s[i] += start_weight * stop_weight * A_seg * f[i]^seg.β
        end
    end

    # * Add Gaussian peaks
    for peak in peaks
        f_peak = exp10(peak.log_f)
        A_peak = exp10(peak.log_A)
        # log_σ gives a constant width in log_f space. Floor σ at a tiny positive fraction of the
        # centre so a degenerate log_σ ≤ 0 (only reachable for hand-built models — the fit bounds
        # keep log_σ > 0) cannot divide by zero and produce NaN at the peak centre.
        σ_peak = max(f_peak * tanh(peak.log_σ), f_peak * 1.0e-6)

        @inbounds @fastmath for i in eachindex(f, s)
            Δf = f[i] - f_peak
            s[i] += A_peak * exp(-Δf^2 / (2 * σ_peak^2))
        end
    end
    return
end

# Classic median absolute deviation and its normal-consistent scale estimate. Used as a
# robust spread for the peak-detection threshold and the lower-envelope background trim;
# unlike `std`, neither is inflated by the peaks we are trying to detect.
_mad(x) = median(abs.(x .- median(x)))
_rstd(x) = 1.4826 * _mad(x)

# Keep points that lie on the lower envelope of `resid`, i.e. drop the positive excursions
# that are peaks. Returns a boolean mask; used to fit the background through the troughs
# rather than the (peak-biased) full spectrum.
_lower_envelope_mask(resid; k = 2) = resid .≤ (median(resid) + k * _rstd(resid))

# Shared robust lower-envelope loop. Starting from `refit(all-kept)`, iteratively drop the
# positive excursions (peaks) and refit through the remaining troughs. `refit(mask)` returns the
# fitted log-space trend over the FULL grid given a boolean keep-mask; `minpoints` stops the loop
# before the fit becomes underdetermined. Both background estimators (least squares here, Optim
# in OptimExt) route through this loop so the detrending behaves identically regardless of which
# inner fit is used — only the inner model differs, not the robust-iteration scheme.
function _robust_envelope(log_s, refit; iters = 3, minpoints = 0)
    trend = refit(trues(length(log_s)))
    mask = trues(length(log_s))
    for _ in 1:iters
        newmask = _lower_envelope_mask(log_s .- trend)
        (count(newmask) ≤ minpoints || newmask == mask) && break
        mask = newmask
        trend = refit(mask)
    end
    return trend
end

# Estimate the log-space background trend for peak detection: a continuous piecewise-linear
# (hard-break) power law fitted by least squares over the component breakpoints — the actual
# background shape, so peaks on a steep, curved background stand out on the residual. Fitted
# robustly through `_robust_envelope`. This least-squares estimate is used for DETECTION whether
# or not Optim is loaded (a single linear solve, ~1000× cheaper than a bounded Optim background
# fit, and reproducible across load state); the smooth Optim model is reserved for the final fit.
function _background_trend(log_f, log_s, components)
    breaks = [c.log_f_stop for c in components][1:(end - 1)]
    basis = hcat(ones(length(log_f)), log_f, (max.(log_f .- b, 0.0) for b in breaks)...)
    refit(mask) = basis * (basis[mask, :] \ log_s[mask])
    return _robust_envelope(log_s, refit; iters = 3, minpoints = size(basis, 2) - 1)
end

"""
    fit_mapple(log_f, log_s; components, peaks = :auto, peak_threshold = 5.0, max_n_peaks = 8, peak_width_limits, w, kwargs...)

Rough MAPPLE initialisation from log-10 frequencies/spectral density: fit `components`
broken-power-law segments by regression, then detect peaks on the detrended residual.

Peak count:
- `peaks = :auto` (default) keeps every detection whose prominence clears `peak_threshold`
  robust standard deviations of the residual, up to `max_n_peaks` strongest. The default
  threshold is deliberately conservative — noise produces local-maxima prominences several
  times the point-noise scale, so a low threshold over-detects. Detection is best-effort;
  pass an explicit `peaks::Integer` when the count is known.
- `peaks::Integer` takes the strongest that many detections (no threshold, no padding);
  fewer are returned, with a warning, if detection finds fewer.

`peak_width_limits = (wmin, wmax)` (log-frequency units) rejects implausibly narrow
(sub-resolution) or wide detections before counting. `w` is the peak-finder smoothing
window; pass an explicit `minprom` to override the `:auto` threshold.
"""
function fit_mapple(
        log_f, log_s;
        w = max(1, length(log_f) ÷ 100),
        peaks = :auto,
        components,
        peak_threshold = 5.0,
        max_n_peaks = 8,
        peak_width_limits = nothing,
        minprom = nothing,
        kwargs...
    )
    # Validate before touching the (lazily-defaulted) width limits, whose default reduces over
    # `diff(log_f)` and would otherwise crash cryptically on degenerate input.
    components isa Integer && components ≥ 1 ||
        throw(ArgumentError("`components` must be a positive integer (got $components)"))
    length(log_f) ≥ 2 ||
        throw(ArgumentError("need at least two frequencies to fit (got $(length(log_f)))"))
    minimum(log_f) < maximum(log_f) ||
        throw(ArgumentError("frequencies must span a nonzero log range"))
    isnothing(peak_width_limits) &&
        (peak_width_limits = (minimum(diff(log_f)), (maximum(log_f) - minimum(log_f)) / 2))
    first(peak_width_limits) ≤ last(peak_width_limits) ||
        throw(ArgumentError("peak_width_limits must be (wmin, wmax) with wmin ≤ wmax (got $peak_width_limits)"))

    log_A = first(log_s) # Estimate of amplitude

    β = last([ones(length(log_f)) log_f] \ log_s) # Simple linear regression. Start by guessing all components have the same exponent, and evenly distribute the breaks
    log_f_stop = range(extrema(log_f)..., length = components + 1)[2:end]
    transition_width = (maximum(log_f) - minimum(log_f)) / (20)

    components = map(1:components) do i
        ComponentArray(; log_f_stop = log_f_stop[i], β = β)
    end

    # * Find peaks on the residual after removing the robust least-squares piecewise-linear
    #   background. Peaks that ride a steep, curved background — and so are not local maxima of
    #   the raw spectrum — stand out cleanly on the residual.
    trend = _background_trend(log_f, log_s, components)
    residual = ToolsArray(log_s .- trend, Log10𝑓(log_f))

    # For a fixed count, take the strongest peaks with no threshold (the count is the
    # selection). For `:auto`, gate on a data-driven prominence threshold. A (near-)flat residual
    # — a smooth spectrum the background already explains — has no peaks; skip detection, since
    # findpeaks' width estimate is NaN on a flat signal and would throw an InexactError.
    auto = !(peaks isa Integer)
    resvec = parent(residual)
    spread = maximum(resvec) - minimum(resvec)
    if spread ≤ 1.0e-10 * max(one(spread), maximum(log_s) - minimum(log_s))
        proms, bounds = Float64[], Any[]
    else
        auto && isnothing(minprom) && (minprom = peak_threshold * _rstd(resvec))
        _, proms, bounds = findpeaks(residual, w; minprom, kwargs...)
    end

    # Reject implausibly narrow (sub-resolution spike) or wide detections before counting.
    wmin, wmax = peak_width_limits
    keep = [wmin ≤ (maximum(b) - minimum(b)) ≤ wmax for b in bounds]
    proms, bounds = proms[keep], bounds[keep]

    # Select the count: a fixed `peaks` takes the top-N; `:auto` keeps survivors up to
    # `max_n_peaks`. Either way, never fabricate peaks to hit a target.
    ncap = auto ? min(max_n_peaks, length(proms)) : peaks
    if ncap < length(proms)
        ps = sortperm(proms; rev = true)[1:ncap]
        proms, bounds = proms[ps], bounds[ps]
    elseif !auto && ncap > length(proms)
        @warn "Requested up to $peaks peaks but only $(length(proms)) were detected"
    end

    # Initialise each peak's amplitude from the LOCAL background under it, not the global spectral
    # floor. A Gaussian peak adds linear power `A_peak` on top of the background `B`, so at the
    # centre the data is `B + A_peak`; with the (log-10) prominence measuring the rise above the
    # local background trend `t`, `A_peak = 10^t·(10^prom − 1)`, i.e. `log_A = t + log10(10^prom−1)`.
    # Anchoring to the global floor instead (the old behaviour) starts any peak whose local
    # background sits well above the high-frequency floor orders of magnitude too low; at that
    # amplitude the peak is invisible to the loss (∂loss/∂log_A ≈ 0) and the refine cannot grow it,
    # so only the single tallest peak survives. Using the local trend gives every detected peak a
    # start near its true height — `log_A < log_s` at the centre always holds (since
    # `log10(10^prom−1) < prom`), so the start cannot overshoot the data and hijack β.
    # NB: use fresh names inside this closure — assigning `log_A`/`log_σ` here would leak to the
    # enclosing `log_A = first(log_s)` (a `do` block shares enclosing locals).
    peakparams = map(proms, bounds) do prom, bound
        pf = mean(bound)
        pσ = (maximum(bound) - minimum(bound)) / 2
        pσ ≤ 0 && (pσ = transition_width) # A guess
        bg = trend[argmin(abs.(log_f .- pf))]   # local background (log-10) under the peak
        pA = bg + log10(max(expm1(prom * log(10)), eps()))
        return ComponentArray(; log_f = pf, log_σ = pσ, log_A = pA)
    end

    return ComponentArray(; log_A, peaks = _emptylike(peakparams), components, transition_width)
end
