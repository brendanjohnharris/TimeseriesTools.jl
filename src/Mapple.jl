export MAPPLE, mapple, fit_mapple
"""
    MAPPLE(params::ComponentArray)
An MAPPLE (Adaptive Peaks and Power-Law Exponents) model for fitting power spectra.
`params` consists of:
- `log_A`: Base log-10 amplitude of the spectrum.
- `components`: An array of components, each with:
    - `log_f_stop`: Log-10 frequency where the component transitions to the next.
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

    # Components (Red)
    printstyled(
        io, "\nComponents ($(length(m.params.components))):\n", color = :red,
        bold = true
    )
    for (i, comp) in enumerate(m.params.components)
        printstyled(io, "  Component $i:\n")
        printstyled(io, "    log_f_stop: ")
        println(io, "$(round(comp.log_f_stop, digits = 4))  (≈ $(round(exp10(comp.log_f_stop), sigdigits = 4)) Hz)")
        printstyled(io, "    β:          ")
        println(io, round(comp.β, digits = 4))
    end

    # Peaks (Green)
    printstyled(io, "\nPeaks ($(length(m.params.peaks))):\n", color = :green, bold = true)
    for (i, peak) in enumerate(m.params.peaks)
        fc = exp10(peak.log_f)
        fwhm = 2 * sqrt(2 * log(2)) * fc * tanh(peak.log_σ)
        printstyled(io, "  Peak $i:\n")
        printstyled(io, "    log_f:  ")
        println(io, "$(round(peak.log_f, digits = 4))  (centre ≈ $(round(fc, sigdigits = 4)) Hz)")
        printstyled(io, "    log_σ:  ")
        println(io, "$(round(peak.log_σ, digits = 4))  (FWHM ≈ $(round(fwhm, sigdigits = 4)) Hz)")
        printstyled(io, "    log_A:  ")
        println(io, "$(round(peak.log_A, digits = 4))  (height ≈ $(round(exp10(peak.log_A), sigdigits = 4)))")
    end
    return
end

# An empty block with a *concrete* ComponentArray eltype (e.g. `Vector{PeakT}()`) makes
# ComponentArrays' `make_idx` divide element-size by a zero count, throwing DivideError.
# An `Any`-eltype empty avoids that path; coerce empty blocks to one.
_emptylike(x) = isempty(x) ? collect(Any, x) : x

function MAPPLE(; peaks, components, log_A, transition_width)
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

Base.sort!(m::MAPPLE; kwargs...) = mapple_sort!(m.params)
function Base.sort(m::MAPPLE; kwargs...)
    m = deepcopy(m)
    sort!(m; kwargs...)
    return m
end

function StatsAPI.predict(m::MAPPLE, freqs)
    return mapple(freqs, m.params)
end
function StatsAPI.predict(m::MAPPLE, freqs::AbstractDimVector)
    return set(freqs, mapple(lookup(freqs, 1), m.params))
end

function frequency_check(f, log_f)
    if first(f) ≤ 0
        @warn "Frequencies should be positive"
    end
    df = diff(log_f)
    return if maximum(df) > minimum(df) * 10
        @warn "Frequencies should be evenly spaced in log-space"
    end
end

# Bayesian information criterion for a fitted model, treating the log-spectral residuals as
# Gaussian: BIC = n·log(RSS/n) + k·log(n). Lower is better; used to pick the component count.
# `k` counts effective free parameters. A component's breakpoint `log_f_stop` is a *free knot*:
# a nonlinear parameter that hunts for structure (including noise), so its effective complexity
# is ~2 dof, not 1 — standard BIC would otherwise under-penalise and over-select components. We
# therefore charge 3 dof per component (β + a 2-dof knot); peaks keep their 3 literal params.
function _bic(params, log_f, log_s)
    n = length(log_f)
    rss = sum(abs2, log_s .- map(log10, mapple(map(exp10, log_f), params)))
    k = 2 + 3 * length(params.components) + 3 * length(params.peaks)
    return n * log(rss / n) + k * log(n)
end

# Refine `params` with Optim when OptimExt is loaded, else return unchanged. Component-count
# selection relies on this: the rough fit gives every component the same slope, so only after
# refinement do different counts produce different fits (without Optim the sweep collapses to a
# single component).
function _refine_if_possible(log_f, log_s, params)
    Base.get_extension(@__MODULE__, :OptimExt) === nothing && return params
    return Base.invokelatest(fit_mapple, log_f, log_s, params)
end

# Pick the component count over `candidates` by minimum BIC, refining each candidate so the
# slopes (and fit quality) actually differ.
function _select_mapple(log_f, log_s, candidates; kwargs...)
    best, best_bic = nothing, Inf
    for nc in candidates
        params = _refine_if_possible(log_f, log_s, fit_mapple(log_f, log_s; components = nc, kwargs...))
        bic = _bic(params, log_f, log_s)
        if bic < best_bic
            best, best_bic = params, bic
        end
    end
    return best
end

"""
    fit(::Type{MAPPLE}, spectrum::UnivariateSpectrum; components = :auto, max_components = 3, kwargs...)

Fit an MAPPLE model to linear frequencies and linear spectral density using peak-finding and
linear regression.

With `components = :auto` (default) the number of broken-power-law segments is chosen by
minimum BIC over `1:max_components`; when Optim is loaded each candidate is refined (so the
selection compares fully-fitted models and the returned model is already refined). Pass an
`Integer` `components` to fix the count and get the rough fit only, then call [`fit!`](@ref) to
refine. `kwargs` (e.g. `peaks`, `w`, `peak_threshold`) are forwarded to the peak-finding init.
Please consider using 'logsample'd spectra for a fit that is less sensitive to high-frequency
noise.
"""
function StatsAPI.fit(
        ::Type{MAPPLE}, spectrum::AbstractDimVector;
        components = :auto, max_components = 3, kwargs...
    )
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(log10, parent(spectrum))

    frequency_check(lookup(spectrum, 1), log_f)

    params = if components isa Integer
        fit_mapple(log_f, log_s; components, kwargs...)
    else
        _select_mapple(log_f, log_s, 1:max_components; kwargs...)
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

betas(m::MAPPLE) = _field(m.params.components, :β)
breakpoints(m::MAPPLE) = _field(m.params.components, :log_f_stop)
peakfreqs(m::MAPPLE) = _field(m.params.peaks, :log_f)
peaksigmas(m::MAPPLE) = _field(m.params.peaks, :log_σ)
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

export mapple_residuals, mapple_loss, rsquared

"""
    mapple_residuals(m::MAPPLE, spectrum)
Log-10 residuals between model `m` and a measured `spectrum` (linear frequency lookup,
linear spectral density), matching the space in which the fit is performed.
"""
function mapple_residuals(m::MAPPLE, spectrum::AbstractDimVector)
    f = lookup(spectrum, 1)
    log_s = log10.(parent(spectrum))
    return log_s .- log10.(mapple(f, m.params))
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
    log_s = log10.(parent(spectrum))
    ss_res = sum(abs2, mapple_residuals(m, spectrum))
    ss_tot = sum(abs2, log_s .- (sum(log_s) / length(log_s)))
    iszero(ss_tot) && return NaN   # flat spectrum: R² is undefined rather than ±Inf
    return 1 - ss_res / ss_tot
end

function mapple(f::AbstractVector, model::ComponentArray{T}; log_f = log10.(f)) where {T}
    components = model[[:log_A, :transition_width, :components]]
    peaks = model[[:peaks]]
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, components, peaks; log_f)
    return s
end
function mapple(
        f, component_params::ComponentArray{T},
        peaks::ComponentArray{F}; log_f = log10.(f)
    ) where {T, F <: AbstractFloat}
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, component_params, peaks; log_f)
    return s
end
function mapple(
        f, component_params::ComponentArray{F},
        peaks::ComponentArray{T}; log_f = log10.(f)
    ) where {T, F <: AbstractFloat}
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, component_params, peaks; log_f)
    return s
end
function mapple(
        f, component_params::ComponentArray{F1},
        peaks::ComponentArray{F2}; log_f = log10.(f)
    ) where {F1 <: AbstractFloat, F2 <: AbstractFloat}
    ElType = promote_type(eltype(f), F1, F2)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
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

    # * Evaluate the model
    log_f_min = minimum(log_f) - 5.0

    @inbounds @fastmath for i in eachindex(f, s)
        # * Add contribution from each component
        for j in 1:n_components
            idx = sorted_indices[j]
            seg = components[idx]
            A_seg = component_amplitudes[idx]

            # * Determine component boundaries without Inf
            log_f_start = if j == 1
                log_f_min
            else
                components[sorted_indices[j - 1]].log_f_stop
            end

            log_f_stop = if j == n_components
                seg.log_f_stop # Inf # No transition width for final component
            else
                seg.log_f_stop
            end

            # * Calculate smooth window weight
            if n_components == 1 # No transition width
                weight = one(El)
            else
                start_weight = (one(El) + tanh((log_f[i] - log_f_start) / width)) / 2
                stop_weight = (one(El) + tanh((log_f_stop - log_f[i]) / width)) / 2
                weight = start_weight * stop_weight
            end

            # * Add weighted contribution
            s[i] += weight * A_seg * f[i]^seg.β
        end
    end

    # * Add Gaussian peaks
    for peak in peaks
        f_peak = exp10(peak.log_f)
        A_peak = exp10(peak.log_A)
        σ_peak = f_peak * tanh(peak.log_σ) # log_σ gives a constant width in log_f space

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

# Estimate the log-space background trend for peak detection. Default: a continuous
# piecewise-linear (hard-break) power law fitted by least squares over the component
# breakpoints — the actual background shape, so peaks on a steep, curved background stand out
# on the residual. Fitted robustly: peaks only push power up, so the background is traced by
# iteratively refitting through the lower envelope (excluding positive excursions). When Optim
# is loaded, OptimExt provides `optim_background_trend`, a full smooth zero-peak fit, used instead.
function _background_trend(log_f, log_s, components, transition_width)
    ext = Base.get_extension(@__MODULE__, :OptimExt)
    if ext !== nothing
        return Base.invokelatest(ext.optim_background_trend, log_f, log_s, components, transition_width)
    end
    breaks = [c.log_f_stop for c in components][1:(end - 1)]
    basis = hcat(ones(length(log_f)), log_f, (max.(log_f .- b, 0.0) for b in breaks)...)
    coef = basis \ log_s
    for _ in 1:3   # iteratively refit through the lower envelope so peaks don't lift the trend
        resid = log_s .- basis * coef
        mask = _lower_envelope_mask(resid)
        count(mask) < size(basis, 2) && break   # need enough points to determine the fit
        coef = basis[mask, :] \ log_s[mask]
    end
    return basis * coef
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
        peak_width_limits = (minimum(diff(log_f)), (maximum(log_f) - minimum(log_f)) / 2),
        minprom = nothing,
        kwargs...
    )
    log_A = first(log_s) # Estimate of amplitude

    β = last([ones(length(log_f)) log_f] \ log_s) # Simple linear regression. Start by guessing all components have the same exponent, and evenly distribute the breaks
    log_f_stop = range(extrema(log_f)..., length = components + 1)[2:end]
    transition_width = (maximum(log_f) - minimum(log_f)) / (20)

    components = map(1:components) do i
        ComponentArray(; log_f_stop = log_f_stop[i], β = β)
    end

    # * Find peaks on the residual after removing the zero-peak background. `_background_trend`
    #   returns a robust least-squares piecewise-linear power law by default, or a full Optim
    #   zero-peak fit when Optim is loaded. Peaks that ride a steep, curved background — and so
    #   are not local maxima of the raw spectrum — stand out cleanly on the residual.
    trend = _background_trend(log_f, log_s, components, transition_width)
    residual = ToolsArray(log_s .- trend, Log10𝑓(log_f))

    # For a fixed count, take the strongest peaks with no threshold (the count is the
    # selection). For `:auto`, gate on a data-driven prominence threshold.
    auto = !(peaks isa Integer)
    auto && isnothing(minprom) && (minprom = peak_threshold * _rstd(parent(residual)))
    _, proms, bounds = findpeaks(residual, w; minprom, kwargs...)

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

    # Initialise peak amplitudes conservatively, anchored to the spectral floor rather than the
    # local background: a prominence-scaled bump above `minimum(log_s)`. On a steep background
    # the local-background lift would make a peak start comparable to the (huge) background and
    # hijack the joint refine — capturing slope curvature and biasing β; a small start lets the
    # background fit first and the refine then grows genuine peaks into the residual.
    # NB: use fresh names inside this closure — assigning `log_A`/`log_σ` here would leak to the
    # enclosing `log_A = first(log_s)` (a `do` block shares enclosing locals).
    floorlevel = minimum(log_s)
    peakparams = map(proms, bounds) do prom, bound
        pf = mean(bound)
        pσ = (maximum(bound) - minimum(bound)) / 2
        pσ ≤ 0 && (pσ = transition_width) # A guess
        pA = floorlevel + log10(max(expm1(prom * log(10)), eps()))
        return ComponentArray(; log_f = pf, log_σ = pσ, log_A = pA)
    end

    return ComponentArray(; log_A, peaks = _emptylike(peakparams), components, transition_width)
end
