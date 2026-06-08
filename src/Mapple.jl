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
        println(io, round(comp.log_f_stop, digits = 4))
        printstyled(io, "    β:          ")
        println(io, round(comp.β, digits = 4))
    end

    # Peaks (Green)
    printstyled(io, "\nPeaks ($(length(m.params.peaks))):\n", color = :green, bold = true)
    for (i, peak) in enumerate(m.params.peaks)
        printstyled(io, "  Peak $i:\n")
        printstyled(io, "    log_f:  ")
        println(io, round(peak.log_f, digits = 4))
        printstyled(io, "    log_σ:  ")
        println(io, round(peak.log_σ, digits = 4))
        printstyled(io, "    log_A:  ")
        println(io, round(peak.log_A, digits = 4))
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

"""
    fit(::Type{MAPPLE}, spectrum::UnivariateSpectrum; kwargs...)

Roughly fit an MAPPLE model to linear frequencies and linear spectral density using
peak-finding and linear regression.
This should be done prior to using [`fit!`](@ref) to refine the parameters with Optim.jl
Please consider using 'logsample'd spectra for a fit that is less sensitive to
high-frequency noise
"""
function StatsAPI.fit(::Type{MAPPLE}, spectrum::AbstractDimVector; kwargs...)
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(log10, parent(spectrum))

    frequency_check(lookup(spectrum, 1), log_f)

    params = fit_mapple(log_f, log_s; kwargs...)
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

# Estimate the log-space background trend for peak detection. Default: a continuous
# piecewise-linear (hard-break) power law fitted by least squares over the component
# breakpoints — the actual background shape, so peaks on a steep, curved background stand out
# on the residual. When Optim is loaded, OptimExt provides `optim_background_trend`, a full
# smooth zero-peak fit, which is used instead.
function _background_trend(log_f, log_s, components, transition_width)
    ext = Base.get_extension(@__MODULE__, :OptimExt)
    if ext !== nothing
        return Base.invokelatest(ext.optim_background_trend, log_f, log_s, components, transition_width)
    end
    breaks = [c.log_f_stop for c in components][1:(end - 1)]
    basis = hcat(ones(length(log_f)), log_f, (max.(log_f .- b, 0.0) for b in breaks)...)
    return basis * (basis \ log_s)
end

function fit_mapple(
        log_f, log_s;
        w = max(1, length(log_f) ÷ 100),
        peaks,
        components,
        minprom = (maximum(log_s) - minimum(log_s)) / 50,
        kwargs...
    )
    logspectrum = ToolsArray(log_s, Log10𝑓(log_f))
    log_A = first(log_s) # Estimate of amplitude

    β = last([ones(length(log_f)) log_f] \ log_s) # Simple linear regression. Start by guessing all components have the same exponent, and evenly distribute the breaks
    log_f_stop = range(extrema(log_f)..., length = components + 1)[2:end]
    transition_width = (maximum(log_f) - minimum(log_f)) / (20)
    # transition_width = max(transition_width, minimum(diff(log_f)))

    components = map(1:components) do i
        ComponentArray(; log_f_stop = log_f_stop[i], β = β)
    end

    # * Find peaks on the residual after removing the zero-peak background. `_background_trend`
    #   returns a least-squares piecewise-linear power law by default, or a full Optim
    #   zero-peak fit when Optim is loaded. Peaks that ride a steep, curved background — and so
    #   are not local maxima of the raw spectrum — stand out cleanly on the residual.
    trend = _background_trend(log_f, log_s, components, transition_width)
    residual = ToolsArray(log_s .- trend, Log10𝑓(log_f))
    _, proms, bounds = findpeaks(residual, w; minprom, kwargs...)

    if !isnothing(peaks) && !isempty(proms)
        if peaks > length(proms)
            proms = vcat(proms, [mean(log_s) for _ in 1:(peaks - length(proms))])
            bounds = vcat(
                bounds,
                [
                    deepcopy(first(bounds))
                        for _ in 1:(peaks - length(bounds))
                ]
            )
        end

        ps = sortperm(proms; rev = true)[1:peaks]
        proms = proms[ps]
        bounds = bounds[ps]
    elseif peaks > 0
        @warn "Number of guessed peaks ($(length(proms))) does not match expected peaks ($peaks)"
        proms = zeros(peaks)
        df = maximum(log_f) - minimum(log_f)
        bounds = [(minimum(log_f) + df / 3) .. (maximum(log_f) - df / 3) for _ in 1:peaks]
    end

    peaks = map(proms, bounds) do prom, bound
        log_f = mean(bound)
        log_σ = (maximum(bound) - minimum(bound)) / 2
        s_f = logspectrum[Near(maximum(bound) + log_σ)]
        log_A = prom + s_f
        if log_σ ≤ 0
            log_σ = transition_width # A guess
        end
        return ComponentArray(; log_f, log_σ, log_A)
    end

    return ComponentArray(; log_A, peaks, components, transition_width)
end
