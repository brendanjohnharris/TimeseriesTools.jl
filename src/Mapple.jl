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
"""
struct MAPPLE <: StatsAPI.RegressionModel
    params::ComponentArray
end

function MAPPLE(; peaks, components, log_A)
    ComponentArray(; log_A, peaks, components) |> MAPPLE
end

function mapple_sort!(params)
    sort!(params.peaks; by = p -> p.log_f)
    sort!(params.components; by = c -> c.log_f_stop)
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
    if maximum(df) > minimum(df) * 10
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
    spectrum = ToolsArray(s, f)
    StatsAPI.fit(MAPPLE, spectrum; kwargs...)
end

StatsAPI.params(m::MAPPLE) = m.params
peaks(m::MAPPLE) = m.params.peaks
components(m::MAPPLE) = m.params.components

function mapple(f::AbstractVector, model::ComponentArray{T}) where {T}
    components = model[[:log_A, :components]]
    peaks = model[[:peaks]]
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, components, peaks)
    return s
end
function mapple(f, component_params::ComponentArray{T},
                peaks::ComponentArray{F}) where {T, F <: AbstractFloat}
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, component_params, peaks)
    return s
end
function mapple(f, component_params::ComponentArray{F},
                peaks::ComponentArray{T}) where {T, F <: AbstractFloat}
    ElType = promote_type(eltype(f), T)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, component_params, peaks)
    return s
end
function mapple(f, component_params::ComponentArray{F},
                peaks::ComponentArray{F}) where {F <: AbstractFloat}
    ElType = promote_type(eltype(f), F)
    s = similar(f, ElType)
    fill!(s, zero(ElType))
    mapple!(s, f, component_params, peaks)
    return s
end

function mapple!(s::AbstractVector{El}, f, component_params, peaks) where {El}
    fill!(s, zero(El))

    components = component_params.components
    peaks = peaks.peaks

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
    log_f = log10.(f)
    log_f_stops = components.log_f_stop
    @assert issorted(log_f)

    @inbounds @fastmath for i in eachindex(f, s)
        idx = searchsortedfirst(log_f_stops, log_f[i])
        idx = idx > n_components ? n_components : idx
        s[i] = component_amplitudes[idx] * f[i]^components.β[idx]
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
end

function fit_mapple(log_f, log_s;
                    w = max(1, length(log_f) ÷ 100),
                    peaks,
                    components,
                    minprom = (maximum(log_s) - minimum(log_s)) / 50,
                    kwargs...)
    logspectrum = ToolsArray(log_s, Log10𝑓(log_f))
    log_A = first(log_s) # Estimate of amplitude

    β = last([ones(length(log_f)) log_f] \ log_s) # Simple linear regression. Start by guessing all components have the same exponent, and evenly distribute the breaks
    log_f_stop = range(extrema(log_f)..., length = components + 1)[2:end]
    # transition_width = (maximum(log_f) - minimum(log_f)) / (20)
    # transition_width = max(transition_width, minimum(diff(log_f)))

    components = map(1:components) do i
        ComponentArray(; log_f_stop = log_f_stop[i], β = β)
    end

    # * Find peaks by looking for local maxima
    _, proms, bounds = findpeaks(logspectrum, w; minprom, kwargs...)

    if !isnothing(peaks) && !isempty(proms)
        if peaks > length(proms)
            proms = vcat(proms, [mean(log_s) for _ in 1:(peaks - length(proms))])
            bounds = vcat(bounds,
                          [deepcopy(first(bounds))
                           for _ in 1:(peaks - length(bounds))])
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
            log_σ = (maximum(log_f) - minimum(log_f)) / (20) # A guess
        end
        return ComponentArray(; log_f, log_σ, log_A)
    end

    return ComponentArray(; log_A, peaks, components)
end
