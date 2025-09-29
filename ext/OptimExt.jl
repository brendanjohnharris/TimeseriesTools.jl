module OptimExt
using Optim
using DimensionalData
using StatsAPI
using StatsBase
import TimeseriesTools: mapple, fit_mapple, MAPPLE, UnivariateSpectrum, Log10𝑓,
                        frequency_check, mapple_sort

function mapple_bounds(log_f, log_s, initial_params)
    lower = deepcopy(initial_params)
    upper = deepcopy(initial_params)

    lower.log_A = -Inf
    upper.log_A = 100 * (maximum(log_s) - minimum(log_s))

    lower.transition_width = minimum(diff(log_f)) / 4
    upper.transition_width = (maximum(log_f) - minimum(log_f))

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
function mapple_loss(params; f, log_s)
    components = params.components
    peaks = params.peaks

    # * Sort components and peaks by frequency
    cidxs = sortperm(components, by = c -> c.log_f_stop)
    pidxs = sortperm(peaks, by = p -> p.log_f)

    # Main prediction error
    pred = mapple(f, params)
    pred_log = map(log10, pred)
    loss = sum((log_s .- pred_log) .^ 2)

    # # Penalty for f_stops being too close together
    # separation_penalty = zero(eltype(params))
    # for i in eachindex(cidxs)[2:end]
    #     gap = components[cidxs[i]].log_f_stop - components[cidxs[i - 1]].log_f_stop
    #     if gap < min_log_f_separation
    #         # Quadratic penalty that increases as gap decreases
    #         separation_penalty += (min_log_f_separation - gap)^2
    #     end
    # end

    # # * Penalty for overlapping Gaussian peaks
    # overlap_penalty = zero(eltype(params))
    # n_peaks = length(pidxs)

    # # Pre-compute all peak parameters once
    # f_peaks = [exp10(params.peaks[pidxs[i]].log_f) for i in 1:n_peaks]
    # σ_peaks = [f_peaks[i] * tanh(params.peaks[pidxs[i]].log_σ) for i in 1:n_peaks]
    # A_peaks = [exp10(params.peaks[pidxs[i]].log_A) for i in 1:n_peaks]

    # # Pre-compute norms (sqrt(2π) cancels in normalized overlap)
    # norms = A_peaks .* σ_peaks

    # @inbounds @fastmath for i in 2:n_peaks
    #     # Use pre-computed values
    #     Δf = f_peaks[i] - f_peaks[i - 1]
    #     σ_i, σ_j = σ_peaks[i - 1], σ_peaks[i]

    #     # Early exit if peaks are far apart (>5σ separation = negligible overlap)
    #     if abs(Δf) > 5 * (σ_i + σ_j)
    #         continue
    #     end

    #     σ_combined_sq = σ_i^2 + σ_j^2
    #     σ_combined = sqrt(σ_combined_sq)

    #     # Simplified formula (constants cancel)
    #     overlap_integral = (norms[i - 1] * norms[i] / σ_combined) *
    #                        exp(-Δf^2 / (2 * σ_combined_sq))

    #     normalized_overlap = 2 * overlap_integral / (norms[i - 1] + norms[i])

    #     # Apply penalty
    #     if normalized_overlap > overlap_threshold
    #         overlap_penalty += (normalized_overlap - overlap_threshold)^2
    #     end
    # end

    # amplitude_penalty = zero(eltype(params))
    # for peak in peaks
    #     if peak.log_A < min_peak_log_A
    #         # Quadratic penalty that increases as amplitude decreases
    #         amplitude_penalty += (min_peak_log_A - peak.log_A)^2
    #     end
    # end

    # # * Penalty for similar slopes (β values)
    # slope_penalty = zero(eltype(params))
    # if length(cidxs) > 1
    #     βs = [params.components.β[i] for i in cidxs]
    #     β_scale = mean(abs.(βs))

    #     # Only check adjacent pairs
    #     for i in 2:length(cidxs)
    #         β_diff = abs(βs[i] - βs[i - 1])
    #         relative_diff = β_diff / (β_scale + 1e-6)
    #         if relative_diff < min_β_separation
    #             slope_penalty += (min_β_separation - relative_diff)^2
    #         end
    #     end
    # end

    return loss
end

mapple_loss(; kwargs...) = params -> mapple_loss(params; kwargs...)

function fit_mapple(log_f, log_s, initial_params;
                    algorithm = LBFGS(), autodiff = :forward, kwargs...) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)

    lower, upper = mapple_bounds(log_f, log_s, initial_params)

    min_peak_log_A = log10((exp10(maximum(log_s)) - exp10(minimum(log_s))) / 50)

    objective = mapple_loss(; f, log_s)

    result = optimize(objective, lower, upper, initial_params, Fminbox(algorithm),
                      Optim.Options(; kwargs...); autodiff)
    return Optim.minimizer(result)
end

"""
    fit!(m::MAPPLE, logspectrum; kwargs...)
Refine the parameters of a MAPPLE model `m` to fit the provided `spectrum`.
"""
function StatsAPI.fit!(m::MAPPLE, spectrum::AbstractDimVector{T, D};
                       kwargs...) where {T, d, D <: Tuple{<:d}}
    log_f = map(log10, lookup(spectrum, 1))
    log_s = map(log10, parent(spectrum))

    frequency_check(lookup(spectrum, 1), log_f)

    params = fit_mapple(log_f, log_s, m.params; kwargs...)
    m.params .= params
    sort!(m)
end
function StatsAPI.fit!(m::Type{MAPPLE}, f::AbstractVector, s::AbstractVector; kwargs...)
    spectrum = ToolsArray(s, f)
    StatsAPI.fit!(m, spectrum; kwargs...)
end

end
