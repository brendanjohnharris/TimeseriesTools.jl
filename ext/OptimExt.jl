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

function huber(x, y; δ)
    a = abs(x - y)
    if a > δ
        return δ * a - 0.5 * δ^2
    else
        return 0.5 * a^2
    end
end

function mapple_loss(params; f, log_s, min_log_f_separation, overlap_threshold, λ, δ)

    # * Sort components and peaks by frequency
    params = mapple_sort(params)
    components = params.components
    peaks = params.peaks

    # Main prediction error
    pred = mapple(f, params)
    pred_log = map(log10, pred)
    # mse_loss = sum((log_s .- pred_log) .^ 2)
    loss = sum(huber.(log_s, pred_log; δ))

    # Penalty for f_stops being too close together
    separation_penalty = zero(eltype(params))
    for i in eachindex(components)[2:end]
        gap = components[i].log_f_stop - components[i - 1].log_f_stop
        if gap < min_log_f_separation
            # Quadratic penalty that increases as gap decreases
            separation_penalty += (min_log_f_separation - gap)^2
        end
    end

    # Penalty for overlapping Gaussian peaks
    overlap_penalty = zero(eltype(params))
    n_peaks = length(peaks)
    for i in eachindex(peaks)[2:end]
        peak_i = params.peaks[i - 1]
        peak_j = params.peaks[i]

        # Calculate overlap between two Gaussians
        f_i = exp10(peak_i.log_f)
        f_j = exp10(peak_j.log_f)
        σ_i = f_i * tanh(peak_i.log_σ)
        σ_j = f_j * tanh(peak_j.log_σ)
        A_i = exp10(peak_i.log_A)
        A_j = exp10(peak_j.log_A)

        σ_combined = sqrt(σ_i^2 + σ_j^2)
        overlap_integral = (A_i * A_j * sqrt(2π) * σ_i * σ_j / σ_combined) *
                           exp(-(f_i - f_j)^2 / (2 * σ_combined^2))

        norm_i = A_i * σ_i * sqrt(2π)
        norm_j = A_j * σ_j * sqrt(2π)
        normalized_overlap = 2 * overlap_integral / (norm_i + norm_j)

        # Apply penalty if overlap exceeds threshold
        if normalized_overlap > overlap_threshold
            overlap_penalty += (normalized_overlap - overlap_threshold)^2
        end
    end

    # Combined loss with regularization terms
    total_loss = loss + λ * (separation_penalty + overlap_penalty)

    return total_loss
end

mapple_loss(; kwargs...) = params -> mapple_loss(params; kwargs...)

function fit_mapple(log_f, log_s, initial_params;
                    algorithm = LBFGS(), autodiff = :forward, kwargs...) # If you have ForwardDiff loaded, you can pass autodiff=:forward
    f = map(exp10, log_f)

    lower, upper = mapple_bounds(log_f, log_s, initial_params)

    objective = mapple_loss(; f, log_s,
                            min_log_f_separation = maximum(diff(log_f)),
                            overlap_threshold = 0.2,
                            λ = 10.0,
                            δ = 3 * median(abs.(diff(log_s))))

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
