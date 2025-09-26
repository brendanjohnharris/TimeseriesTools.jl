using TimeseriesTools
using Optim
using Statistics
using Test
using Statistics
using CairoMakie
import CairoMakie: Axis
using ComponentArrays
using ForwardDiff

begin # * Generate test data
    components = ComponentArray([
                                    ComponentArray(; β = 2.0, log_f_stop = 1.5),
                                    ComponentArray(; β = 4.0, log_f_stop = 5.0)])
    peaks = ComponentArray([
                               ComponentArray(; log_f = 0.7, log_σ = 0.1, log_A = 1.5),
                               ComponentArray(; log_f = 2.5, log_σ = 0.1, log_A = 0.5)])
    params = ComponentArray(; components, peaks, transition_width = 0.2, log_A = 1.0)
end

begin # * Generate test data
    log_f = range(0, 3, length = 500)
    f = map(exp10, log_f)
    s = apple(f, params)
    log_s = log10.(s) .+ 0.1 * randn(length(s))
    s = exp10.(log_s)
    lines(f, s; axis = (; xscale = log10, yscale = log10))
end

begin # * Try fit
    log_s = log10.(s)
    fitted_params = fit_apple(log_f, log_s; components = 2, peaks = 2, w = 50)
    refined_params = fit_apple(log_f, log_s, fitted_params; autodiff = :forward)

    a = @timed fit_apple(log_f, log_s, fitted_params; autodiff = :finite)
    b = @timed fit_apple(log_f, log_s, fitted_params; autodiff = :forward) # So much faster
    @test b.time < a.time / 3
    @test b.bytes < a.bytes

    fig = Figure()
    ax = Axis(fig[1, 1];
              xlabel = "Frequency (Hz)",
              ylabel = "Power",
              title = "Spectrum Fitting with apple Model",
              xscale = log10, yscale = log10)
    lines!(ax, f, s; label = "Data", color = :black, linewidth = 2, alpha = 0.4)
    fitted_s = apple(f, fitted_params)
    lines!(ax, f, fitted_s; label = "Fit", color = :crimson, linewidth = 2)
    refined_s = apple(f, refined_params)
    lines!(ax, f, refined_s;
           label = "Refined Fit", color = :black, linewidth = 2)
    axislegend(ax, position = :rt)
    display(fig)

    @test cor(refined_s, s) > 0.9
end

begin # * Test parameter agreement
    @test params.log_A≈refined_params.log_A atol=0.05
    @test params.transition_width≈refined_params.transition_width atol=0.05
    @test params.components.β≈refined_params.components.β atol=0.05
    @test sort(params.components.log_f_stop)[1]≈sort(refined_params.components.log_f_stop)[1] atol=0.05
    @test sort(params.peaks.log_f)≈sort(refined_params.peaks.log_f) atol=0.05
end
