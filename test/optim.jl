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
    params = ComponentArray(; components, peaks, transition_width = 0.15, log_A = 1.0)
end

begin # * Generate test data
    log_f = range(0, 3, length = 500)
    f = map(exp10, log_f)
    s = mapple(f, params)
    log_s = log10.(s) .+ 0.1 * randn(length(s))
    s = exp10.(log_s)
    lines(f, s; axis = (; xscale = log10, yscale = log10))
end

begin # * Try fit
    log_s = log10.(s)
    fitted_params = fit_mapple(log_f, log_s; components = 2, peaks = 0, w = 50)
    refined_params = fit_mapple(log_f, log_s, fitted_params;
                                autodiff = Optim.ADTypes.AutoForwardDiff())

    a = @timed fit_mapple(log_f, log_s, fitted_params;
                          autodiff = Optim.ADTypes.AutoFiniteDiff())
    b = @timed fit_mapple(log_f, log_s, fitted_params;
                          autodiff = Optim.ADTypes.AutoForwardDiff()) # So much faster
    @test b.time < a.time / 2
    @test b.bytes < a.bytes

    fig = Figure()
    ax = Axis(fig[1, 1];
              xlabel = "Frequency (Hz)",
              ylabel = "Power",
              title = "Spectrum Fitting with mapple Model",
              xscale = log10, yscale = log10)
    lines!(ax, f, s; label = "Data", color = :black, linewidth = 2, alpha = 0.4)
    fitted_s = mapple(f, fitted_params)
    lines!(ax, f, fitted_s; label = "Fit", color = :crimson, linewidth = 2)
    refined_s = mapple(f, refined_params)
    lines!(ax, f, refined_s;
           label = "Refined Fit", color = :black, linewidth = 2)
    axislegend(ax, position = :lb)
    display(fig)

    @test cor(refined_s, s) > 0.9
end

begin # * Test parameter agreement
    @test params.log_A≈refined_params.log_A atol=0.05
    # @test params.transition_width≈refined_params.transition_width atol=0.15
    @test params.components.β≈refined_params.components.β atol=0.1
    @test sort(params.components.log_f_stop)[1]≈sort(refined_params.components.log_f_stop)[1] atol=0.05
    # @test sort(params.peaks.log_f)≈sort(refined_params.peaks.log_f) atol=0.05
end

begin # * Test with overlapping peaks
    components = ComponentArray([
                                    ComponentArray(; β = -4.0, log_f_stop = 7.0)])
    peaks = ComponentArray([
                               ComponentArray(; log_f = 1.3, log_σ = 0.1, log_A = 2.0),
                               ComponentArray(; log_f = 1.2, log_σ = 1.0, log_A = 0.05)])
    params = ComponentArray(; components, peaks, transition_width = 0.1, log_A = 1.0)
    log_f = range(0, 3, length = 500)
    f = map(exp10, log_f)
    s = mapple(f, params)
    log_s = log10.(s) .+ 0.1 * randn(length(s))
    s = exp10.(log_s)
    fitted_params = fit_mapple(log_f, log_s; components = 1, peaks = 2, w = 50)
    refined_params = fit_mapple(log_f, log_s, fitted_params;
                                autodiff = Optim.ADTypes.AutoForwardDiff())

    fig = Figure()
    ax = Axis(fig[1, 1];
              xlabel = "Frequency (Hz)",
              ylabel = "Power",
              title = "Spectrum Fitting with mapple Model",
              xscale = log10, yscale = log10)
    lines!(ax, f, s; label = "Data", color = :black, linewidth = 2, alpha = 0.4)
    fitted_s = mapple(f, fitted_params)
    lines!(ax, f, fitted_s; label = "Fit", color = :crimson, linewidth = 2)
    refined_s = mapple(f, refined_params)
    lines!(ax, f, refined_s;
           label = "Refined Fit", color = :black, linewidth = 2)
    axislegend(ax, position = :lb)
    display(fig)

    begin # * Compare individual params
        # Peak frequencies
        @test sort(params.peaks.log_f)≈sort(refined_params.peaks.log_f) atol=0.05

        # * Peak amplitudes
        @test sort(params.peaks.log_A)≈sort(refined_params.peaks.log_A) atol=0.05

        # * Peak widths
        @test sort(params.peaks.log_σ)≈sort(refined_params.peaks.log_σ) atol=0.05

        # * Component slopes
        @test params.components.β≈refined_params.components.β atol=0.1
    end
end
