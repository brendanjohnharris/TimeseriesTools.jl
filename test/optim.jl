using TimeseriesTools
using Optim
using Statistics
using Test
using Statistics
using CairoMakie
import CairoMakie: Axis
using ComponentArrays
using ForwardDiff

begin # * Try construction
    log_f = 0.001:0.001:3
    peaks = [
        (; log_f = log10(20), σ̃ = 0.1, log_A = -0.5),
        (; log_f = log10(100), σ̃ = 0.05, log_A = -1.0),
        (; log_f = log10(500), σ̃ = 0.05, log_A = -2.0)
    ]
    p = (; log_b = 1, β = 1.5, log_k = 1, log_c = -3,
         peaks)
    log_s = @inferred oneoneff(log_f, p)

    logspectrum = ToolsArray(log_s, Log10𝑓(log_f)) .+ 0.05 .* randn(length(log_s))

    n_peaks = length(p.peaks)

    lines(logspectrum)
end

begin # * First, estimate initial parameters from data
    params = fit_oneoneff(logspectrum; n_peaks = 3, w = 100)

    function plot_fit(logspectrum, params)
        fig = Figure(size = (800, 600))
        ax = Axis(fig[1, 1],
                  xlabel = "log₁₀(frequency)",
                  ylabel = "log₁₀(power)",
                  title = "Spectrum Fitting with 1/f Model")

        # Plot original data
        lines!(ax, logspectrum,
               label = "Data", color = :black, linewidth = 2, alpha = 0.5)

        # Plot initial guess
        fitted_log_s = oneoneff(lookup(logspectrum, 1), params)
        lines!(ax, lookup(logspectrum, 1), fitted_log_s, color = :crimson, linewidth = 2)

        axislegend(ax, position = :rt)
        return fig
    end
    plot_fit(logspectrum, params) |> display
end

begin # * Use Optim to refine these initial guesses
    a = @timed fit_oneoneff(logspectrum, params; autodiff = :finite)
    b = @timed fit_oneoneff(logspectrum, params; autodiff = :forward) # So much faster
    @test b.time < a.time / 3
    @test b.bytes < a.bytes

    fitted_params = @inferred fit_oneoneff(logspectrum, params; autodiff = :forward)
    # Calculate fitted spectrum
    fitted_log_s = @inferred oneoneff(log_f, fitted_params)

    @test cor(fitted_log_s, logspectrum) > 0.99

    # Plot results
    plot_fit(logspectrum, fitted_params) |> display
end
