using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2
using Printf
using Plots
using LaTeXStrings

pgfplotsx()
include(srcdir("helpers.jl"))

default(
   fontfamily = "Computer Modern",
   legendfontsize = 17,
   guidefontsize = 22,
   tickfontsize = 20,
   lw = 2,
   markersize = 5,
   size = (600, 400),
   grid = true
)


N =2000 # Must match the imported files
q_vals = [0.9995, 0.999,0.9985, 0.998, 0.9975, 0.997, 0.9965, 0.996, 0.9955, 0.995, 0.994, 0.993, 0.992, 0.991, 0.99]

#xis arrays contain correlation lengths
xis_F, slope_vals_F = compute_correlation_length_MPS(N, q_vals, "Fredkin")
xis_M, slope_vals_M = compute_correlation_length_MPS(N, q_vals, "Motzkin")

log_q_vals = log.(q_vals)
log_xis_F = log.(xis_F)
log_xis_M = log.(xis_M)

# Perform linear fit on the log-log data
log_xi_q_func_coeffs_F = fit_line(log.(-log_q_vals), log_xis_F)
log_xi_q_func_coeffs_M = fit_line(log.(-log_q_vals), log_xis_M)

slope_F, intercept_F = log_xi_q_func_coeffs_F
slope_M, intercept_M = log_xi_q_func_coeffs_M


# Round to 1 decimal
slope_F_string = @sprintf("%.2f", abs(slope_F))
intercept_F_string = @sprintf("%.1f", abs(intercept_F))
slope_F_sgn = slope_F < 0 ? L"-" : L"+"
slope_F_string = slope_F_sgn * slope_F_string
intercept_F_sgn = intercept_F < 0 ? L"-" : L"+"

slope_M_string = @sprintf("%.2f", abs(slope_M))
intercept_M_string = @sprintf("%.1f", abs(intercept_M))
slope_M_sgn = slope_M < 0 ? L"-" : L"+"
slope_M_string = slope_M_sgn * slope_M_string
intercept_M_sgn = intercept_M < 0 ? L"-" : L"+"


# Build LaTeX-style label: log10(ξ) = m * log10(1 - q) + b
fit_label_F = latexstring(L"y\propto\;", slope_F_string, L"\log(-\tau) ")
fit_label_M = latexstring(L"y\propto\;", slope_M_string, L"\log(\tau)")

# Fitted values
fit_vals_F = slope_F .* log.(-log_q_vals) .+ intercept_F
fit_vals_M = slope_M .* log.(-log_q_vals) .+ intercept_M

# Plot Fredkin 
scatter(log.(-log_q_vals), log_xis_F,
     xlabel = L"\log(\tau)",
     ylabel = L"\log(\xi(\tau))",
     label = L"\log(\xi(\tau))",
     marker = :circle,
     color = :black,
     legend = :topright)

plot!(log.(-log_q_vals), fit_vals_F,
      label = fit_label_F,
      linestyle = :dash,
      color = :gray,
      lw = 2)



savefig(plotsdir("Correlation_length_MPS", savename("Fredkin_correlation_length_MPS", (; N = N), "pdf")))

# Plot Motzkin
scatter(log.(-log_q_vals), log_xis_M,
     xlabel = L"\log(\tau)",
     ylabel = L"\log(\xi(\tau))",
     label = L"\log(\xi(\tau))",
     marker = :circle,
     color = :black,
     legend = :topright)

plot!(log.(-log_q_vals), fit_vals_M,
      label = fit_label_M,
      linestyle = :dash,
      color = :gray,
      lw = 2)

savefig(plotsdir("Correlation_length_MPS", savename("Motzkin_correlation_length_MPS", (; N = N), "pdf")))