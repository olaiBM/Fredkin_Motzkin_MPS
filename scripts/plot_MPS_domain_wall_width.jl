using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2
using Plots
using LaTeXStrings
using Printf

pgfplotsx()
include(srcdir("helpers.jl"))

default(
   fontfamily = "Computer Modern",
   legendfontsize = 17,
   guidefontsize = 22,
   tickfontsize = 20,
   lw = 1.8,
   markersize = 5,
   size = (600, 400),
   grid = true
)


N = 2000
r_vals = collect(1:2000)
qg1_vals = [1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.0035, 1.004, 1.0045, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01]

qg1_data_dict_F = load_q_data(N, qg1_vals, "Fredkin")
qg1_curve_fit_F = fermi_dirac_fit(r_vals, qg1_data_dict_F, "Fredkin")
qg1_data_dict_M = load_q_data(N, qg1_vals, "Motzkin")
qg1_curve_fit_M = fermi_dirac_fit(r_vals, qg1_data_dict_M, "Motzkin")


### Fit 

qg1_xi_pairs_F = collect(pairs(qg1_curve_fit_F))
slope_F, intercept_F = fit_line(log.(log.(first.(qg1_xi_pairs_F))), log.(last.(qg1_xi_pairs_F)))
qg1_xi_pairs_M = collect(pairs(qg1_curve_fit_M))
slope_M, intercept_M = fit_line(log.(log.(first.(qg1_xi_pairs_M))), log.(last.(qg1_xi_pairs_M)))

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
F_fit_label = latexstring(L"y\propto\;", slope_F_string, L"\log(-\tau) ")
F_fit_vals = slope_F .* log.(log.(qg1_vals)) .+ intercept_F
M_fit_label = latexstring(L"y\propto\;", slope_M_string, L"\log(-\tau) ")
M_fit_vals = slope_M .* log.(log.(qg1_vals)) .+ intercept_M



# Plot Fredkin
scatter(log.(log.(first.(qg1_xi_pairs_F))), log.(last.(qg1_xi_pairs_F)),
    xlabel = L"\log(-\tau)",
    ylabel = L"\log(w(\tau))",
    label = L"\log(w(\tau))",
    marker = :circle,
    color = :black,
    legend = :topright,
    lw = 2)

plot!(log.(log.(qg1_vals)), F_fit_vals,
    label = F_fit_label,
    linestyle = :dash,
    color = :gray,
    lw = 2)

savefig(plotsdir("Correlation_length_MPS", savename("Fredkin_domain_wall_width_MPS", (; N = N), "pdf")))

# Plot Motzkin
scatter(log.(log.(first.(qg1_xi_pairs_M))), log.(last.(qg1_xi_pairs_M)),
    xlabel = L"\log(-\tau)",
    ylabel = L"\log(w(\tau))",
    label = L"\log(w(\tau))",
    marker = :circle,
    color = :black,
    legend = :topright,
    lw = 2)

plot!(log.(log.(qg1_vals)), M_fit_vals,
    label = M_fit_label,
    linestyle = :dash,
    color = :gray,
    lw = 2)

savefig(plotsdir("Correlation_length_MPS", savename("Motzkin_domain_wall_width_MPS", (; N = N), "pdf")))