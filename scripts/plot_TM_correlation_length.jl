using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2

using Plots
using LaTeXStrings

pgfplotsx()

default(
   fontfamily = "Computer Modern",
   legendfontsize = 17,
   guidefontsize = 22,
   tickfontsize = 20,
   lw = 1.8,
   markersize = 6.5,
   size = (600, 400),
   grid = true
)


L = 4000
q_vals_qg1 = [1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.0035, 1.004, 1.0045, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01]
q_vals_ql1 = [0.9995, 0.999,0.9985, 0.998, 0.9975, 0.997, 0.9965, 0.996, 0.9955, 0.995, 0.994, 0.993, 0.992, 0.991, 0.99]


@load datadir("sims", savename((; L=L), "qg1_Fredkin_correlation_lengths.jld2")) xi_vals_qg1_F
@load datadir("sims", savename((; L=L), "ql1_Fredkin_correlation_lengths.jld2")) xi_vals_ql1_F
@load datadir("sims", savename((; L=L), "qg1_Motzkin_correlation_lengths.jld2")) xi_vals_qg1_M
@load datadir("sims", savename((; L=L), "ql1_Motzkin_correlation_lengths.jld2")) xi_vals_ql1_M


### Prepare for plotting 
x_qg1 = log.(log.(q_vals_qg1))
x_ql1 = log.(abs.(log.(q_vals_ql1)))

y_F_qg1 = log.(xi_vals_qg1_F)
y_F_ql1 = log.(xi_vals_ql1_F)

y_M_qg1 = log.(xi_vals_qg1_M)
y_M_ql1 = log.(xi_vals_ql1_M)

# Design matrix: column of ones (intercept) and x
X_qg1 = [ones(length(x_qg1)) x_qg1]
X_ql1 = [ones(length(x_ql1)) x_ql1]

# Solve least squares: β = (X'X)⁻¹ X'y
β_F_qg1 = X_qg1 \ y_F_qg1
β_F_ql1 = X_ql1 \ y_F_ql1

β_M_qg1 = X_qg1 \ y_M_qg1
β_M_ql1 = X_ql1 \ y_M_ql1


intercept_F_qg1, slope_F_qg1 = β_F_qg1
intercept_F_ql1, slope_F_ql1 = β_F_ql1

intercept_M_qg1, slope_M_qg1 = β_M_qg1
intercept_M_ql1, slope_M_ql1 = β_M_ql1


println("Fit: y_F_qg1 ≈ $slope_F_qg1 * x_F_qg1 + $intercept_F_qg1")
println("Fit: y_F_ql1 ≈ $slope_F_ql1 * x_F_ql1 + $intercept_F_ql1")
println("Fit: y_M_qg1 ≈ $slope_M_qg1 * x_M_qg1 + $intercept_M_qg1")
println("Fit: y_M_ql1 ≈ $slope_M_ql1 * x_M_ql1 + $intercept_M_ql1")

# Define fit line values
x_fit_qg1 = range(minimum(x_qg1), maximum(x_qg1), length=200)
x_fit_ql1 = range(minimum(x_ql1), maximum(x_ql1), length=200)

F_y_fit_qg1 = intercept_F_qg1 .+ slope_F_qg1 .* x_fit_qg1
F_y_fit_ql1 = intercept_F_ql1 .+ slope_F_ql1 .* x_fit_ql1

M_y_fit_qg1 = intercept_M_qg1 .+ slope_M_qg1 .* x_fit_qg1
M_y_fit_ql1 = intercept_M_ql1 .+ slope_M_ql1 .* x_fit_ql1


################################### Fredkin plot ###################################
plot(
    x_qg1, y_F_qg1,
    seriestype = :scatter, label = L"q>1", 
    xlabel = L"\log(|\tau|)", 
    ylabel = L"\log(\xi_{\mathrm{TM}}(\tau))", 
    color = :red,
    legend = :bottomleft, markershape = :circle
)


# Plot together
plot!(x_ql1,y_F_ql1,
    seriestype = :scatter, label = L"q<1",
    color = :blue, markershape = :utriangle
)


plot!(x_fit_qg1, F_y_fit_qg1, label = L"\mathrm{Fit},~\nu = %$(abs(round(slope_F_qg1, digits=2)))", linestyle = :dot, color = :gray, lw = 2)
savefig(plotsdir("Correlation_length_TM", savename((; L=L), "_correlation_length_TM_Fredkin.pdf")))

################################### Motzkin plot ###################################


plot(
    x_qg1, y_M_qg1,
    seriestype = :scatter, label = L"q>1", 
    xlabel = L"\log_{10}(|\tau|)", 
    ylabel = L"\log_{10}(\xi_{\mathrm{TM}}(\tau))", 
    color = :red,
    legend = :bottomleft, markershape = :circle
)


# Plot together
plot!(x_ql1,y_M_ql1,
    seriestype = :scatter, label = L"q<1",
    color = :blue, markershape = :utriangle
)


plot!(x_fit_qg1, M_y_fit_qg1, label = L"\mathrm{Fit},~\nu = %$(round(slope_M_qg1, digits=2))", linestyle = :dot, color = :gray, lw = 2)
savefig(plotsdir("Correlation_length_TM", savename((; L=L), "_correlation_length_TM_Motzkin.pdf")))

#########################################################################################################