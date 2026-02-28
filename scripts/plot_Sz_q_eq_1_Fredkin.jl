using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"

using Plots
using LaTeXStrings
using JLD2

pgfplotsx()

include(srcdir("helpers.jl"))


default(
    fontfamily = "Computer Modern",
    legendfontsize = 17,
    guidefontsize = 22,
    tickfontsize = 20,
    lw = 1.8,
    markersize = 5,
    size = (600, 250),
    grid = true
)


N = 1000
r_vals = collect(1:N)
q = 1.0
@load datadir("sims", savename((; N = N, q = q), "_Fredkin_Sz.jld2")) Sz_expect_data

Sz_expect_analytic = [expect_Sz_analytic_F(r_vals[i]) for i in 1:N]

####################### Plot of 50 first spins  ####################### 
plot(r_vals[1:50], Sz_expect_data[1:50],
     label = "MPS Fredkin",
     xlabel = L"r",
     ylabel = L"\langle S^{z}_{r} \rangle",
     color = :blue,
     linestyle = :solid,
     legend = :topright)

plot!(r_vals[1:50], Sz_expect_analytic[1:50],
     label = "TM Fredkin",
     color = :black,
     linestyle = :dot)


savefig(plotsdir("Sz_Fredkin", savename((; N=N, q=q), "_[1:50]_Fredkin_Sz.pdf")))



####################### Plot of first half  ####################### 
plot(r_vals[1:500], Sz_expect_data[1:500],
     label = "MPS Fredkin",
     xlabel = L"r",
     ylabel = L"\langle S^{z}_{r} \rangle",
     color = :blue,
     linestyle = :solid,
     legend = :topright)

plot!(r_vals[1:500], Sz_expect_analytic[1:500],
     label = "TM Fredkin",
     color = :black,
     linestyle = :dot)


savefig(plotsdir("Sz_Fredkin", savename((; N=N, q=q), "_[1:500]_Fredkin_Sz.pdf")))








