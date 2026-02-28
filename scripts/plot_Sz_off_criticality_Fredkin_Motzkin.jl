using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2
using Plots
using LaTeXStrings

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
r_vals = collect(1:N)
ql1_vals = [0.998, 0.99, 0.98]
qg1_vals = [1.002, 1.01, 1.02]

ql1_linestyle_dict = Dict(ql1_vals[1] => :dash, ql1_vals[2] => :dot, ql1_vals[3] => :solid)
qg1_linestyle_dict = Dict(qg1_vals[1] => :dash, qg1_vals[2] => :dot, qg1_vals[3] => :solid)


# Contains start and stop indices for N = 2000 Frekdin and Motzkin, used to determine the window of values used for the line fitting, determined by plotting and testing
ql1_start_stop_index_F_N_2000 = Dict(ql1_vals[1] => [6, 400], ql1_vals[2] => [6, 125], ql1_vals[3] => [6, 75]) 
ql1_start_stop_index_M_N_2000 = Dict(ql1_vals[1] => [6, 300], ql1_vals[2] => [6, 110], ql1_vals[3] => [6, 50]) 



### Loading, filtering and fitting for q<1 data
ql1_filtered_data_dict_F = load_and_filter_ql1_data(N, ql1_vals, "Fredkin") ### We are doing a semilog plot, so negative values cannot be included
ql1_slope_and_intercept_dict_F = return_slope_and_intercept_ql1(ql1_filtered_data_dict_F, ql1_start_stop_index_F_N_2000)
ql1_filtered_data_dict_M = load_and_filter_ql1_data(N, ql1_vals, "Motzkin") ### We are doing a semilog plot, so negative values cannot be included
ql1_slope_and_intercept_dict_M = return_slope_and_intercept_ql1(ql1_filtered_data_dict_M, ql1_start_stop_index_M_N_2000)



### loading and fitting for q>1 data


qg1_data_dict_F = load_q_data(N, qg1_vals, "Fredkin")
qg1_curve_fit_F = fermi_dirac_fit(r_vals, qg1_data_dict_F, "Fredkin")
qg1_data_dict_M = load_q_data(N, qg1_vals, "Motzkin")
qg1_curve_fit_M = fermi_dirac_fit(r_vals, qg1_data_dict_M, "Motzkin")







#### PLOTTING Fredkin q<1 ####
p = plot( xlim = (0, 500),ylim = (-20, 0),  xlabel = L"r", ylabel = L"\log\langle S^{z}_{r} \rangle",legend = :topright)
for q in ql1_vals
    slope, intercept= ql1_slope_and_intercept_dict_F[q]
    linear_fit_vals = [linear_fit(r, slope, intercept) for r in r_vals]

    plot!(r_vals, linear_fit_vals,
        label = false,
        color = :gray,
        linestyle = :solid,
        lw = 1)
    
        fit_label_F = latexstring(L"q = ", q, L"\;\mathrm{F}")
        if ql1_linestyle_dict[q] == :dash
            plot!(ql1_filtered_data_dict_F[q][1], log.(ql1_filtered_data_dict_F[q][2]), label = fit_label_F, color = :red, linestyle = ql1_linestyle_dict[q], dash_pattern= "on 2mm off 2mm")
        else
            plot!(ql1_filtered_data_dict_F[q][1], log.(ql1_filtered_data_dict_F[q][2]), label = fit_label_F, color = :red, linestyle = ql1_linestyle_dict[q])
        end
    
    
end

#### PLOTTING Motzkin q<1 ####
for q in ql1_vals
    slope, intercept= ql1_slope_and_intercept_dict_M[q]
    linear_fit_vals = [linear_fit(r, slope, intercept) for r in r_vals]
    
    plot!(r_vals, linear_fit_vals,
        label = false,
        color = :gray,
        linestyle = :solid,
        lw = 1)
    
        fit_label_M = latexstring(L"q = ", q, L"\;\mathrm{M}")
        if ql1_linestyle_dict[q] == :dash
            plot!(ql1_filtered_data_dict_M[q][1], log.(ql1_filtered_data_dict_M[q][2]), label = fit_label_M, color = :blue, linestyle = ql1_linestyle_dict[q], dash_pattern= "on 2mm off 2mm")
        else
            plot!(ql1_filtered_data_dict_M[q][1], log.(ql1_filtered_data_dict_M[q][2]), label = fit_label_M, color = :blue, linestyle = ql1_linestyle_dict[q])
        end
    
    
end


savefig(plotsdir("Correlation_length_MPS", savename("Sz_ql1", (; N = N), "pdf")))



#### PLOTTING Fredkin q>1 ####
p = plot(ylim = (-1.1, 1.1),  xlabel = L"r", ylabel = L"\langle S^{z}_{r} \rangle",legend = :topright)
for q in qg1_vals
    domain_wall_width = qg1_curve_fit_F[q]
    fermi_dirac_fits = fermi_dirac_model_F(r_vals, [domain_wall_width])
    plot!(r_vals, fermi_dirac_fits,
        label = false,
        color = :gray,
        linestyle = :solid,
        lw = 1)
    
        fit_label_F = latexstring(L"q = ", q, L"\;\mathrm{F}")
        if qg1_linestyle_dict[q] == :dash
            plot!(r_vals, qg1_data_dict_F[q], label = fit_label_F, color = :red, linestyle = qg1_linestyle_dict[q], dash_pattern= "on 2mm off 2mm")
        else
            plot!(r_vals, qg1_data_dict_F[q], label = fit_label_F, color = :red, linestyle = qg1_linestyle_dict[q])
        end
    
    
end


#### PLOTTING Motzkin q>1 ####
for q in qg1_vals
    domain_wall_width = qg1_curve_fit_M[q]
    fermi_dirac_fits = fermi_dirac_model_M(r_vals, [domain_wall_width])
    plot!(r_vals, fermi_dirac_fits,
        label = false,
        color = :gray,
        linestyle = :solid,
        lw = 1)
    
        fit_label_M = latexstring(L"q = ", q, L"\;\mathrm{M}")
        if qg1_linestyle_dict[q] == :dash
            plot!(r_vals, qg1_data_dict_M[q], label = fit_label_M, color = :blue, linestyle = qg1_linestyle_dict[q], dash_pattern= "on 2mm off 2mm")
        else
            plot!(r_vals, qg1_data_dict_M[q], label = fit_label_M, color = :blue, linestyle = qg1_linestyle_dict[q])
        end
    
    
end
savefig(plotsdir("Correlation_length_MPS", savename("Sz_qg1", (; N = N), "pdf")))