using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2

include(srcdir("transfer_matrix.jl"))





L = 4000
q_vals_qg1 = [1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.0035, 1.004, 1.0045, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01]
q_vals_ql1 = [0.9995, 0.999,0.9985, 0.998, 0.9975, 0.997, 0.9965, 0.996, 0.9955, 0.995, 0.994, 0.993, 0.992, 0.991, 0.99]


### Arrays for correlation lengths xi
xi_vals_qg1_F = Vector{Float64}(undef, length(q_vals_qg1))
xi_vals_qg1_M = Vector{Float64}(undef, length(q_vals_qg1))
xi_vals_ql1_F = Vector{Float64}(undef, length(q_vals_qg1))
xi_vals_ql1_M = Vector{Float64}(undef, length(q_vals_qg1))



#Compute correlation lengths
@time begin
    for i in eachindex(q_vals_qg1)
        Fredkin_T = tridiag_transfer_matrix_F(L, q_vals_qg1[i])
        Motzkin_T = tridiag_transfer_matrix_M(L, q_vals_qg1[i])
        xi_vals_qg1_F[i]= correlation_length(Fredkin_T)
        xi_vals_qg1_M[i]= correlation_length(Motzkin_T)
       
    end
end


@time begin
    for i in eachindex(q_vals_ql1)
        Fredkin_T = tridiag_transfer_matrix_F(L, q_vals_ql1[i])
        Motzkin_T = tridiag_transfer_matrix_M(L, q_vals_ql1[i])
        xi_vals_ql1_F[i]= correlation_length(Fredkin_T)
        xi_vals_ql1_M[i]= correlation_length(Motzkin_T)
       
    end
end


@save datadir("sims", savename((; L=L), "q>1_Fredkin_correlation_lengths.jld2")) xi_vals_qg1_F
@save datadir("sims", savename((; L=L), "q<1_Fredkin_correlation_lengths.jld2")) xi_vals_ql1_F
@save datadir("sims", savename((; L=L), "q>1_Motzkin_correlation_lengths.jld2")) xi_vals_qg1_M
@save datadir("sims", savename((; L=L), "q<1_Motzkin_correlation_lengths.jld2")) xi_vals_ql1_M