using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"
using JLD2

include(srcdir("transfer_matrix.jl"))


# Parameters
L_start = 1000
L_final = 10000


Ls = L_start:1000:L_final
qg1s =[1.001, 1.1, 1.2,  1.5] #[1.001, 1.1,  1.5]
ql1s = [0.999, 0.9, 0.8, 0.5]#[0.999, 0.9,0.5]
params_qg1 = (
    regime = "q>1",
    L_start = L_start,
    L_final = L_final,
)

params_ql1 = (
    regime = "q<1",
    L_start = L_start,
    L_final = L_final,
)
# q > 1 gaps
@time begin
    F_gap_dict_qg1 = Dict{Float64, Vector{Float64}}()
    M_gap_dict_qg1 = Dict{Float64, Vector{Float64}}()
    for q in qg1s
        F_gap_dict_qg1[q] = [compute_gap(tridiag_transfer_matrix_F(L, q)) for L in Ls]
        M_gap_dict_qg1[q] = [compute_gap(tridiag_transfer_matrix_M(L, q)) for L in Ls]
    end
    @save datadir("sims",savename("Fredkin_gap", params_qg1, "jld2")) Ls qg1s F_gap_dict_qg1
    @save datadir("sims",savename("Motzkin_gap", params_qg1, "jld2")) Ls qg1s M_gap_dict_qg1
end

# q < 1 gaps
@time begin
    F_gap_dict_ql1 = Dict{Float64, Vector{Float64}}()
    M_gap_dict_ql1 = Dict{Float64, Vector{Float64}}()
    for q in ql1s
        F_gap_dict_ql1[q] = [compute_gap(tridiag_transfer_matrix_F(L, q)) for L in Ls]
        M_gap_dict_ql1[q] = [compute_gap(tridiag_transfer_matrix_M(L, q)) for L in Ls]
    end
    @save datadir("sims",savename("Fredkin_gap", params_ql1, "jld2")) Ls ql1s F_gap_dict_ql1
    @save datadir("sims",savename("Motzkin_gap", params_ql1, "jld2")) Ls ql1s M_gap_dict_ql1
end