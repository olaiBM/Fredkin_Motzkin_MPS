using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"

using Plots
using LaTeXStrings
using JLD2

pgfplotsx()

default(
   fontfamily = "Computer Modern",
   legendfontsize = 17,
   guidefontsize = 22,
   tickfontsize = 20,
   lw = 1.8,
   markersize = 4.5,
   size = (600, 400),
   grid = true
)


# Parameters
L_start = 1000
L_final = 10000


Ls = L_start:1000:L_final
qg1s =[1.001, 1.1, 1.2,  1.5] #[1.001, 1.1,  1.5]
ql1s = [0.999, 0.9, 0.8, 0.5]#[0.999, 0.9,0.5]

params_qg1 = (
    regime = "qg1",
    L_start = L_start,
    L_final = L_final,
)

params_ql1 = (
    regime = "ql1",
    L_start = L_start,
    L_final = L_final,
)

@load datadir("sims", savename("Fredkin_gap", params_qg1, "jld2")) Ls qg1s F_gap_dict_qg1
@load datadir("sims", savename("Motzkin_gap", params_qg1, "jld2")) Ls qg1s M_gap_dict_qg1
@load datadir("sims", savename("Fredkin_gap", params_ql1, "jld2")) Ls ql1s F_gap_dict_ql1
@load datadir("sims", savename("Motzkin_gap", params_ql1, "jld2")) Ls ql1s M_gap_dict_ql1

# For q > 1
markerdict_qg1 = Dict(
    1.001 => :dtriangle,
    1.1 => :utriangle,
    1.2   => :circle,
    1.5   => :square,
)

# For q < 1
markerdict_ql1 = Dict(
    0.999 => :dtriangle,
    0.9 => :utriangle,
    0.8   => :circle,
    0.5   => :square,
)


# end
# display(current())

# ===========================
# Left part: (a) Fredkin + Motzkin eigenvalue gaps
# ===========================
p1 = plot(layout=(1,2), ylim=(-0.2,1),
          xticks=([1e3, 5e3, 1e4],
                  [L"10^{3}", L"5\times10^{3}", L"10^{4}"])
)

xlabel!(p1[1], L"L")
ylabel!(p1[1], L"\Delta(q, L)")

xlabel!(p1[2], L"L")
yticks!(p1[2], yticks(p1[1])[1], ["", "", "", "", ""])
# --- Fredkin ---
for current_q in qg1s
    F_gap_q = F_gap_dict_qg1[current_q]
    plot!(p1[1], Ls, F_gap_q,
          color=:red, marker=markerdict_qg1[current_q], label="")
end
for current_q in ql1s
    F_gap_q = F_gap_dict_ql1[current_q]
    plot!(p1[1], Ls, F_gap_q,
          color=:blue, linestyle=:dot, marker=markerdict_ql1[current_q], label="")
end
title!(p1[1], "Fredkin chain")

# --- Motzkin ---
for current_q in qg1s
    M_gap_q = M_gap_dict_qg1[current_q]
    plot!(p1[2], Ls, M_gap_q,
          color=:red, marker=markerdict_qg1[current_q],
          label="q = $current_q")
end
for current_q in ql1s
    M_gap_q = M_gap_dict_ql1[current_q]
    plot!(p1[2], Ls, M_gap_q,
          color=:blue, linestyle=:dot, marker=markerdict_ql1[current_q],
          label="q = $current_q")
end
title!(p1[2], "Motzkin chain")

savefig(plotsdir("Correlation_length_TM", savename("gaps", (L_start = L_start, L_final = L_final), "pdf")))
