using DrWatson
@quickactivate "Fredkin_Motzkin_MPS"

using ITensors, ITensorMPS
using JLD2
include(srcdir("Fredkin_MPS.jl"))
include(srcdir("helpers.jl"))

### Variables
N = 200
sites_F = [Index(QN("Sz", 1) => 1, QN("Sz", -1) => 1; tags="S=1/2, Site,n=$i") for i in 1:N]


q_vals = [0.9995, 0.999,0.9985, 0.998, 0.9975, 0.997, 0.9965, 0.996, 0.9955, 0.995, 0.994, 0.993, 0.992, 0.991, 0.99, 
1.0, 1.0005, 1.001, 1.0015, 1.002, 1.0025, 1.003, 1.0035, 1.004, 1.0045, 1.005, 1.006, 1.007, 1.008, 1.009, 1.01]


println("################################ START ################################ \n")
for q in q_vals
    mps_F= FredkinMPS(N, sites_F, q)
    lognorm_mps_F = []
    normalize!(mps_F; (lognorm!) = lognorm_mps_F)
    @time begin
        println("q = ",q, "\t computation: \n")
        ### OBS: converting to dense mps speeds up calculation but uses  about 5 times as much memory
        Sz_expect_data = expect(dense_mps(mps_F), "Sz")
        @save datadir("sims", savename("Fredkin_Sz", (N=N, q=q), "jld2"; digits = 4)) Sz_expect_data
    end



end
println("################################ END ################################\n ")