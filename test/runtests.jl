using DrWatson, Test
@quickactivate "Fredkin_Motzkin_MPS"

# Here you include files using `srcdir`
using ITensors, ITensorMPS
include(srcdir("Fredkin_MPS.jl"))
include(srcdir("Motzkin_MPS.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Fredkin_Motzkin_MPS tests" begin
    N = 10
    q = 1.0
    sites_F = [Index(QN("Sz", 1) => 1, QN("Sz", -1) => 1; tags="S=1/2, Site,n=$i") for i in 1:N]
    mps_F = FredkinMPS(N, sites_F, q)
    #println(expect(mps_F, "Sz"))


    sites_M = [Index(QN("Sz", 1) => 1, QN("Sz", 0) => 1, QN("Sz", -1) => 1; tags="S=1, Site,n=$i") for i in 1:N]
    mps_M = MotzkinMPS(N, sites_M, q)
    println(expect(mps_M, "Sz"))

    # MPS_ket = MPS(mps.tensors)
    # # --- Compute expectation values of Sz on all sites ---
    # exp_vals = expect(MPS_ket, "Sz"; sites=1:N)

    # # --- Plot ---
    # plot(1:N, exp_vals,
    #      xlabel="Site index",
    #      ylabel="⟨Sz⟩",
    #      ylims = (-1.5, 1.5),
    #      title="Expectation value of Sz across the chain",
    #      lw=2, marker=:circle)
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
