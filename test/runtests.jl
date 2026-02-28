using DrWatson, Test
@quickactivate "Fredkin_Motzkin_MPS"

# Here you include files using `srcdir`
using ITensors, ITensorMPS
include(srcdir("Fredkin_MPS.jl"))
include(srcdir("Motzkin_MPS.jl"))
include(srcdir("helpers.jl"))

# Run test suite
println("Starting tests")
ti = time()

@testset "Fredkin_Motzkin_MPS tests" begin
    N = 10
    q = 1.0

    ### Check that the norms are correct (at q = 1 they should be Catalan and Motzkin numbers)
    sites_F = [Index(QN("Sz", 1) => 1, QN("Sz", -1) => 1; tags="S=1/2, Site,n=$i") for i in 1:N]
    mps_F = FredkinMPS(N, sites_F, q)
    @test catalan(div(N, 2)) == norm(mps_F)^(2)


    sites_M = [Index(QN("Sz", 1) => 1, QN("Sz", 0) => 1, QN("Sz", -1) => 1; tags="S=1, Site,n=$i") for i in 1:N]
    mps_M = MotzkinMPS(N, sites_M, q)
    @test motzkin(N) == round(norm(mps_M)^(2), digits = 3)

end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
