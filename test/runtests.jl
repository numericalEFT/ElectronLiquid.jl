using Test
include("../common/counterterm.jl")
using .CounterTerm

@testset "CounterTerm" begin
    mu = Dict()
    mu[(1, 0, 0)] = rand()
    mu[(2, 0, 0)] = rand()
    mu[(1, 1, 0)] = rand()
    mu[(1, 0, 1)] = rand()

    δμ, _ = CounterTerm.derive_onebody_parameter_from_sigma(2, mu)
    @test δμ[1] ≈ -mu[(1, 0, 0)]
    @test δμ[2] ≈ -(mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 1, 0)] * δμ[1])

    z = Dict()
    z[(1, 0, 0)] = rand()
    z[(2, 0, 0)] = rand()
    z[(1, 1, 0)] = rand()
    z[(1, 0, 1)] = rand()

    δμ, δz = CounterTerm.derive_onebody_parameter_from_sigma(2, mu, z)
    @test δμ[1] ≈ -mu[(1, 0, 0)]
    @test δz[1] ≈ z[(1, 0, 0)]
    zR = Dict()
    zR[(2, 0)] = z[(2, 0, 0)] + z[(1, 0, 1)] + z[(1, 0, 0)] * δz[1]
    zR[(1, 1)] = z[(1, 1, 0)]
    δz2 = zR[(2, 0)] + zR[(1, 1)] * δμ[1]

    muR = Dict()
    muR[(2, 0)] = mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 0, 0)] * δz[1]
    muR[(1, 1)] = mu[(1, 1, 0)]
    δmu2 = -(muR[(2, 0)] + muR[(1, 1)] * δμ[1])
    @test δμ[2] ≈ δmu2

end