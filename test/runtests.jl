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
    @assert δμ[1] ≈ -mu[(1, 0, 0)]
    @assert δμ[2] ≈ -(mu[(2, 0, 0)] + mu[(1, 0, 0)] + mu[(1, 0, 1)] * δμ[1])
end