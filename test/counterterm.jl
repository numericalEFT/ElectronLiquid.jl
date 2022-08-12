@testset "CounterTerm" begin
    mu = Dict()
    mu[(1, 0, 0)] = rand()
    mu[(2, 0, 0)] = rand()
    mu[(1, 1, 0)] = rand()
    mu[(1, 0, 1)] = rand()

    sw = Dict()
    sw[(1, 0, 0)] = rand()
    sw[(2, 0, 0)] = rand()
    sw[(1, 1, 0)] = rand()
    sw[(1, 0, 1)] = rand()

    δzi, δμ, δz = CounterTerm.sigmaCT(2, mu, sw)
    @test δμ[1] ≈ -mu[(1, 0, 0)]
    @test δμ[2] ≈ -(mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 1, 0)] * δμ[1])

    rmu = CounterTerm.chemicalpotential_renormalization(2, mu, δμ)
    @test δμ ≈ rmu .* (-1)

    @test δzi[1] ≈ sw[(1, 0, 0)]
    @test δzi[2] ≈ (sw[(2, 0, 0)] + sw[(1, 0, 1)] + sw[(1, 1, 0)] * δμ[1])

    ######## 1+δz[1]*x+δz[2]*x^2+ ... ≈ 1/(1+δzi[1]*x+δzi[2]*x^2+...) ########
    @test δz[1] ≈ -δzi[1]
    @test δz[2] ≈ (δzi[1]^2 - δzi[2])
    # @test δz[3] ≈ (-δzi[1]^3 + δzi[1] * δzi[2] - δzi[3])

    # data = CounterTerm.renormalization(2, sw, δμ, δz; nbody=1, zrenorm=true)
    swR = CounterTerm.chemicalpotential_renormalization(2, sw, δμ)
    swR = [1.0, swR...]
    # println(swR)
    # println(δzi)
    swRR = CounterTerm.z_renormalization(2, swR, δz, 1)
    # expect swRR = [1.0, 0.0, 0.0, ...]
    @test swRR[1] ≈ swR[1]
    for o in 2:2
        @test swRR[o] ≈ 0.0
    end


    # δμ, δz = CounterTerm.sigmaCT(2, mu, z)
    # @test δμ[1] ≈ -mu[(1, 0, 0)]
    # @test δz[1] ≈ -z[(1, 0, 0)]
    # zR = Dict()
    # zR[(2, 0)] = z[(2, 0, 0)] + z[(1, 0, 1)] + z[(1, 0, 0)] * δz[1]
    # zR[(1, 1)] = z[(1, 1, 0)]
    # δz2 = zR[(2, 0)] + zR[(1, 1)] * δμ[1]
    # @test δz[2] ≈ δz2

    # muR = Dict()
    # muR[(2, 0)] = mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 0, 0)] * δz[1]
    # muR[(1, 1)] = mu[(1, 1, 0)]
    # δmu2 = -(muR[(2, 0)] + muR[(1, 1)] * δμ[1])
    # @test δμ[2] ≈ δmu2

end