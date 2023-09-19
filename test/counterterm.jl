@testset verbose=true "CounterTerm basic" begin
    #### test inverse #################
    z = [0.1, 0.21, 0.321, 0.4321, 0.54321] # 1+0.1x+0.21x^2+0.321x^3+0.4321x^4+0.54321x^5
    zi = CounterTerm._inverse(z)
    expect = [-0.1, -0.2, -0.28, -0.33, -0.344]
    @test zi ≈ expect

    zm = [[zi zi; zi zi] for zi in z]
    zim = CounterTerm._inverse(zm)
    @test [e[1, 1] for e in zim] ≈ expect


    ######## test renormalization ########
    mu = [0.1, 0.21, 0.321, 0.4321, 0.54321]
    inv_mu = [-0.1, -0.2, -0.28, -0.33, -0.344]
    partition = UEG.partition(5)

    z = Dict()
    z[(5, 0)] = 0.5
    z[(4, 0)], z[(4, 1)] = 0.4, 0.41
    z[(3, 0)], z[(3, 1)], z[(3, 2)] = 0.3, 0.31, 0.32
    z[(2, 0)], z[(2, 1)], z[(2, 2)], z[(2, 3)] = 0.2, 0.21, 0.22, 0.23
    z[(1, 0)], z[(1, 1)], z[(1, 2)], z[(1, 3)], z[(1, 4)] = 0.1, 0.11, 0.12, 0.13, 0.14


    mu = [0.1, 0.2, 0.3, 0.4]
    zr = CounterTerm.chemicalpotential_renormalization(5, z, mu)
    expect = [0.1, 0.211, 0.3442, 0.51313, 0.735024]
    @test zr ≈ expect


    for k in keys(z)
        z[k] = [z[k] z[k]; z[k] z[k]]
    end
    zr = CounterTerm.chemicalpotential_renormalization(5, z, mu)
    @test zr ≈ [[e e; e e] for e in expect]


    ######## test z_renormalization #########
    z = Dict()
    z[(5, 0)] = 0.5
    z[(4, 0)], z[(4, 1)] = 0.4, 0.41
    z[(3, 0)], z[(3, 1)], z[(3, 2)] = 0.3, 0.31, 0.32
    z[(2, 0)], z[(2, 1)], z[(2, 2)], z[(2, 3)] = 0.2, 0.21, 0.22, 0.23
    z[(1, 0)], z[(1, 1)], z[(1, 2)], z[(1, 3)], z[(1, 4)] = 0.1, 0.11, 0.12, 0.13, 0.14

    dz = [0.1, 0.21, 0.321, 0.4321, 0.54321] # 1+0.1x+0.21x^2+0.321x^3+0.4321x^4+0.54321x^5
    dzi = [-0.1, -0.2, -0.28, -0.33, -0.344] #1/(1+dz)

    zr = CounterTerm.z_renormalization(5, z, dz, 1) # z*dz
    zrr = CounterTerm.z_renormalization(5, zr, dzi, 1) # z*dz*dzi -> z

    for k in keys(z)
        @test z[k] ≈ zrr[k]
    end

    @testset "I/O" begin
        # Test saving and loading to/from a CSV tempfile
        mktemp() do tmp, _
            # Split tmp into root directory and base names
            tmpdir, tmpname = dirname(tmp), basename(tmp)

    for k in keys(z)
        z[k] = [z[k] z[k]; z[k] z[k]]
    end
    zr = CounterTerm.z_renormalization(5, z, dz, 1) # z*dz
    zrr = CounterTerm.z_renormalization(5, zr, dzi, 1) # z*dz*dzi -> z

    for k in keys(z)
        @test zrr[k] ≈ z[k]
    end

end

@testset verbose = true "CounterTerm" begin
    # Setup counterterms
    order = 2
    partitions = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1)]
    mu = Dict()
    sw = Dict()
    # Fill mu and sw with random measurements
    for P in partitions
        mu[P] = rand() ± rand()
        sw[P] = rand() ± rand()
    end
    δzi, δμ, δz = CounterTerm.sigmaCT(2, mu, sw)

    @testset "Renormalization" begin
        @test δμ[1] ≈ -mu[(1, 0, 0)]
        @test δμ[2] ≈ -(mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 1, 0)] * δμ[1])

        rmu = CounterTerm.chemicalpotential_renormalization(order, mu, δμ)
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

    @testset "I/O" begin
        # Test saving and loading to/from a CSV tempfile
        mktemp() do tmp, _
            # Split tmp into root directory and base names
            tmpdir, tmpname = dirname(tmp), basename(tmp)

            # Test parameters
            para = UEG.ParaMC(rs=1.0, beta=40.0, order=order, mass2=2.0, isDynamic=false)
            paraid = UEG.paraid(para)

            # Test saving to CSV file
            save_successful = true
            try
                df = nothing
                for P in partitions
                    df = CounterTerm.appendDict(
                        df,
                        paraid,
                        Dict(
                            "partition" => P,
                            "μ" => mu[P].val,
                            "μ.err" => mu[P].err,
                            "Σw" => sw[P].val,
                            "Σw.err" => sw[P].err,
                        );
                        replace=true,
                        verbose=false
                    )
                end
                CounterTerm.toFile(df, tmpname; root_dir=tmpdir, verbose=false)
            catch
                save_successful = false
            end
            @test save_successful

            # Test loading from CSV file
            local df
            load_successful = true
            try
                df = CounterTerm.fromFile(tmpname; root_dir=tmpdir, verbose=false)
            catch
                load_successful = false
            end
            @test load_successful

            # Test getSigma using preloaded DataFrame
            μ, Σw = CounterTerm.getSigma(df, paraid, order)
            for P in partitions
                @test μ[P] ≈ mu[P]
                @test Σw[P] ≈ sw[P]
            end

            # Test getSigma using CSV file
            μ, Σw = CounterTerm.getSigma(para; parafile=tmpname, root_dir=tmpdir)
            for P in partitions
                @test μ[P] ≈ mu[P]
                @test Σw[P] ≈ sw[P]
            end
        end
    end
end

@testset verbose = true "CounterTerm with Array" begin
    function compare(a, b)
        @test a[1, 1] ≈ b[1, 1]
    end

    # Setup counterterms
    order = 2
    partitions = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1)]
    mu = Dict()
    sw = Dict()
    # Fill mu and sw with random measurements
    for P in partitions
        # 2x2 matrix
        # mu[P] = rand() ± rand()*0.001
        # sw[P] = rand() ± rand()*0.001
        mu[P] = [rand() ± rand()*0.001 rand() ± rand()*0.001; rand() ± rand()*0.001 rand() ± rand()*0.001]
        sw[P] = [rand() ± rand()*0.001 rand() ± rand()*0.001; rand() ± rand()*0.001 rand() ± rand()*0.001]
        # mu[P] = [rand() ± rand()*0.001,]
        # sw[P] = [rand() ± rand()*0.001,]
    end
    δzi, δμ, δz = CounterTerm.sigmaCT(2, mu, sw)

    @testset "Renormalization" begin
        @test δμ[1] ≈ -mu[(1, 0, 0)]
        # @test compare(δμ[2], -(mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 1, 0)] * δμ[1]))
        @test δμ[2] ≈ -(mu[(2, 0, 0)] + mu[(1, 0, 1)] + mu[(1, 1, 0)] .* δμ[1])

        # exit(0)

        rmu = CounterTerm.chemicalpotential_renormalization(order, mu, δμ)
        @test δμ ≈ rmu .* (-1)

        @test δzi[1] ≈ sw[(1, 0, 0)]
        @test δzi[2] ≈ (sw[(2, 0, 0)] + sw[(1, 0, 1)] + sw[(1, 1, 0)] .* δμ[1])

        ######## 1+δz[1]*x+δz[2]*x^2+ ... ≈ 1/(1+δzi[1]*x+δzi[2]*x^2+...) ########
        @test δz[1] ≈ -δzi[1]
        @test δz[2] ≈ (δzi[1].^2 .- δzi[2])
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
            @test swRR[o] ≈ [0.0 0.0; 0.0 0.0]
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
end
