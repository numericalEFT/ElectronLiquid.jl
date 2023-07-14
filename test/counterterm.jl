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
