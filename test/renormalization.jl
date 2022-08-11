using TaylorSeries

@testset "Renormalization" begin
    Renorm.set_params(order=2, numct=2)
    x1, x2, x3 = get_variables()

    keylist = [(1, 0, 0), (2, 0, 0), (1, 1, 0), (1, 0, 1)]

    ########################################################################
    ############################ check number ############################
    ####################################################################
    data = Dict{Any,Any}()
    size = [2, 3]
    for k in keylist
        data[k] = rand()
    end

    ######### check dict to series #############
    series = Renorm.dict2series(data)
    for k in keylist
        @test data[k] ≈ getcoeff(series, k)
    end
    # println(data)
    # println(series)

    ####### check merge ##################
    mseries = Renorm.merge(series, 3; target=1)
    # println(mseries)
    @test data[(2, 0, 0)] + data[(1, 0, 1)] ≈ getcoeff(mseries, (2, 0, 0))
    # println(mseries)

    ######## check renormalization #########
    a, b = rand(), rand()
    rseries = Renorm.renormalize(mseries, 2, a * x1 + b * x1^2)
    @test data[(2, 0, 0)] + data[(1, 0, 1)] + data[(1, 1, 0)] * a ≈ getcoeff(rseries, (2, 0, 0))
    @test data[(1, 0, 0)] ≈ getcoeff(rseries, (1, 0, 0))
    # x1, x2, x3 = get_variables()
    # println(rseries)

    #####################################################################
    ############################# check array first ######################
    ####################################################################
    data = Dict{Any,Any}()
    size = [2, 3]
    for k in keylist
        data[k] = rand(Float64, size...)
    end

    ######### check dict to series #############
    series = Renorm.dict2series(data)
    for k in keylist
        @test data[k][1] ≈ getcoeff(series[1], k)
    end

    ####### check merge ##################
    mseries = Renorm.merge(series, 3; target=1)
    @test data[(2, 0, 0)][1] + data[(1, 0, 1)][1] ≈ getcoeff(mseries[1], (2, 0, 0))

    ######## check renormalization #########
    rseries = Renorm.renormalize(mseries, 2, a * x1 + b * x1^2)
    @test (data[(2, 0, 0)]+data[(1, 0, 1)]+data[(1, 1, 0)]*a)[1] ≈ getcoeff(rseries[1], (2, 0, 0))
    @test data[(1, 0, 0)][1] ≈ getcoeff(rseries[1], (1, 0, 0))




end