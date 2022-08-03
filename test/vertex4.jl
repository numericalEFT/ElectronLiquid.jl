function compare(data, expect)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=5 * data.err)
end

@testset "Exchange interaction" begin
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=1.0, isDynamic=false)
    kF = para.kF
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=kF^2, isDynamic=false)
    Wp, Wm, θgrid = Ver4.exchange_Coulomb(para)
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    # (1/2)* int dx 1/(4*sin(x/2)^2+1)*sin(x) = log(5)/4
    @test Fp ≈ 4π * para.e0^2 / kF^2 * para.NF * (log(5) / 2) / 2

end

@testset "Vertex4" begin
    ### test Yukawa interaction ###########
    p = (1, 0, 0)
    mass2 = 1.0
    ####################### STATIC #########################
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=false)
    diagram = Ver4.diagram(para, [p,])

    ########## test Yukawa interaction ############
    data, result = Ver4.KW(para, diagram; neval=1e6, print=-1)
    obs = data[p]
    expect = 0.0
    compare(real(obs[1]), expect)
    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)

    ########## test l=0 PH averged Yukawa interaction ############
    data, result = Ver4.PH(para, diagram; neval=1e6, print=-1, l=[0,])
    obs = data[p]
    println(obs)
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_Coulomb(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    # expect = 0.0
    # compare(real(obs[1]), expect)
    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)

    ####################### DYNAMIC #########################
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    diagram = Ver4.diagram(para, [p,])

    ########## test l=0 RPA interaction  (Ω = 0, q->0) ######################
    data, result = Ver4.PH(para, diagram; neval=1e6, print=-1, l=[0,], n=[0, 0, 0])
    obs = data[p]
    println(obs)
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2 + 4π * para.e0^2 * para.NF) * para.NF
    compare(real(obs[2]), expect)


    ########## test l=0 RPA interaction  (Ω -> 0, q=0) ######################
    data, result = Ver4.PH(para, diagram; neval=1e6, print=-1, l=[0,], n=[-1, 0, 0])
    obs = data[p]
    println(obs)
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)

    ########## test l=0 RPA interaction  (Ω -> 0, q = 0) ######################

end