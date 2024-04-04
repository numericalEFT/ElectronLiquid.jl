seed = 1234

@testset "Exchange interaction" begin
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=1.0, isDynamic=false)
    kF = para.kF
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=kF^2, isDynamic=false)
    Wp, Wm, θgrid = Ver4.exchange_Coulomb(para)
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    # (1/2)* int dx 1/(4*sin(x/2)^2+1)*sin(x) = log(5)/4
    @test Fp ≈ 4π * para.e0^2 / kF^2 * para.NF * (log(5) / 2) / 2

end

@testset "Vertex4 Tree-level" begin
    ### test Yukawa interaction ###########
    p = (0, 0, 0)
    mass2 = 1.0
    ####################### STATIC #########################
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=false)
    # diagram = Ver4.diagram(para, [p,])
    diagram = Diagram.diagram_parquet_response(:vertex4, para, [p,])

    ########## test Yukawa interaction ############
    # data, result = Ver4.KW(para, diagram; neval=1e5, print=-1, seed=seed)
    # obs = data[p]
    # expect = 0.0
    # compare(real(obs[1]), expect)
    # expect = -4π * para.e0^2 / (mass2) * para.NF
    # compare(real(obs[2]), expect)

    ########## test l=0 PH averged Yukawa interaction ############
    data, result = Ver4.lavg(para, diagram; neval=1e5, print=-1, l=[0,], seed=seed)
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
    # diagram = Ver4.diagram(para, [p,])
    diagram = Diagram.diagram_parquet_response(:vertex4, para, [p,])

    ########## test l=0 RPA interaction  (Ω = 0, q->0) ######################
    data, result = Ver4.lavg(para, diagram; neval=1e5, print=-1, l=[0,], n=[0, 0, 0], seed=seed)
    obs = data[p]
    println(obs)
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2 + 4π * para.e0^2 * para.NF) * para.NF
    compare(real(obs[2]), expect)


    ########## test l=0 RPA interaction  (Ω -> 0, q=0) ######################
    data, result = Ver4.lavg(para, diagram; neval=1e5, print=-1, l=[0,], n=[-1, 0, 0], seed=seed)
    obs = data[p]
    println(obs)
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)


    ############################ generic PH one-angle average ###########################
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[-1, 0, 0], [0, 0, 0]], :PH, 0),]
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=1e5, print=-1, seed=seed)
    obs = data[p]
    println("obs 1:", obs[:, 1, 1])
    println("obs 2:", obs[:, 2, 1])

    #############   (Ω -> 0, q=0)   #####################
    exchange = obs[1, 1, 1] - obs[2, 1, 1] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2, 1, 1]), expect)

    #############   (Ω = 0, q->0)   #####################
    exchange = obs[1, 2, 1] - obs[2, 2, 1] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_interaction(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    expect = -4π * para.e0^2 / (mass2 + 4π * para.e0^2 * para.NF) * para.NF
    compare(real(obs[2, 2, 1]), expect)

    include("vertex4PP.jl")
end

@testset "Vertex4 One-loop" begin
    ########## test l=0 one-loop interaction  (Ω -> 0, q = 0) ######################
    p = (1, 0, 0)
    mass2 = 0.01
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    diagram = Ver4.diagram(para, [p,]; filter=[Proper, NoBubble], channel=[PHr, PHEr, PPr])
    data, result = Ver4.lavg(para, diagram; neval=1e6, print=0, l=[0,], n=[-1, 0, 0], seed=seed)
    obs = data[p]

    expect = 1.0314 # +- 0.0044
    compare(real(obs[1]), expect)

    expect = 0.4723 # +- 0.0041
    compare(real(obs[2]), expect)

    ############################ generic PH one-angle average ###########################
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[-1, 0, 0],], :PH, 0),]
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=1e5, print=-1, seed=seed)
    obs = data[p]

    expect = 1.0314 # +- 0.0044
    compare(real(obs[1, 1, 1]), expect)

    expect = 0.4723 # +- 0.0041
    compare(real(obs[2, 1, 1]), expect)
end