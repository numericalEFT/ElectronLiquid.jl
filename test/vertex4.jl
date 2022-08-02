function compare(data, expect)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=5 * data.err)
end

@testset "Vertex4" begin
    ### test Yukawa interaction ###########
    p = (1, 0, 0)
    mass2 = 1.0
    para = ElectronLiquid.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=2, mass2=mass2, isDynamic=false)
    diagram = Ver4.diagram(para, [p,])
    data, result = Ver4.KW(para, diagram; neval=1e5, print=-1)
    obs = data[p]
    expect = 0.0
    compare(real(obs[1]), expect)
    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)

    ########## test l=0 PH averged Yukawa interaction ############
    data, result = Ver4.PH(para, diagram; neval=1e5, print=-1, l=[0,])
    obs = data[p]
    exchange = obs[1] - obs[2] # exchange = upup - updn
    Wp, Wm, θgrid = Ver4.exchange_Coulomb(para) # Wp = exchanged Coulomb interaction, Wm = 0
    Fp = Ver4.Legrendre(0, Wp, θgrid)
    compare(real(exchange), Fp)

    # expect = 0.0
    # compare(real(obs[1]), expect)
    expect = -4π * para.e0^2 / (mass2) * para.NF
    compare(real(obs[2]), expect)

end