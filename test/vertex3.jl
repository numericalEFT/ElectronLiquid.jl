function compare(data, expect)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=5 * data.err)
end

@testset "Vertex3" begin
    ### test Yukawa interaction ###########
    p = (1, 0, 0)
    mass2 = 0.01
    ####################### STATIC #########################
    para = UEG.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    diagram = Ver3.diagram(para, [p,])
    kin = [[para.kF, 0.0, 0.0],]
    qout = [[0.0, 0.0, 0.0],]
    data, result = Ver3.KW(para, diagram; neval = 1e6, kin=kin, qout=qout, nkin=[0,], nqout=[1, ])
    obs = data[p]
    # for the one-loop vertex3 diagram, we expect
    # \gamma_3(q=0, w->0)  = - dIm\Sigma)/dw_n
    expect = -0.44275
    compare(real(obs[1]), expect)
end