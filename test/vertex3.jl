
import FeynmanDiagram.FrontEnds: Filter, NoHartree, NoFock, DirectOnly, Wirreducible, Girreducible, NoBubble, Proper

function compare(data, expect)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=5 * data.err)
end

@testset "Vertex3" begin

    # @testset "Vertex3 init" begin
    #     para = UEG.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=2, mass2=0.01, isDynamic=false)
    #     ver3, result = Ver3.MC_KW_angle(para)
    #     println(ver3[(2, 0, 0)])
    # end

    @testset "Vertex3 Dynamic O(1)" begin
        para = UEG.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=0.01, isDynamic=true)
        kin = [para.kF,]
        # Nth = 4
        # theta = [(i) / (Nth * π) for i in 0:Nth] # N+1 points
        theta = [0.0,]
        qout = [[para.kF * (1 - cos(θ)), -para.kF * sin(θ), 0.0] for θ in theta]
        ver3, result = Ver3.MC_KW_angle(para; kin=kin, qout=qout, nkin=[0,], nqout=[-1,])
        obs = ver3[(1, 0, 0)][:, 1, 1, :, 1]
        println(obs)
        # println(ver3[(2, 0, 0)][:, 1, 1, :, 1])

        # for the one-loop vertex3 diagram, we expect
        # \gamma_3(q=0, w->0)  = - dIm\Sigma)/dw_n
        # using ElectronGas
        # para = Parameter.rydbergUnit(1/25, 5.0, 3, Λs=0.01)
        # sigma=ElectronGas.SelfEnergy.G0W0(para)
        # z=SelfEnergy.zfactor(para, sigma[1]; ngrid=[-1,0])[1]
        # ngrid=[-1, 0] -> -0.508
        # ngrid=[0,1] -> -0.44
        expect = -0.44275
        compare(real(obs[1]), expect)

    end

    # @testset "Vertex3 static" begin
    #     para = UEG.ParaMC(rs=1.0, beta=25.0, Fs=0.0, order=3, mass2=3.5, isDynamic=false)
    #     kin = [para.kF,]
    #     Nth = 4
    #     theta = [(i) / (Nth * π) for i in 0:Nth] # N+1 points
    #     # theta = [0.0,]
    #     qout = [[para.kF * (1 - cos(θ)), -para.kF * sin(θ), 0.0] for θ in theta]
    #     transferLoop = [1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    #     ver3, result = Ver3.MC_KW_angle(para;
    #         kin=kin, qout=qout, nkin=[0,], nqout=[0,],
    #         filter=[NoHartree, Proper], transferLoop=transferLoop)
    #     println(ver3)
    #     # obs = ver3[(1, 0, 0)][:, 1, 1, :, 1]
    #     # println(obs)
    #     println(ver3[(1, 0, 0)][:, 1, 1, :, 1] + ver3[(2, 0, 0)][:, 1, 1, :, 1] + ver3[(3, 0, 0)][:, 1, 1, :, 1])

    # end

end

# @testset "Vertex3" begin
#     ### test Yukawa interaction ###########
#     p = (1, 0, 0)
#     mass2 = 0.01
#     ####################### STATIC #########################
#     para = UEG.ParaMC(rs=5.0, beta=25.0, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
#     diagram = Ver3.diagram(para, [p,])
#     kin = [[para.kF, 0.0, 0.0],]
#     qout = [[0.0, 0.0, 0.0],]
#     data, result = Ver3.KW(para, diagram; neval = 1e6, kin=kin, qout=qout, nkin=[0,], nqout=[1, ])
#     obs = data[p]
#     # for the one-loop vertex3 diagram, we expect
#     # \gamma_3(q=0, w->0)  = - dIm\Sigma)/dw_n
#     expect = -0.44275
#     compare(real(obs[1]), expect)
# end