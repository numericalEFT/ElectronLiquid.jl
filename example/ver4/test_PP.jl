using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG


function compare(data, expect, ratio=5)
    # println(data, ", ", expect)
    @test isapprox(data.val, expect, atol=ratio * data.err)
end

function PP_interaction_dynamic(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KO_W(q, n, para)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
end

function yukawa_pp(para)
    factor = para.e0^2 / para.ϵ0 * para.NF / 2
    kF2 = para.kF^2
    return factor / 2 / kF2 * log((para.mass2 + 4 * kF2) / para.mass2)
end

@testset "PP" begin
    seed = 1234
    # p = (1, 0, 0)
    order = 2
    # p = (order, 0, 0)
    # p = (1, 0, 0)
    rs = 1.0
    beta = 25
    mass2 = 1.0
    neval = 1e6
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=order, mass2=mass2, isDynamic=true)
    partition = UEG.partition(para.order)
    println(partition)
    # partition = [(2, 0, 0), (1, 0, 1), (1, 0, 0)]
    UEG.MCinitialize!(para)
    println(para)

    ############################ generic PH one-angle average ###########################
    # nlist = [0, 1, 2]
    # paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[0, nlist[1], -1], [0, nlist[2], -1], [0, nlist[3], -1]], :PP, 0),]
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[0, 0, -1],], :PP, 0),]

    # diagram = Ver4.diagram(para, [p, (2, 0, 0)]; channel=[PPr,], filter=[NoFock, NoBubble])
    # diagram = Ver4.diagramParquet(para, [p,]; channel=[PHr, PHEr, PPr,], filter=[NoFock,])
    # data, result = Ver4.one_angle_averaged_ParquetAD(paras, diagram; neval=neval, print=-1, seed=seed)

    # diagram = Ver4.diagram(para, [p,]; channel=[PHr, PHEr, PPr,], filter=[NoHartree, NoBubble])
    # diagram = Ver4.diagram(para, partition; channel=[PHr, PHEr, PPr,], filter=[NoHartree, NoBubble])
    # data, result = Ver4.one_angle_averaged(paras, diagram; neval=neval, print=-1, seed=seed)
    # println(data[(1, 0, 0)], data[(2, 0, 0)])
    # obs2 = data[p]
    # println(obs2)

    # println(yukawa_pp(para))
    # println("obs 1:", obs[:, 1, 1])
    # println("obs 2:", obs[:, 2, 1])
    # println("obs 3:", obs[:, 3, 1])

    # println(PP_interaction_dynamic(nlist[1], para) / 2)
    # println(PP_interaction_dynamic(nlist[2], para) / 2)
    # println(PP_interaction_dynamic(nlist[3], para) / 2)
    # for i in 1:length(nlist)
    #     println(real(obs[:, i, 1][2]), ", ", -PP_interaction_dynamic(nlist[i], para) / 2)
    #     # compare(real(obs[:, i, 1][2]), -PP_interaction_dynamic(nlist[i], para) / 2)
    # end

    # diagram = Ver4.diagramParquet(para, partition; channel=[PHr, PHEr, PPr,], filter=[NoHartree, NoBubble])
    # data, result = Ver4.one_angle_averaged_ParquetAD(paras, diagram; neval=neval, print=-1, seed=seed)

    diagram = Ver4.diagramParquet_load(para, partition; filter=[NoHartree, NoBubble])
    data, result = Ver4.one_angle_averaged_ParquetAD_Clib(paras, diagram; neval=neval, print=-1, seed=seed)
    println(data[(1, 0, 0)], data[(2, 0, 0)])
    # obs2 = data[p]
    # println(obs2)
end