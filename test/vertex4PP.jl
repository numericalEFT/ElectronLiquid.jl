# using Test
# using ElectronLiquid
# using FeynmanDiagram
# using FiniteDifferences
# using Lehmann
# using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG


# function compare(data, expect, ratio=5)
#     # println(data, ", ", expect)
#     @test isapprox(data.val, expect, atol=ratio * data.err)
# end

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

@testset "PP" begin
    seed = 1234
    p = (0, 0, 0)
    rs = 5.0
    beta = 25
    mass2 = 1e-2
    neval = 1e6
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para)
    println(para)
    # diagram = Ver4.diagram(para, [p,]; channel=[], filter=[])
    diagram = Diagram.diagram_parquet_response(:vertex4, para, [p,]; channels=[], filter=[])

    ############################ generic PH one-angle average ###########################
    nlist = [0, 1, 2]
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[0, nlist[1], -1], [0, nlist[2], -1], [0, nlist[3], -1]], :PP, 0),]
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=neval, print=-1, seed=seed)
    # data, result = Ver4.one_angle_averaged_ParquetAD(paras, diagram; neval=neval, print=-1, seed=seed)
    obs = data[p]
    # println("obs 1:", obs[:, 1, 1])
    # println("obs 2:", obs[:, 2, 1])
    # println("obs 3:", obs[:, 3, 1])

    # println(PP_interaction_dynamic(nlist[1], para) / 2)
    # println(PP_interaction_dynamic(nlist[2], para) / 2)
    # println(PP_interaction_dynamic(nlist[3], para) / 2)
    for (i, n) in enumerate(nlist)
        compare(real(obs[:, i, 1][2]), -PP_interaction_dynamic(n, para) / 2)
    end
end