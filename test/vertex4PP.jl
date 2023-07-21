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

    # θgrid = CompositeGrid.LogDensedGrid(:gauss, [0.0, π], [0.0, π], 16, 0.001, 16)
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    # qs = [2 * kamp * sin(θ / 2) for θ in θgrid.grid]
    # qs = [sqrt(kamp^2 + kamp2^2 - 2 * cos(θ) * kamp * kamp2) for θ in θgrid.grid]
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    # Wm = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        # Wp[qi] = KOstatic(q, para; ct=ct)
        Wp[qi] = UEG.KO_W(q, n, para)
        # Wp[qi] = UEG.KO_W(q, n, para) * sin(θgrid[qi])
        # Wm[qi] = UEG.KOstatic_spin(q, para; ct=ct)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    # Wm *= para.NFstar
    # return Wp, Wm, θgrid
    # return Interp.integrate1D(Wp, θgrid)
    return Interp.integrate1D(Wp, xgrid)
end

@testset "PP" begin
    p = (1, 0, 0)
    rs = 1.0
    beta = 25
    mass2 = 0.1
    neval = 4e6
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para)
    println(para)
    diagram = Ver4.diagram(para, [p,]; channel=[], filter=[])

    ############################ generic PH one-angle average ###########################
    nlist = [0, 5, 10]
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], [[0, nlist[1], -1], [0, nlist[2], -1], [0, nlist[3], -1]], :PP, 0),]
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=neval, print=0)
    obs = data[p]
    println("obs 1:", obs[:, 1, 1])
    println("obs 2:", obs[:, 2, 1])
    println("obs 3:", obs[:, 3, 1])

    println(PP_interaction_dynamic(nlist[1], para) / 2)
    println(PP_interaction_dynamic(nlist[2], para) / 2)
    println(PP_interaction_dynamic(nlist[3], para) / 2)

end