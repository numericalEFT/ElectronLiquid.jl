using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG

using ElectronLiquid.ElectronGas
using ElectronLiquid.ElectronGas.GreenFunc

function DCkernel(param, nlist)
    kF = param.kF
    Euv, rtol = 100 * param.EF, 1e-10
    Nk, minK, order, maxK = 8, 1e-7kF, 8, 10kF
    int_type = :rpa

    #--- prepare kernel ---
    W = LegendreInteraction.DCKernel0(param;
        Euv=Euv, rtol=rtol, Nk=Nk, maxK=maxK, minK=minK, order=order,
        int_type=int_type)
    kgrid = W.kgrid
    qgrids = W.qgrids
    ik = locate(kgrid, kF)
    iq = locate(qgrids[ik], kF)
    Wfreq = W.kernel[ik, iq, :] .+ W.kernel_bare[ik, iq]
    dlr = W.dlrGrid
    Wdlr = Lehmann.matfreq2dlr(dlr, Wfreq)
    Wn = Lehmann.dlr2matfreq(dlr, Wdlr, nlist)
    return Wn ./ kF^2
end

function KO_Ws(q, n::Integer, para::ParaMC; Pi=polarKW(q, n, para))
    if abs(q) < 1e-6
        q = 1e-6
    end
    invKOinstant = 1.0 / KOinstant(q, para)
    Rs = 1.0 ./ (invKOinstant - Pi) - para.fs
    return Rs
end

function PPs(n, para::ParaMC, kamp=para.kF, kamp2=para.kF; kawargs...)
    kF = para.kF
    xgrid = CompositeGrid.LogDensedGrid(:gauss, [-1.0, 1.0], [-1.0, 1.0], 16, 0.001, 16)
    qs = [sqrt(kamp^2 + kamp2^2 - 2 * x * kamp * kamp2) for x in xgrid]

    Wp = zeros(Float64, length(qs))
    for (qi, q) in enumerate(qs)
        Wp[qi] = UEG.KO_Ws(q, n, para)
    end
    Wp *= para.NFstar  # additional minus sign because the interaction is exchanged
    return Interp.integrate1D(Wp, xgrid)
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

@testset "PP" begin
    seed = 1234
    p = (1, 0, 0)
    rs = 3.0
    beta = 100
    mass2 = 1e-2
    neval = 1e6
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para)

    nlist = [0, 1, 2, 3, 10, 100, 1000]
    for i in nlist
        W1 = PP_interaction_dynamic(i, para) / para.NFstar
        println(W1)
    end
    println(DCkernel(para.basic, nlist))
end