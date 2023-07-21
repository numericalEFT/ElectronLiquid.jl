# compute full gamma4 averaged in pp channel

using Test
using ElectronLiquid
using FeynmanDiagram
using FiniteDifferences
using Lehmann
using Measurements
using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG
using ElectronLiquid.ElectronGas.Parameters

@inline negative_n(n) = -1 - n
function full_nlist(nmax)
    n = []
    for i in 1:(nmax+1)^2*2
        n1, n2 = floor(Int, (i - 1) / (2nmax + 2)), (i - 1) % (2nmax + 2) - nmax - 1
        push!(n, [n1, n2, negative_n(n1)])
    end
    return n
end

function gamma4pp(para, order;
    ω_c=0.1para.EF, neval=1e6, kwargs...)
    @unpack β, EF = para
    ps = [(i, 0, 0) for i in 1:order]
    channel = [PHr, PHEr, PPr]
    filter = [NoHartree, NoBubble]
    diagram = Ver4.diagram(para, ps; channel=channel, filter=filter)
    nmax = floor(Int, ω_c / 2 / π * β + 0.5)
    nlist = full_nlist(nmax)
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], nlist, :PP, 0),]
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=neval, print=-1, kwargs...)
    return data, result
end

@testset "PP" begin
    rs = 5.0
    beta = 25
    mass2 = 1e-2
    neval = 1e6
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para)
    data, result = gamma4pp(para, 2)
    println(data)
end