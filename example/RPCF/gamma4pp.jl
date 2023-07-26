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
using GreenFunc
using JLD2

@inline negative_n(n) = -1 - n
function full_nlist(para; ω_c=0.1para.EF)
    @unpack β, EF = para
    nmax = floor(Int, ω_c / 2 / π * β + 0.5)
    n1grid = GreenFunc.ImFreq(β, FERMION; grid=[i for i in 0:nmax])
    n2grid = GreenFunc.ImFreq(β, FERMION; grid=[i for i in -nmax-1:nmax])
    n = Vector{Int}[]
    for n1 in n1grid.grid
        for n2 in n2grid.grid
            push!(n, [n1, n2, negative_n(n1)])
        end
    end
    return n, n1grid, n2grid
end

function gamma4pp(para, order;
    ω_c=0.1para.EF, neval=1e6, kwargs...)
    @unpack β, EF = para
    ps = [(i, 0, 0) for i in 1:order]
    channel = [PHr, PHEr, PPr]
    filter = [NoHartree, NoBubble]
    diagram = Ver4.diagram(para, ps; channel=channel, filter=filter)

    # nlist = full_nlist(nmax)
    nlist, n1grid, n2grid = full_nlist(para)
    paras = [Ver4.OneAngleAveraged(para, [para.kF, para.kF], nlist, :PP, 0),]
    gamma4 = MeshArray(n1grid, n2grid, [1, 2], [i for i in 1:order]; dtype=Measurement{Float64})
    data, result = Ver4.one_angle_averaged(paras, diagram; neval=neval, print=-1, kwargs...)

    for (i, n) in enumerate(nlist)
        n1, n2 = n[1], n[2]
        i1, i2 = locate(n1grid, n1), locate(n2grid, n2)
        for p in ps
            obs = data[p]
            gamma4.data[i1, i2, 1, p[1]] = real(obs[1, i, 1])
            gamma4.data[i1, i2, 2, p[1]] = real(obs[2, i, 1])
        end
    end

    return data, result, gamma4
end

function save_results(fname, data, result, gamma4)
    f = jldopen(fname, "w")
    f["data"] = data
    f["result"] = result
    f["gamma4"] = gamma4
    close(f)
end

@testset "PP" begin
    rs = 3.0
    beta = 100
    mass2 = 1e-6
    neval = 1e7
    para = ElectronLiquid.ParaMC(rs=rs, beta=beta, Fs=0.0, order=1, mass2=mass2, isDynamic=true)
    UEG.MCinitialize!(para)
    data, result, gamma4 = gamma4pp(para, 2; neval=neval)
    save_results("./run/gamma4pp.jld2", data, result, gamma4)
    println(gamma4)
    println(gamma4.data)
end