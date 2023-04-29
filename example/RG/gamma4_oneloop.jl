using ElectronLiquid
using ElectronGas
using CompositeGrids
using FiniteDifferences
using GreenFunc
using MCIntegration
using FeynmanDiagram
using JLD2
using Printf

# rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
rs = [4.0,]
mass2 = [1e-5,]
beta = [100.0,]
order = [1,]
neval = 1e6

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true, isFock=false)
    Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)
    ngrid = [0, 1]

    f = jldopen("data_f.jld2", "r")
    key = "$(UEG.short(para))"
    para, Λgrid, fs, us = f[key]

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=para.order,
        mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

    partition = [(2, 0, 0),]
    channel = [
        PHr,
        PHEr,
        PPr
    ]

    neighbor = UEG.neighbor(partition)
    filter = [NoHartree,
        NoBubble,
        Proper
    ]

    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter, dR=true)

    # reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    n = [-1, 0, 0, -1] # q=0 and w -> 0 to measure F
    # n = [0, 0, 0, 0] # q -> 0 and w = 0 to measure A
    lgrid = [0,] # angular momentum

    # sigma, result = Sigma.KW_df(paras, diagram;
    #     neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
    #     kgrid=Λgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    ver4, result = Ver4.PH_df(paras, diagram;
        kamp=Λgrid, n=n, l=lgrid,
        neighbor=neighbor,
        neval=neval)

    if isnothing(ver4) == false
        for (p, data) in ver4
            printstyled("permutation: $p\n", color=:yellow)
            for (li, l) in enumerate(lgrid)
                printstyled("l = $l\n", color=:green)
                @printf("%12s    %16s    %16s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "di", "ex", "symmetric", "asymmetric")
                for (ki, k) in enumerate(kgrid)
                    factor = 1.0
                    d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
                    s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
                    di, ex = (s - a), (a) * 2.0
                    @printf("%12.6f    %16s    %16s    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$di", "$ex", "$s", "$a")
                end
            end
        end

        jldopen("ver4_PH.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (Λgrid, lgrid, ver4)
        end
    end
end