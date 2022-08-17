using ElectronLiquid
using CompositeGrids
using MCIntegration
using FeynmanDiagram
using Printf
using JLD2

rs = [5.0,]
mass2 = [0.01, 0.001]
Fs = [-1.0,]
beta = [25,]
order = [2,]
neval = 1e9

mission = ARGS[1]
println("mission (L or K): ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
    kF = para.kF

    # partition = UEG.partition(_order)
    partition = [(2, 0, 0),]
    channel = [
        PHr,
        PHEr,
        PPr
    ]
    neighbor = UEG.neighbor(partition)
    filter = [NoHatree,
        NoBubble,
        Proper
    ]
    start = time()
    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter)
    println("diagram generation takes: ", time() - start)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    if mission == "L"
        ######### calcualte Z factor ######################
        Nk, korder = 4, 4
        minK = 0.2kF
        kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.1 * kF, 2kF], [kF,], Nk, minK, korder).grid
        # kgrid = [0.5 * kF, kF, 1.5 * kF]
        # kgrid = [kF, ]
        # n = [0, 1, 1, 0] # q=0 and w -> 0 to measure F
        n = [-1, 0, 0, -1] # q=0 and w -> 0 to measure F
        # n = [0, 0, 0, 0] # q -> 0 and w = 0 to measure A
        lgrid = [0, 1] # angular momentum

        ver4, result = Ver4.PH(para, diagram;
            kamp=kgrid, n=n, l=lgrid,
            neval=neval, print=0)

        # MCIntegration.summary(result)

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

            jldopen("ver4_$(mission).jld2", "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (para, kgrid, lgrid, ver4)
            end
        end

        # elseif mission == "K"
        #     ######### calculate K dependence #####################
        #     Nk, korder = 4, 4
        #     minK = 0.2kF
        #     kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder).grid
        #     ngrid = [0,]

        #     ver4, result = Ver4.PH(para, diagram;
        #         kamp=kamp, n=[-1, 0, 1],
        #         neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        #         kgrid=kgrid, ngrid=ngrid, neval=neval)
    else
        error("unknown mission")
    end
end

