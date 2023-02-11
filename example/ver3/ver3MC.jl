using ElectronLiquid
using CompositeGrids
using MCIntegration
using FeynmanDiagram
using Printf
using JLD2

rs = [5.0,]
mass2 = [0.01, ]
Fs = [-0.0, ]
beta = [25,]
order = [2,]
neval = 1e6

mission = ARGS[1]
println("mission (Z or A): ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
    kF, dim = para.kF, para.dim

    partition = UEG.partition(_order)
    # partition = [(1, 0, 0), (1, 1, 0), (2, 0, 0)]
    neighbor = UEG.neighbor(partition)
    filter = [NoHartree,
        # NoBubble,
        Proper
    ]
    start = time()
    diagram = Ver3.diagram(para, partition; filter=filter)
    println("diagram generation takes: ", time()-start)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    if mission == "Z" || mission=="A"
        ######### calcualte Z factor ######################
        Nk, korder = 4, 4
        minK = 0.2kF
        # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.1 * kF, 2kF], [kF,], Nk, minK, korder).grid
        # kgrid = [0.5*kF, kF, 1.5*kF]
        kgrid = [kF, ]
        qgrid = [0.0, ]
        nin = [0, ]
        nqout = (mission == "A") ? [0,] : [1, ]

        ver3, result = Ver3.KW(para, diagram; neighbor=neighbor,
            kin=[getK(k, dim, 1) for k in kgrid], qout =[getK(q, dim, 1) for q in qgrid],
            nin=nin, nqout = nqout,
            neval=neval, print=0)

        # MCIntegration.summary(result)

        if isnothing(ver3) == false
            for (p, data) in ver3
                printstyled("permutation: $p\n", color=:yellow)
                for (qi, q) in enumerate(qgrid)
                    printstyled("q = $q\n", color=:green)
                    @printf("%12s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "sum")
                    for (ki, k) in enumerate(kgrid)
                        d1, d2 = real(data[1, ki, qi, 1, 1]), real(data[2, ki, qi, 1, 1])
                        s = d1 + d2
                        @printf("%12.6f    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$s")
                    end
                end
            end

            jldopen("ver3_$(mission).jld2", "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (para, kgrid, qgrid, nin, nqout, ver3)
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

