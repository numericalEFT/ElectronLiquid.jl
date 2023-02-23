using ElectronLiquid
using CompositeGrids
using JLD2

rs = [4.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [25.0,]
order = [1,]
neval = 1e6

isDynamic = false
isFock = true

# mission = :Z
mission = :K
# mission = ARGS[1]
# println("mission: ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock)
    kF = para.kF

    ######### calculate K dependence #####################
    Nk, korder = 4, 4
    minK = 0.2kF
    kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder).grid
    tgrid = [_beta - 1e-8,]

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    @time diagram = Green.diagram(para, partition)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    g, result = Green.KT(para, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=kgrid, tgrid=tgrid, neval=neval, paralell=:thread)

    if isnothing(g) == false
        for (ki, k) in enumerate(kgrid)
            println("k = $(k/para.kF), nk = $(g[partition[1]][1, ki])")
        end
        jldopen("data_$(mission).jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, tgrid, kgrid, g)
        end
    end
end

