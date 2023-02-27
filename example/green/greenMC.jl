using ElectronLiquid
using CompositeGrids
using JLD2

rs = [1.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [40.0,]
order = [3,]
neval = 1e8

isDynamic = false
isFock = false

# mission = :Z
# mission = :K
mission = :n
# mission = ARGS[1]
# println("mission: ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, isFock=isFock)
    kF = para.kF

    ######### calculate K/T dependence #####################
    Nk, korder = 4, 4
    minK = 0.2kF
    kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder).grid
    tgrid = [_beta - 1e-8,]

    lowest_loop_order = 0  # Green's functions start at zeroth loop order
    partition = UEG.partition(_order; offset=lowest_loop_order)
    neighbor = UEG.neighbor(partition)
    @time diagram = Green.diagram(para, partition)
    valid_partition = diagram[1]

    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]
    reweight_pad = repeat([2.0], max(0, length(valid_partition) - length(reweight_goal) + 1))
    reweight_goal = [reweight_goal; reweight_pad]
    @assert length(reweight_goal) ≥ length(valid_partition) + 1

    # Density or green integration (either integrate or bin ExtK)
    if mission == :n
        g, result = Green.densityKT(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal[1:length(valid_partition)+1],
            neval=neval, parallel=:thread)
    else
        g, result = Green.KT(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal[1:length(valid_partition)+1],
            kgrid=kgrid, tgrid=tgrid, neval=neval, parallel=:thread)
    end

    if isnothing(g) == false
        if mission != :n
            for (ki, k) in enumerate(kgrid)
                println("k = $(k/para.kF), nk = $(g[partition[1]][1, ki])")
            end
        end
        jldopen("data_$(mission).jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = mission == :n ? (para, g) : (para, tgrid, kgrid, g)
        end
    end
end

