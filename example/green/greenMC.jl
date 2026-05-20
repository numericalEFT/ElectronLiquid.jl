using ElectronLiquid
using CompositeGrids
using JLD2

rs = [1.0,]
mass2 = [0.01,]
Fs = [-0.0,]
beta = [25.0,]
order = [3,]
neval = 1e6

# isDynamic = true
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
    Nk, korder = 5, 5
    minK = 0.02kF
    kgrid = CompositeGrid.LogDensedGrid(:cheb, [0.0, 3kF], [kF,], Nk, minK, korder).grid
    tgrid = [para.β - 1e-8,]

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
    elseif mission == :K
        g, result = Green.KT(para, diagram;
            neighbor=neighbor, reweight_goal=reweight_goal[1:length(valid_partition)+1],
            kgrid=kgrid, tgrid=tgrid, neval=neval, parallel=:thread)
    else
        error("unknown mission: $(mission)")
    end

    if isnothing(g) == false
        if mission == :K
            for (ki, k) in enumerate(kgrid)
                println("k = $(k/para.kF), $(g[(0, 1, 0)][1, ki]), $(g[(1, 0, 0)][1, ki])")
            end
        elseif mission == :n
            for (ip, p) in enumerate(valid_partition)
                println("n$(p) = $(g[p]/g[(0, 0, 0)])")
            end
        else
            error("unknown mission: $(mission)")
        end
        jldopen("data_$(mission).jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = mission == :n ? (g,) : (tgrid, kgrid, g)
        end
    end
end

