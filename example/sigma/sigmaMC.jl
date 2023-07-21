using ElectronLiquid
using CompositeGrids
using JLD2, Printf
# using MPI

# MPI.Init()
dim = 2
rs = [0.5,]
mass2 = [4.0,]
Fs = [-0.0,]
beta = [50.0]
order = [2,]
neval = 4e7
isDynamic = false
isFock = false

# mission = :Z
# mission = :K
mission = ARGS[1]
println("mission: ", mission)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=isDynamic, dim=dim, isFock=isFock)
    kF = para.kF

    if mission == "Z"
        ######### calcualte Z factor ######################
        kgrid = [kF,]
        ngrid = [-1, 0]
        # ngrid = [-1, 0, 1]
    elseif mission == "K"
        ######### calculate K dependence #####################
        Nk, korder = 4, 4
        minK = 0.2kF
        # kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
        # kgrid = kF .+ [-0.1, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1] * kF
        kgrid = kF .+ [-0.1, -0.05, 0, 0.05, 0.1] * kF
        ngrid = [0,]
    else
        error("unknown mission")
    end

    # partition = [(1, 0, 0), (2, 0, 0), (3, 0, 0)]
    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    # @time diagram = Sigma.diagram(para, partition)
    @time diagram = Sigma.diagramGV(para, partition)
    reweight_goal = Float64[]
    for (order, sOrder, vOrder) in partition
        push!(reweight_goal, 4.0^(order + vOrder - 1))
    end
    push!(reweight_goal, 2.0)

    # sigma, result = Sigma.KW(para, diagram;
    sigma, result = Sigma.GV(para, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    if isnothing(sigma) == false
        jldopen("data_$(mission)_GV.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (ngrid, kgrid, sigma)
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
end

