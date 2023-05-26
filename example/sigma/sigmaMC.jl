using ElectronLiquid
using CompositeGrids
using JLD2
# using MPI

# MPI.Init()
rs = [1.0,]
mass2 = [1.0,]
Fs = [-0.0,]
beta = [40.0,]
order = [2,]
neval = 1e8

# mission = :Z
# mission = :K
mission = ARGS[1]
println("mission: ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=true)
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
        kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 3kF], [kF,], Nk, minK, korder).grid
        ngrid = [0,]
    else
        error("unknown mission")
    end

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    @time diagram = Sigma.diagram(para, partition)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    sigma, result = Sigma.KW(para, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    if isnothing(sigma) == false 
        jldopen("data_$(mission).jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, ngrid, kgrid, sigma)
        end
    end
end

