using ElectronLiquid
using CompositeGrids
using JLD2
# using MPI

# MPI.Init()
dim = 2
rs = [0.5,]
mass2 = [2.0,]
Fs = [-0.0,]
beta = [50.0]
order = [3,]
neval = 2e9

# mission = :Z
# mission = :K
mission = ARGS[1]
println("mission: ", mission)
# exit(0)

for (_rs, _mass2, _F, _beta, _order) in Iterators.product(rs, mass2, Fs, beta, order)
    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=_F, order=_order, mass2=_mass2, isDynamic=false, dim=dim)
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
        kgrid = CompositeGrid.LogDensedGrid(:uniform, [0.0, 2.2kF], [kF,], Nk, minK, korder).grid
        # kgrid = kF .+ [-0.1, -0.05, -0.03, -0.01, -0.005, -0.001, 0, 0.001, 0.005, 0.01, 0.03, 0.05, 0.1] * kF
        # ngrid = [0,]
        ngrid = [-1, 0]
    else
        error("unknown mission")
    end

    partition = UEG.partition(_order)
    neighbor = UEG.neighbor(partition)
    @time diagram = Sigma.diagram(para, partition)
    # reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]
    reweight_goal = [1.0, 1.0, 1.0, 1.0,
        2.0, 2.0, 2.0, 4.0, 4.0, 8.0, 2.0, 2.0, 2.0,
        4.0, 4.0, 8.0, 4.0, 4.0, 8.0, 8.0, 2.0]

    sigma, result = Sigma.KW(para, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    if isnothing(sigma) == false
        jldopen("data_$(mission)_all1.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, ngrid, kgrid, sigma)
        end
    end
end

