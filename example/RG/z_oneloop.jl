using ElectronLiquid
using ElectronGas
using CompositeGrids
using FiniteDifferences
using GreenFunc
using MCIntegration
using JLD2

rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
mass2 = [1e-5,]
beta = [100.0,]
order = [1,]
neval = 1e7

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true)
    Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 20 * para.kF], [para.kF,], 4, 0.01 * para.kF, 4)
    ngrid = [0, 1]

    f = jldopen("data_f.jld2", "r")
    key = "$(UEG.short(para))"
    para, Λgrid, fs, us = f[key]

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=para.order,
        mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

    partition = UEG.partition(para.order)
    neighbor = UEG.neighbor(partition)
    @time diagram = Sigma.diagram(para, partition; dR=true)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    sigma, result = Sigma.KW_df(paras, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=Λgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    if isnothing(sigma) == false
        jldopen("data_Z.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (para, ngrid, Λgrid, sigma)
        end
    end
end