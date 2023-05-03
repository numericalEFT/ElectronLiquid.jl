using ElectronLiquid
using ElectronGas
using CompositeGrids
using FiniteDifferences
using GreenFunc
using MCIntegration
using JLD2

include("gamma4_treelevel.jl")

# rs = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
rs = [4.0,]
mass2 = [1e-5,]
beta = [50.0,]
order = [1,]
neval = 1e7

scheme = :KO
# scheme = :CS

for (_rs, _mass2, _beta, _order) in Iterators.product(rs, mass2, beta, order)

    para = UEG.ParaMC(rs=_rs, beta=_beta, Fs=-0.0, order=_order, mass2=_mass2, isDynamic=true, isFock=false)
    Λgrid = CompositeGrid.LogDensedGrid(:gauss, [1.0 * para.kF, 100 * para.kF], [para.kF,], 8, 0.01 * para.kF, 8)
    ngrid = [-1, 0, 1]

    if scheme == :CS
        f = jldopen("data_f.jld2", "r")
        key = "$(UEG.short(para))"
        para, Λgrid, fs, us = f[key]
    elseif scheme == :KO
        fs = [KO(para, lambda, lambda; verbose=0)[1] for (li, lambda) in enumerate(Λgrid)]
    else
        error("scheme not implemented")
    end
    println(fs)

    paras = [UEG.ParaMC(rs=para.rs, beta=para.beta, Fs=fs[li], Fa=0.0, order=para.order,
        mass2=para.mass2, isDynamic=true, isFock=false) for li in eachindex(Λgrid)]

    partition = UEG.partition(para.order)
    neighbor = UEG.neighbor(partition)
    diagram = Sigma.diagram(para, partition; dR=true)
    reweight_goal = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]

    # sigma, result = Sigma.KW_df(paras, diagram;
    #     neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
    #     kgrid=Λgrid, ngrid=ngrid, neval=neval, parallel=:thread)

    sigma_df, result = Sigma.KW_df(paras, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        kgrid=Λgrid, ngrid=ngrid, neval=neval)

    plist = Array{Tuple{ParaMC,Float64,Int}}(undef, (3, length(Λgrid)))
    for ip in eachindex(paras)
        plist[1, ip] = (paras[ip], Λgrid[ip], -1)
        plist[2, ip] = (paras[ip], Λgrid[ip], 0)
        plist[3, ip] = (paras[ip], Λgrid[ip], 1)
    end

    diagram = Sigma.diagram(para, partition; dR=false)
    sigma, result = Sigma.Generic(plist, diagram;
        neighbor=neighbor, reweight_goal=reweight_goal[1:length(partition)+1],
        neval=neval)

    if isnothing(sigma) == false
        jldopen("data_Z_$scheme.jld2", "a+") do f
            key = "$(UEG.short(para))"
            if haskey(f, key)
                @warn("replacing existing data for $key")
                delete!(f, key)
            end
            f[key] = (ngrid, Λgrid, sigma, sigma_df)
        end
    end
end