function integrandKT(idx, vars, config)
    varK, varT, ExtTidx, ExtKidx = vars
    para, kgrid, tgrid, maxMomNum, extT_labels = config.userdata[1:5]
    leafStat, leaf_maps, momLoopPool, root = config.userdata[6:9]
    graphfuncs! = config.userdata[10][idx]
    isLayered2D = config.userdata[11]

    varK.data[1, 1] = kgrid[ExtKidx[1]]
    varT.data[2] = tgrid[ExtTidx[1]]
    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])

    leafval = _eval_leafvalues!(idx, varK, varT, para, leafStat, leaf_maps, momLoopPool, isLayered2D)
    graphfuncs!(root, leafval)

    weight = sum(root[i] for i in eachindex(extT_labels[idx]))
    loopNum = config.dof[idx][1]
    return weight / (2π)^(para.dim * loopNum)
end

function measureKT(idx, vars, obs, weight, config)
    t = vars[3][1]
    k = vars[4][1]
    obs[idx][t, k] += weight
end

function KT(para::ParaMC, diagram;
    kgrid=[para.kF,],
    tgrid=[para.β - 1e-8,], # must be (0, β)
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Green.KT"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, FeynGraphs, extT_labels = diagram
    maxMomNum, funcGraphs!, leafStat, leaf_maps, momLoopPool, root = _prepare_parquetad(para, diagram)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    T = MCIntegration.Continuous(0.0, β; offset=2, alpha=alpha, adapt=true)
    T.data[1] = 0.0
    T.data[2] = tgrid[1]
    ExtTidx = MCIntegration.Discrete(1, length(tgrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 2, 1, 1] for p in diagpara] # K, T, ExtTidx, ExtKidx
    obs = [zeros(Float64, length(tgrid), length(kgrid)) for o in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, ExtTidx, ExtKidx),
            dof=dof,
            type=Float64,
            obs=obs,
            userdata=(para, kgrid, tgrid, maxMomNum, extT_labels,
                leafStat, leaf_maps, momLoopPool, root, funcGraphs!, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrandKT; config=config, measure=measureKT, print=print, neval=neval, solver=solver, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end

        if print >= -2
            println(result)
        end

        datadict = _result_dict(partition, result, (avg, std) -> measurement.(avg, std))
        return datadict, result
    else
        return nothing, nothing
    end
end
