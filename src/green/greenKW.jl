function integrandKW(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, maxMomNum, extT_labels = config.userdata[1:5]
    leafStat, leaf_maps, momLoopPool, root = config.userdata[6:9]
    graphfuncs! = config.userdata[10][idx]
    isLayered2D = config.userdata[11]

    varK.data[1, 1] = kgrid[ExtKidx[1]]
    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])

    leafval = _eval_leafvalues!(idx, varK, varT, para, leafStat, leaf_maps, momLoopPool, isLayered2D)
    graphfuncs!(root, leafval)

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, para.β) for (i, extT) in enumerate(extT_labels[idx]))
    loopNum = config.dof[idx][1]
    return weight / (2π)^(para.dim * loopNum)
end

function measureKW(idx, vars, obs, weight, config)
    n = vars[3][1]
    k = vars[4][1]
    obs[idx][n, k] += weight
end

function KW(para::ParaMC, diagram;
    kgrid=[para.kF],
    ngrid=[0],
    neval=1e6,
    print=0,
    alpha=3.0,
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Green.KW"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, FeynGraphs, extT_labels = diagram
    maxMomNum, funcGraphs!, leafStat, leaf_maps, momLoopPool, root = _prepare_parquetad(para, diagram)

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    T = MCIntegration.Continuous(0.0, β; offset=1, alpha=alpha, adapt=true)
    T.data[1] = 0.0
    varN = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara]
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, varN, ExtKidx),
            dof=dof,
            type=ComplexF64,
            obs=obs,
            userdata=(para, kgrid, ngrid, maxMomNum, extT_labels,
                leafStat, leaf_maps, momLoopPool, root, funcGraphs!, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrandKW; config=config, measure=measureKW, print=print, neval=neval, solver=solver, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end
        if print >= -2
            println(result)
        end

        datadict = _result_dict(partition, result, (avg, std) ->
            Complex.(measurement.(real(avg), real(std)), measurement.(imag(avg), imag(std)))
        )
        return datadict, result
    else
        return nothing, nothing
    end
end
