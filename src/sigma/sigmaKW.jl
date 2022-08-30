
function integrandKW(config)
    para, diag, extT, kgrid, ngrid = config.para
    diagram = diag[config.curr]
    weight = diagram.node.current
    l = config.var[3][1]
    k = config.var[4][1]
    varK, varT = config.var[1], config.var[2]
    varK.data[1, 1] = kgrid[k]
    wn = ngrid[l]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, extT[config.curr][ri], wn, para.β) for (ri, r) in enumerate(diagram.root))

    loopNum = config.dof[config.curr][1]
    factor = 1.0 / (2π)^(para.dim * loopNum)
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measureKW(config)
    l = config.var[3][1]
    k = config.var[4][1]
    o = config.curr
    config.observable[o, l, k] += config.relativeWeight
end

function KW(para::ParaMC, diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, diag, root, extT = diagram

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), length(ngrid), length(kgrid)) # observable for the Fock diagram 

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            obs=obs,
            para=(para, diag, extT, kgrid, ngrid),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = MCIntegration.sample(config, integrandKW, measureKW; print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        avg, std = result.mean, result.stdev
        if print >= 0
            MCIntegration.summary(result.config)
            println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))
        end

        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        data = Complex.(r, i)
        datadict = Dict{eltype(partition),typeof(data[1, :, :])}()
        for i in 1:length(dof)
            datadict[partition[i]] = data[i, :, :]
        end
        return datadict, result
    else
        return nothing, nothing
    end
end