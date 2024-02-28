function integrandKW(idx, var, config)
    para, diag, root, extT, kin, qout, nkin, nqout = config.userdata

    kF, β = para.kF, para.β
    varK, varT = var[1], var[2]
    loopNum = config.dof[idx][1]
    # error(loopNum)
    _kin, _qout = kin[var[3][1]], qout[var[4][1]]
    _nkin, _nqout = nkin[var[5][1]], nqout[var[6][1]]

    # println(kinL, ", ", koutL, ", ", kinR)
    varK.data[:, 1] .= _qout
    varK.data[:, 2] .= _kin

    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    # varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    if !isempty(rootuu)
        wuu = sum(weight[root] * phaseC(varT, extTuu[ri], _nkin, _nqout, β) for (ri, root) in enumerate(rootuu))
    else
        wuu = zero(ComplexF64)
    end
    if !isempty(rootud)
        wud = sum(weight[root] * phaseC(varT, extTud[ri], _nkin, _nqout, β) for (ri, root) in enumerate(rootud))
    else
        wud = zero(ComplexF64)
    end
    # factor = para.NF / (2π)^(para.dim * loopNum)
    factor = 1.0 / (2π)^(para.dim * loopNum)

    # if isdefined(Main, :Infiltrator)
    #     Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
    # end

    return Weight(wuu * factor, wud * factor)
    # return Weight(zero(ComplexF64), wud * para.NF)
end

function measureKW(idx, var, obs, weight, config)
    kin, qout = var[3][1], var[4][1]
    nkin, nqout = var[5][1], var[6][1]
    obs[idx][1, kin, qout, nkin, nqout] += weight.d
    obs[idx][2, kin, qout, nkin, nqout] += weight.e
end

function KW(para::ParaMC, diagram;
    kin=[getK(para.kF, para.dim, 1),],
    nkin=[0,],
    qout=[getK(0.0, para.dim, 1),],
    nqout=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF

    partition, diagpara, diag, root, extT = diagram
    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    # X = MCIntegration.Continuous(-1.0, 1.0) #x=cos(θ)

    Nkin, Nqout = length(kin), length(qout)
    Nwin, Nwqout = length(nkin), length(nqout)

    vKin = MCIntegration.Discrete(1, Nkin, alpha=alpha)
    vQout = MCIntegration.Discrete(1, Nqout, alpha=alpha)

    vWin = MCIntegration.Discrete(1, Nwin, alpha=alpha)
    vWqout = MCIntegration.Discrete(1, Nwqout, alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nkin, Nqout, Nwin, Nwqout) for p in diagpara] # observable for the Fock diagram 

    if isnothing(config)
        config = MCIntegration.Configuration(; var=(K, T, vKin, vQout, vWin, vWqout),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(para, diag, root, extT, kin, qout, nkin, nqout),
            kwargs...)
    end
    result = integrate(integrandKW; config=config, measure=measureKW, solver=:mcmc, neval=neval, print=print, kwargs...)

    if isnothing(result) == false
        # if print >= 0
        #     report(result.config)
        #     report(result; pick=o -> (real(o[1, 1, 1, 1, 1])), name="uu")
        #     report(result; pick=o -> (real(o[2, 1, 1, 1, 1])), name="ud")
        # end

        datadict = Dict{eltype(partition),Any}()
        for k in 1:length(dof)
            avg, std = result.mean[k], result.stdev[k]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[k]] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end

end