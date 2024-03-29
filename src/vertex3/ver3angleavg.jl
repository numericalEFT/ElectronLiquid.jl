function integrandAA(idx, var, config)
    para, diag, root, extT, kamp, kamp2, nkin, nqout = config.userdata

    kF, β = para.kF, para.β
    varK, varT, varX = var[1], var[2], var[7]
    loopNum = config.dof[idx][1]
    # error(loopNum)
    _kin, _kout = kamp[var[3][1]], kamp2[var[4][1]]
    _nkin, _nqout = nkin[var[5][1]], nqout[var[6][1]]

    # println(kinL, ", ", koutL, ", ", kinR)
    x = varX[1][1]
    varK.data[1, 1] = _kout * x - _kin
    varK.data[2, 1] = _kout * sqrt(1 - x^2)
    varK.data[1, 2] = _kin

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

function measureAA(idx, var, obs, weight, config)
    kin, kout = var[3][1], var[4][1]
    nkin, nqout = var[5][1], var[6][1]
    obs[idx][1, kin, kout, nkin, nqout] += weight.d
    obs[idx][2, kin, kout, nkin, nqout] += weight.e
end

function AA(para::ParaMC, diagram;
    kamp=[para.kF,], kamp1=[kamp[1],],
    nkin=[0,],
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
    X = MCIntegration.Continuous(-1.0, 1.0) #x=cos(θ)

    Nkin, Nkout = length(kamp), length(kamp1)
    Nwin, Nwqout = length(nkin), length(nqout)

    vKin = MCIntegration.Discrete(1, Nkin, alpha=alpha)
    vKout = MCIntegration.Discrete(1, Nkout, alpha=alpha)

    vWin = MCIntegration.Discrete(1, Nwin, alpha=alpha)
    vWqout = MCIntegration.Discrete(1, Nwqout, alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nkin, Nkout, Nwin, Nwqout) for p in diagpara] # observable for the Fock diagram 

    if isnothing(config)
        config = MCIntegration.Configuration(; var=(K, T, vKin, vKout, vWin, vWqout, X),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(para, diag, root, extT, kamp, kamp1, nkin, nqout),
            kwargs...)
    end
    result = integrate(integrandAA; config=config, measure=measureAA, solver=:mcmc, neval=neval, print=print, kwargs...)

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