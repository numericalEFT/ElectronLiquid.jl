"""
Calculate vertex4 averged on the Fermi surface
"""
function integrandPH_df(idx, vars, config)
    paras, diag, root, extT, kampgrid, lgrid, n = config.userdata

    kF, β = para.kF, para.β
    varK, varT, varX, varL, extKidx = var
    x = varX[1]
    # error("$(varK.data[:, 1])")
    l = lgrid[varL[1]]
    loopNum = config.dof[idx][1]
    kidx = extKidx[1]
    para = paras[kidx]
    kamp = kampgrid[kidx]
    # varK.data[1, 1], varK.data[1, 2] = kamp, kamp
    varK.data[1, 1], varK.data[1, 2] = kF, kF
    varK.data[:, 3] = [kamp * x, kamp * sqrt(1 - x^2), 0.0]


    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)

    factor = 1.0
    if l == 0
        factor = 1.0 / 2
    elseif l == 1
        factor = x / 2.0
    else
        error("not implemented")
    end

    if !isempty(rootuu)
        wuu = factor * sum(weight[root] * phase(varT, extTuu[ri], n, β) for (ri, root) in enumerate(rootuu))
    else
        wuu = zero(ComplexF64)
    end
    if !isempty(rootud)
        wud = factor * sum(weight[root] * phase(varT, extTud[ri], n, β) for (ri, root) in enumerate(rootud))
    else
        wud = zero(ComplexF64)
    end

    factor = para.NF / (2π)^(para.dim * loopNum)
    return Weight(wuu * factor, wud * factor)
end

function measurePH_df(idx, var, obs, weight, config)
    Lidx = var[4][1]
    Kidx = var[5][1]

    obs[idx][1, Lidx, Kidx] += weight.d
    obs[idx][2, Lidx, Kidx] += weight.e
end

function PH_df(paras::Vector{ParaMC}, diagram;
    kamp=[para.kF,], #amplitude of k of four external legs
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)

    para = paras[1]
    if paras[1].isDynamic
        for p in paras
            UEG.MCinitialize!(p)
        end
    end

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    Nl = length(l)
    Nk = length(kamp)

    partition, diagpara, diag, root, extT = diagram
    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    # K.data[:, 1] .= UEG.getK(kamp[1], para.dim, 1)
    # K.data[:, 2] .= UEG.getK(kamp[1], para.dim, 1)
    K.data[:, 1] .= UEG.getK(kF, para.dim, 1)
    K.data[:, 2] .= UEG.getK(kF, para.dim, 1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    L = MCIntegration.Discrete(1, Nl, alpha=alpha) # angular momentum
    AMP = MCIntegration.Discrete(1, Nk, alpha=alpha) # kamp

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nl, Nk) for p in diagpara]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, L, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(paras, diag, root, extT, kamp, l, n),
            kwargs...
        )
    end
    println(config.neighbor)
    result = integrate(integrandPH_df; measure=measurePH_df, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
            report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        end

        datadict = Dict{eltype(partition),Any}()
        if length(dof) == 1
            avg, std = result.mean, result.stdev
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[1]] = data
        else
            for k in 1:length(dof)
                avg, std = result.mean[k], result.stdev[k]
                r = measurement.(real(avg), real(std))
                i = measurement.(imag(avg), imag(std))
                data = Complex.(r, i)
                datadict[partition[k]] = data
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end

end
