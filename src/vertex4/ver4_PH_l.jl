"""
Calculate vertex4 averged on the Fermi surface
"""
function integrandPH(config)
    para, diag, root, extT, kampgrid, lgrid, n = config.para

    kF, β = para.kF, para.β
    idx = config.curr
    var = config.var
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    # error("$(varK.data[:, 1])")
    l = lgrid[var[4][1]]
    loopNum = config.dof[idx][1]
    kamp = kampgrid[var[5][1]]
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

function measurePH(config)
    factor = 1.0 / config.reweight[config.curr]

    Lidx = config.var[4][1]
    Kidx = config.var[5][1]

    o = config.curr
    weight = integrandPH(config)
    config.observable[o, 1, Lidx, Kidx] += weight.d / abs(weight) * factor
    config.observable[o, 2, Lidx, Kidx] += weight.e / abs(weight) * factor
end

function PH(para::ParaMC, diagram;
    kamp=[para.kF,], #amplitude of k of four external legs
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

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
    AMP = MCIntegration.Discrete(1, Nk, alpha=alpha) # angular momentum

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), 2, Nl, Nk)

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration((K, T, X, L, AMP), dof, obs;
            para=(para, diag, root, extT, kamp, l, n),
            kwargs...
        )
    end
    result = MCIntegration.sample(config, integrandPH, measurePH; neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        if print >= 0
            MCIntegration.summary(result.config)
            MCIntegration.summary(result, [o -> (real(o[i, 1, 1, 1])) for i in 1:length(dof)], ["uu$i" for i in 1:length(dof)])
            MCIntegration.summary(result, [o -> (real(o[i, 2, 1, 1])) for i in 1:length(dof)], ["ud$i" for i in 1:length(dof)])
        end

        avg, std = result.mean, result.stdev
        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        data = Complex.(r, i)
        datadict = Dict{eltype(partition),typeof(data[1, :, :, :])}()
        for i in 1:length(dof)
            datadict[partition[i]] = data[i, :, :, :]
        end
        return datadict, result
    else
        return nothing, nothing
    end

end
