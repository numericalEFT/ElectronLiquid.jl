"""
Calculate vertex4 averged on the Fermi surface
"""
function _integrandGeneric(idx, var, config)
    weight, factor = _diagram_weight(idx, var, config)
    return weight
end

function _diagram_weight(idx, var, config)
    paras, diag, root, extT = config.userdata

    varK, varT, varX, varN = var[1], var[2], var[3], var[4]

    x = varX[3][1]
    n = varN[4][1]

    para = paras[n]
    l = para.l
    dim = para.para.dim

    k1, k2 = para.kamp
    if para.channel == :PH
        varK.data[1, 1], varK.data[1, 2] = k1, k1
        varK.data[1, 3] = k2 * x
        varK.data[2, 3] = k2 * sqrt(1 - x^2)
    elseif para.channel == :PP
        varK.data[1, 1], varK.data[1, 3] = k1, -k1
        varK.data[1, 2] = k2 * x
        varK.data[2, 2] = k2 * sqrt(1 - x^2)
    else
        error("not implemented")
    end

    loopNum = config.dof[idx][1]
    diagram = diag[idx]
    weight = diagram.node.current

    ExprTree.evalKT!(diagram, varK.data, varT.data, para.para)

    factor = para.para.NF / (2π)^(dim * loopNum)
    if l == 0
        factor *= 1.0 / 2
    elseif l == 1
        factor *= x / 2.0
    else
        error("not implemented")
    end

    rootuu, rootud = root[1][idx], root[2][idx]
    if !isempty(rootuu)
        wuu = factor * sum(weight[root] for (ri, root) in enumerate(rootuu))
    else
        wuu = zero(ComplexF64)
    end
    if !isempty(rootud)
        wud = factor * sum(weight[root] for (ri, root) in enumerate(rootud))
    else
        wud = zero(ComplexF64)
    end

    return Weight(wuu, wud), factor
end

function _measureGeneric(idx, var, obs, relative_weight, config)
    #relative_weight = abs(w)/(probability of sampling this diagram)
    w, factor = _diagram_weight(idx, var, config)
    paras, diag, root, extT = config.userdata
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]
    diagram = diag[idx]
    weight = diagram.node.current

    varT = var[2]
    n = var[4][1]
    para = paras[n]
    ωn = para.ωn #get the frequency grid to measure
    β = para.para.β

    for i in eachindex(ωn)
        if !isempty(rootuu)
            obs[idx][1, i, n] += factor * sum(weight[root] * phase(varT, extTuu[ri], ωn[i], β) for (ri, root) in enumerate(rootuu)) * relative_weight.d / abs(w)
        end
        if !isempty(rootud)
            obs[idx][2, i, n] += factor * sum(weight[root] * phase(varT, extTud[ri], ωn[i], β) for (ri, root) in enumerate(rootud)) * relative_weight.e / abs(w)
        end
    end
end

function one_angle_averaged(paras::Vector{OneAngleAveraged}, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)

    dim, β, kF = paras[1].para.dim, paras[1].para.β, paras[1].para.kF
    Nw = length(paras[1].ωn)
    for p in paras
        p.para.isDynamic && UEG.MCinitialize!(p.para)
        @assert length(p.ωn) == Nw "All parameters must have the same frequency list"
        @assert p.para.dim == dim "All parameters must have the same dimension"
        @assert p.para.β ≈ β "All parameters must have the same inverse temperature"
        @assert p.para.kF ≈ kF "All parameters must have the same Fermi momentum"
    end

    partition, diagpara, diag, root, extT = diagram
    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    # all momentum will be sampled around the dimensionless Fermi momentum 1.0
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3) # the first three momenta are external
    K.data[:, 1] = UEG.getK(kF, dim, 1)
    K.data[:, 2] = UEG.getK(kF, dim, 1)
    K.data[:, 3] = UEG.getK(kF, dim, 1)
    # all time variables will be sampled within [0.0, 1.0]
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha) # the first one is external
    T.data[1] = 0.0
    X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    N = MCIntegration.Discrete(1, length(paras), alpha=alpha) #index of paras

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nw, length(paras)) for p in diagpara]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration(
            var=(K, T, X, N),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(paras, diag, root, extT),
            kwargs...
        )
    end
    result = integrate(_integrandGeneric; measure=_measureGeneric, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

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
