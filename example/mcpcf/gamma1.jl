using ElectronLiquid
using ElectronLiquid.CompositeGrids
using ElectronLiquid.UEG
using ElectronLiquid.Ver4: phase
using ElectronLiquid.Measurements
using FeynmanDiagram
using MCIntegration

function integrandFG(idx, var, config)
    para, diag, root, extT, L, n = config.userdata

    kF, β = para.kF, para.β
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    # error("$(varK.data[:, 1])")
    l = L[1]
    loopNum = config.dof[idx][1]
    # kamp = kampgrid[var[5][1]]
    # varK.data[1, 1], varK.data[1, 2] = kamp, kamp
    varK.data[1, 1], varK.data[1, 3] = kF, -kF
    # varK.data[:, 3] = [kamp * x, kamp * sqrt(1 - x^2), 0.0]
    varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]


    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)

    factor = 1.0
    if l == 0
        # factor = 1.0 / 2
        factor = 1.0
    elseif l == 1
        # factor = x / 2.0
        factor = x
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

function measureFG(idx, var, obs, weight, config)
    Lidx = var[4][1]
    Kidx = var[5][1]

    obs[idx][1, Lidx, Kidx] += weight.d
    obs[idx][2, Lidx, Kidx] += weight.e
end

function Full_Gamma(para::ParaMC, diagram;
    n=[0, 0, 0],
    l=0,
    neval=1e6,
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...)
    UEG.MCinitialize!(para)

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    Nl = 1
    Nk = 1

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
    obs = [zeros(ComplexF64, 2, Nl, Nk) for p in diagpara]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration(
            var=(K, T, X, L, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(para, diag, root, extT, l, n),
            kwargs...
        )
    end
    result = integrate(integrandFG; measure=measureFG, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

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

rs = 3.0
beta = 25
Fs = -0.433
order = 2
neval = 1e7

paramc = UEG.ParaMC(rs=rs, beta=beta, Fs=Fs, order=order, isDynamic=true)
UEG.MCinitialize!(paramc)

println(paramc)

partition = [(2, 0, 0)]
channel = [PHr, PHEr, PPr]
neighbor = UEG.neighbor(partition)
filter = [NoHartree, NoBubble]

diagram = Ver4.diagram(paramc, partition; channel=channel, filter=filter)

data, result = Full_Gamma(paramc, diagram; n=[0, 0, 0], neval=1e7)