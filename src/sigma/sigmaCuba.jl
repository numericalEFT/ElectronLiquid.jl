"""
This example demonstrate how to use Cuba package to calculate the diagrams
"""

function integrandCuba(x, f, userdata)
    para, diagram, kgrid, ngrid, varK, varT, loopNum, tauNum = userdata
    weight = diagram.node.current
    object = diagram.node.object
    # vK = reshape(view(x, 1:loopNum*para.dim), para.dim, loopNum)
    varT[1] = 0.0
    for i = 1:tauNum
        varT[1+i] = x[i+loopNum*para.dim] * para.β
    end

    varK[1, 1] = kgrid[1]
    factor = (1.0 / (2π)^(para.dim) * 2π^2)^loopNum * (para.β)^tauNum
    for i = 1:loopNum
        shift = 3 * (i - 1)
        r = x[shift+1] / (1 - x[shift+1])
        θ = x[shift+2] * π
        ϕ = x[shift+3] * 2 * π
        varK[:, i+1] .= [r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)]
        factor *= r^2 / (1 - x[shift+1])^2 * sin(θ)
    end
    ExprTree.evalKT!(diagram, varK, varT, para)
    for l in eachindex(ngrid)
        wn = ngrid[l]
        w = zero(ComplexF64)
        for r in diagram.root
            w += weight[r] * phase(varT, object[r].para.extT::Tuple{Int,Int}, wn, para.β)
        end
        w *= factor
        # w = factor * sum(weight[r] * phase(varT, object[r].para.extT, wn, para.β) for r in diagram.root)
        f[2l-1] = real(w)
        f[2l] = imag(w)
    end
end

function cuba(para::ParaMC, diagram;
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

    # K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    # K.data[:, 1] .= 0.0
    # K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    # T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    # T.data[1] = 0.0
    # X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    # ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    for (pidx, p) in enumerate(diagpara)
        loopNum, tauNum = p.innerLoopNum, p.totalTauNum - 1
        varK = zeros(para.dim, loopNum + 1)
        varT = zeros(tauNum + 1)
        ndim = para.dim * loopNum + tauNum
        ncomp = 4
        result = vegas(integrandCuba, ndim, ncomp; maxevals=neval, userdata=(para, diag[pidx], kgrid, ngrid, varK, varT, loopNum, tauNum), kwargs...)
        # result = suave(integrandCuba, ndim, ncomp; maxevals=neval, userdata=(para, diag[pidx], kgrid, ngrid, varK, varT, loopNum, tauNum))
        # result = cuhre(integrandCuba, ndim, ncomp;
        #     maxevals=neval,
        #     userdata=(para, diag[pidx], kgrid, ngrid, varK, varT, loopNum, tauNum)
        # )
        println(result)
        # return result
    end

    # dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    # obs = zeros(ComplexF64, length(dof), length(ngrid), length(kgrid)) # observable for the Fock diagram 

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    # if isnothing(config)
    #     config = MCIntegration.Configuration(; var=(K, T, X, ExtKidx),
    #         dof=dof,
    #         obs=obs,
    #         para=(para, diag, kgrid, ngrid),
    #         kwargs...
    #         # neighbor=neighbor,
    #         # reweight_goal=reweight_goal, kwargs...
    #     )
    # end

    # result = MCIntegration.sample(config, integrandKW, measureKW; print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    # if isnothing(result) == false

    #     avg, std = result.mean, result.stdev
    #     if print >= 0
    #         MCIntegration.summary(result.config)
    #         println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))
    #     end

    #     r = measurement.(real(avg), real(std))
    #     i = measurement.(imag(avg), imag(std))
    #     data = Complex.(r, i)
    #     datadict = Dict{eltype(partition),typeof(data[1, :, :])}()
    #     for i in 1:length(dof)
    #         datadict[partition[i]] = data[i, :, :]
    #     end
    #     return datadict, result
    # else
    #     return nothing, nothing
    # end
end