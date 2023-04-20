
function integrandKW_df(idx, vars, config)
    # function integrandKW(idx, varK, varT, config)
    varK, varT, N, ExtKidx = vars
    paras, diag, extT, kgrid, ngrid = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current
    l = N[1]
    k = ExtKidx[1]
    para = paras[k]
    varK.data[1, 1] = kgrid[k]
    wn = ngrid[l]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, extT[idx][ri], wn, para.β) for (ri, r) in enumerate(diagram.root))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(para.dim * loopNum)
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

#for vegas algorithm
# function integrandKW(varK, varT, N, ExtKidx, config)
#     return integrandKW(1, varK, varT, N, ExtKidx, config)
# end

function measureKW_df(idx, vars, obs, weight, config)
    l = vars[3][1]  #matsubara frequency
    k = vars[4][1]  #K
    obs[idx][l, k] += weight
end

#for vegas algorithm
# function measureKW(obs, weights, config)
#     l = config.var[3][1]  #matsubara frequency
#     k = config.var[4][1]  #K
#     for o in 1:config.N
#         obs[o][l, k] += weights[o]
#     end
# end

function KW_df(paras::Vector{ParaMC}, # we assume that all paras share the same kF and beta 
    diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    kwargs...
)
    # if haskey(kwargs, :solver)
    # @assert kwargs[:solver] == :mcmc "Only :mcmc is supported for Sigma.KW"
    # end
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.KW"
    if paras[1].isDynamic
        for p in paras
            UEG.MCinitialize!(p)
        end
    end

    dim, β, kF = paras[1].dim, paras[1].β, paras[1].kF
    partition, diagpara, diag, root, extT = diagram

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)) for o in 1:length(dof)]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(paras, diag, extT, kgrid, ngrid),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = integrate(integrandKW_df; config=config, measure=measureKW_df, print=print, neval=neval, solver=solver, kwargs...)
    # result = integrate(integrandKW; config=config, print=print, neval=neval, kwargs...)
    # niter=niter, print=print, block=block, kwargs...)

    if isnothing(result) == false

        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end

        if print >= -2
            println(result)
        end

        # datadict = Dict{eltype(partition),Complex{Measurement{Float64}}}()
        datadict = Dict{eltype(partition),Any}()

        if length(dof) == 1
            avg, std = result.mean, result.stdev
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[1]] = data
        else
            for o in 1:length(dof)
                avg, std = result.mean[o], result.stdev[o]
                r = measurement.(real(avg), real(std))
                i = measurement.(imag(avg), imag(std))
                data = Complex.(r, i)
                datadict[partition[o]] = data
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end
end