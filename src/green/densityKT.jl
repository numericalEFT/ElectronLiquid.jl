function integrandDensityKT(idx, vars, config)
    # function integrandKW(idx, varK, varT, config)
    varK, varT = vars
    para, diag = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current

    # (ExtTin, ExtTout) = (0, β⁻)
    @assert varT.data[1] == 0
    @assert varT.data[2] == para.β - 1e-8

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] for r in diagram.root)

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(para.dim * loopNum)
    return w * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation of sigma has an additional minus sign compared to the standard defintion
end


function measureDensityKT(idx, vars, obs, weight, config)
    obs[idx] += weight
end

function densityKT(para::ParaMC, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    kwargs...
)
    # if haskey(kwargs, :solver)
    # @assert kwargs[:solver] == :mcmc "Only :mcmc is supported for Green.densityKT"
    # end
    @assert solver == :mcmc "Only :mcmc is supported for Green.densityKT"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, diag, root = diagram

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF)
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=2, alpha=alpha)
    ExtT = para.β - 1e-8  # normal-ordering, but must be (0, β)
    T.data[1] = 0.0
    T.data[2] = ExtT

    # NOTE: We integrate the external momentum loop in the density calculation
    dof = [[p.totalLoopNum, p.totalTauNum - 2] for p in diagpara] # K, T
    # observable of density diagram of different permutations
    obs = zeros(length(dof))

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(K, T),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, diag),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = integrate(integrandDensityKT; config=config, measure=measureDensityKT, print=print, neval=neval, solver=solver, kwargs...)
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

        for o in 1:length(dof)
            avg, std = result.mean[o], result.stdev[o]
            data = measurement.(avg, std)
            datadict[partition[o]] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end