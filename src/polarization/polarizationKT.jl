function integrandKT(idx, vars, config)
    varK, varT, ExtTidx, ExtKidx = vars
    para, diag, kgrid, tgrid = config.userdata
    diagram = diag[idx]
    weight = diagram.node.current
    t = ExtTidx[1]
    k = ExtKidx[1]
    varK.data[1, 1] = kgrid[k]
    varT.data[2] = tgrid[t]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] for r in diagram.root)

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(para.dim * loopNum)
    return w * factor #the current implementation has an additional minus sign compared to the standard defintion
    # return real(w) * factor #the current implementation has an additional minus sign compared to the standard defintion
end


function measureKT(idx, vars, obs, weight, config)
    t = vars[3][1]  #imaginary time
    k = vars[4][1]  #K
    obs[idx][t, k] += weight
end

function KT(para::ParaMC, diagram;
    kgrid=[0.0,],
    tgrid=[para.β - 1e-8,], # must be (0, β)
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Polarization.KT"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, diag, root, _ = diagram

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=2, alpha=alpha)
    T.data[1] = 0.0
    T.data[2] = tgrid[1]
    ExtTidx = MCIntegration.Discrete(1, length(tgrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 2, 1, 1] for p in diagpara] # K, T, ExtTidx, ExtKidx
    # observable of polarization diagram of different permutations
    obs = [zeros(Float64, length(tgrid), length(kgrid)) for o in eachindex(dof)]

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = Configuration(;
            var=(K, T, ExtTidx, ExtKidx),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, diag, kgrid, tgrid),
            kwargs...
            # neighbor=neighbor,
            # reweight_goal=reweight_goal, kwargs...
        )
    end

    result = integrate(integrandKT; config=config, measure=measureKT, print=print, neval=neval, solver=solver, kwargs...)
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
            data = measurement.(avg, std)
            datadict[partition[1]] = data
        else
            for o in eachindex(dof)
                avg, std = result.mean[o], result.stdev[o]
                data = measurement.(avg, std)
                datadict[partition[o]] = data
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end
end