function integrand_parquetAD(idx, vars, config)
    varK, varT = vars
    para, maxMomNum, leafStat, leaf_maps, momLoopPool, root, funcGraphs!, partition, isLayered2D = config.userdata

    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    loopNum = config.dof[idx][1]
    is_zero_order = partition[idx] == (0, 0, 0) ? true : false

    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])
    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = leafOrders[idx][i][1]
            leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, order)
            if is_zero_order
                leafval[idx][i] *= (ϵ + μ)
            end
        elseif lftype == 2 #bosonic 
            diagid = leaf_maps[idx][i].properties
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            τ2, τ1 = varT[leafτ_o[idx][i]], varT[leafτ_i[idx][i]]
            idorder = leafOrders[idx][i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder;
                idtype=diagid.type, tau_num=interactionTauNum(diagid.type), isLayered=isLayered2D)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    funcGraphs![idx](root, leafval[idx])  # allocations due to run-time variable `idx`
    factor = 1.0 / (2π)^(dim * loopNum)
    return root[1] * factor
end

function measure(idx, var, obs, weight, config) # for the mcmc algorithm
    obs[idx] += weight
end

@inline function interactionTauNum(type::AnalyticProperty)
    if type == Instant
        return 1
    else
        return 2
    end
end

function ParquetAD(para::ParaMC, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    integrand::Function=integrand_parquetAD,
    kwargs...
)
    partition, diagpara, FeynGraphs = diagram

    @assert solver == :mcmc "Only :mcmc is supported for Sigma.GV"
    if para.isDynamic
        if NoBubble in diagpara[1].filter
            UEG.MCinitialize!(para, false)
        else
            UEG.MCinitialize!(para, true)
        end
    end

    for p in diagpara
        @assert diagpara[1].filter == p.filter "filter should be the same"
    end
    @assert length(diagpara) == length(FeynGraphs)

    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF

    maxMomNum = maximum([key[1] for key in partition]) + 1
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()

    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end
    leafStat, loopBasis = FeynmanDiagram.leafstates(leaf_maps, maxMomNum)
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)
    println("static compile has finished!")

    root = zeros(Float64, 1)
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 100.0 * kF)
    T = Continuous(0.0, β; alpha=alpha, adapt=true)

    dof = [[p.innerLoopNum, p.totalTauNum] for p in diagpara] # K, T
    obs = [zero(Float64) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, maxMomNum, leafStat, leaf_maps, momLoopPool, root, funcGraphs!, partition, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measure, print=print, neval=neval, solver=solver, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            if key[1] == 0
                datadict[key] = measurement(avg, std)
            else
                datadict[key] = measurement(avg, std) / β
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end
end
