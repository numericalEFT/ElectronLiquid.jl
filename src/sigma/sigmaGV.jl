function integrandGV(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx = config.userdata[6]
    interactionleaf = config.userdata[7][idx]
    LoopPool, root = config.userdata[8:9]
    graphfuncs! = config.userdata[10][idx]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    # varK.data[:, 2:MaxLoopNum] .= [para.kF, 0.0]
    FrontEnds.update(LoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif isodd(lftype) #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = (lftype - 1) / 2
            if order == 0
                leaf[idx][i] = Propagator.green(τ, ϵ, β)
            elseif order == 1
                leaf[idx][i] = -Spectral.kernelFermiT_dω(τ, ϵ, β)
            elseif order == 2
                leaf[idx][i] = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
            elseif order == 3
                leaf[idx][i] = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
            elseif order == 4
                leaf[idx][i] = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
            elseif order == 5
                leaf[idx][i] = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
            else
                error("not implemented!")
            end
        else
            kq = FrontEnds.loop(LoopPool, leafMomIdx[idx][i])
            order = lftype / 2 - 1
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leaf[idx][i] = e0^2 / ϵ0 * invK * (λ * invK)^order
            elseif dim == 2
                invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                leaf[idx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
            else
                error("not implemented!")
            end
        end
    end

    graphfuncs!(root, leaf[idx], interactionleaf)  # allocations due to run-time variable `idx`

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function measureGV(idx, vars, obs, weight, config) # for the mcmc algorithm
    n = vars[3][1]  #matsubara frequency
    k = vars[4][1]  #K
    obs[idx][n, k] += weight
end

# function measureGV(vars, obs, weight, config) # for vegas and vegasmc algorithms
#     N = length(config.dof)
#     n = vars[3][1]  #matsubara frequency
#     k = vars[4][1]  #K
#     for idx in 1:N
#         obs[idx][n, k] += weight
#     end
# end

function GV(para::ParaMC, diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.GV"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, FeynGraphs, FermiLabel, BoseLabel, mappings = diagram
    MaxLoopNum = maximum([key[1] for key in partition]) + 2
    LoopPool = FermiLabel.labels[3]
    PropagatorMap, InteractionMap = mappings

    PropagatorStat, InteractionStat, extT_labels = FeynmanDiagram.LeavesState(FeynGraphs, FermiLabel, BoseLabel, partition)
    root = zeros(Float64, 24)
    funcGraphs! = Dict{Int,Function}(i => Compilers.compile(FeynGraphs[key][1], PropagatorMap[key], InteractionMap[key]) for (i, key) in enumerate(partition))

    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    # T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=2)
    T.data[1:2] .= 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), adapt=false)
    # ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum, 1, 1] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, ngrid, MaxLoopNum, extT_labels, PropagatorStat, InteractionStat[1], LoopPool, root, funcGraphs!),
            kwargs...
        )
    end

    result = integrate(integrandGV; config=config, measure=measureGV, print=print, neval=neval, solver=solver, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end
        if print >= -2
            println(result)
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            # if length(dof) == 1
            #     avg, std = result.mean, result.stdev
            # else
            avg, std = result.mean[o], result.stdev[o]
            # end
            r = measurement.(real(avg), real(std)) ./ (-2β)
            i = measurement.(imag(avg), imag(std)) ./ (-2β)
            data = Complex.(r, i)
            datadict[key] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end