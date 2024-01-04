function integrandGV(idx, vars, config)
    varK, varT = vars
    para, MaxLoopNum = config.userdata[1:2]
    leaf, leafType, leafτ_i, leafτ_o, leafMomIdx = config.userdata[3]
    LoopPool, root = config.userdata[4:5]
    graphfuncs! = config.userdata[6][idx]
    isLayered2D = config.userdata[end]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

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
            elseif order == -1
                leaf[idx][i] = Propagator.green(τ, ϵ, β) * (ϵ + μ)
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
                if isLayered2D == false
                    invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                    leaf[idx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                else
                    if order == 0
                        q = sqrt(dot(kq, kq) + 1e-16)
                        invK = 1.0 / q
                        leaf[idx][i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leaf[idx][i] = 0.0 # no high-order counterterms
                    end
                end
            else
                error("not implemented!")
            end
        end
    end

    graphfuncs!(root, leaf[idx])  # allocations due to run-time variable `idx`

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return root[1] * factor
end

function measureGV(idx, var, obs, weight, config) # for the mcmc algorithm
    obs[idx] += weight
end

function GV(para::ParaMC, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.GV"
    para.isDynamic && UEG.MCinitialize!(para)

    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF
    # partition, diagpara, FeynGraphs, FermiLabel, BoseLabel, mappings = diagram
    partition, diagpara, FeynGraphs, FermiLabel, BoseLabel, leafMap = diagram
    MaxLoopNum = maximum([key[1] for key in partition]) + 1
    LoopPool = FermiLabel.labels[3]

    leafStat, extT_labels = FeynmanDiagram.leafstates(FeynGraphs, FermiLabel, BoseLabel, partition)
    root = zeros(Float64, 1)
    funcGraphs! = Dict{Int,Function}(i => Compilers.compile(FeynGraphs[key][1], leafMap[key]) for (i, key) in enumerate(partition))

    # K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF)
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 100.0 * kF)
    T = Continuous(0.0, β; alpha=alpha, adapt=true)

    dof = [[p.innerLoopNum, p.totalTauNum] for p in diagpara] # K, T, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zero(Float64) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, MaxLoopNum, leafStat, LoopPool, root, funcGraphs!, isLayered2D),
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
            avg, std = result.mean[o], result.stdev[o]
            datadict[key] = measurement(avg, std) / β
        end
        return datadict, result
    else
        return nothing, nothing
    end
end