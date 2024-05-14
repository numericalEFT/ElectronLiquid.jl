function integrand_Spec_Jl(idx, var, config)
    para, chan,  kampgrid, kamp2grid, Anggrid, phigrid, n = config.userdata[1:7]
    maxMomNum, extT_labels, spin_conventions, leafStat, leaf_maps, momLoopPool, root, funcGraphs! = config.userdata[8:end]

    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
    varK, varT = var[1], var[2]
    # x = config.var[3][1]
    x = Anggrid[var[3][1]]
    y = phigrid[var[4][1]]
    loopNum = config.dof[idx][1]
    extKidx = var[5][1]
    kamp = kampgrid[extKidx]
    kamp2 = kamp2grid[extKidx]
    # qamp = qgrid[extKidx]

    if chan == :PH
        varK.data[1, 1] = kamp * cos(x / 2)
        varK.data[2, 1] = kamp * sin(x / 2)
        varK.data[1, 2] = kamp * cos(x / 2)
        varK.data[2, 2] = kamp * sin(x / 2) * cos(y)
        varK.data[3, 2] = kamp * sin(x / 2) * sin(y)
        varK.data[1, 3] = kamp2 * cos(x / 2)
        varK.data[2, 3] = -kamp2 * sin(x / 2)
    elseif chan == :PP
        varK.data[1, 1] = kamp * cos(x / 2)
        varK.data[2, 1] = kamp * sin(x / 2)
        varK.data[1, 3] = -kamp * cos(x / 2)
        varK.data[2, 3] = -kamp * sin(x / 2) * cos(y)
        varK.data[3, 3] = -kamp * sin(x / 2) * sin(y)
        varK.data[1, 2] = kamp2 * cos(x / 2)
        varK.data[2, 2] = -kamp2 * sin(x / 2)
    else
        error("not implemented")
    end

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
        elseif lftype == 2 #bosonic 
            diagid = leaf_maps[idx][i].properties
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            τ2, τ1 = varT[leafτ_o[idx][i]], varT[leafτ_i[idx][i]]
            idorder = leafOrders[idx][i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder; idtype=diagid.type, tau_num=interactionTauNum(diagid.type))
            # leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder; idtype=Instant, tau_num=1)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs! = funcGraphs![idx]
    graphfuncs!(root, leafval[idx])

    # factor = legendfactor(x, l, dim)
    factor = para.NF / (2π)^(dim * loopNum)

    wuu = zero(ComplexF64)
    wud = zero(ComplexF64)
    for ri in 1:length(extT_labels[idx])
        if spin_conventions[idx][ri] == UpUp
            wuu += root[ri] * phase(varT, extT_labels[idx][ri], n, β)
        elseif spin_conventions[idx][ri] == UpDown
            wud += root[ri] * phase(varT, extT_labels[idx][ri], n, β)
        end
    end

    return Weight(wuu * factor, wud * factor)
end

function Spec_Jl(para::ParaMC, diagram;
    chan=:PH,  # or :PP
    kamp=[para.kF,], #amplitude of k of the left legs
    kamp2=kamp, #amplitude of k of the right leg
    n=[0, 0, 0],
    theta=[0,],
    phi=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    integrand::Function=integrand_Spec_Jl,
    kwargs...
)
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram

    # if para.isDynamic
    #     root_dir = joinpath(@__DIR__, "source_codeParquetAD/dynamic/")
    # end

    # UEG.MCinitialize!(para)
    if NoBubble in diagpara[1].filter
        UEG.MCinitialize!(para, false)
    else
        UEG.MCinitialize!(para, true)
    end

    for p in diagpara
        @assert diagpara[1].filter == p.filter "filter should be the same"
    end

    dim, β, kF = para.dim, para.β, para.kF
    NA = length(theta)
    NP = length(phi)
    Nk = length(kamp)

    @assert length(diagpara) == length(partition) == length(extT_labels) == length(spin_conventions)
    @assert dim == 3

    maxMomNum = maximum([key[1] for key in partition]) + 3

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()

    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end
    leafStat, loopBasis = FeynmanDiagram.leafstates(leaf_maps, maxMomNum)

    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)
    root = zeros(Float64, maximum(length.(extT_labels)))
    println("static compile has finished!")

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    K.data[:, 1] .= UEG.getK(kF, dim, 1)
    K.data[:, 2] .= UEG.getK(kF, dim, 1)
    K.data[:, 3] .= 0.0
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0

    # if dim == 3
    #     X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    # elseif dim == 2
    #     X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
    # end
    Ang = MCIntegration.Discrete(1, NA, alpha=alpha)
    phidx = MCIntegration.Discrete(1, NP, alpha=alpha)
    AMP = MCIntegration.Discrete(1, Nk, alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, NA, NP, Nk) for _ in diagpara]

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, Ang, phidx, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            # type=ComplexF64, # type of the integrand
            userdata=(para, chan, kamp, kamp2, theta, phi, n, maxMomNum, extT_labels,
                spin_conventions, leafStat, leaf_maps, momLoopPool,
                root, funcGraphs!),
            kwargs...
        )
    end
    result = integrate(integrand; measure=measure_Spec, config=config, solver=solver, neval=neval, print=print, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            # report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
            # report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        end

        datadict = Dict{eltype(partition),Any}()
        for k in 1:length(dof)
            avg, std = result.mean[k], result.stdev[k]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[partition[k]] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end

end

function MC_Spec_Jl(para; kamp=[para.kF,], kamp2=kamp, n=[0, 0, 0, 0], theta=[0,], phi=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],    # filter=[NoHartree, NoBubble, Proper],
    channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order),
    transferLoop=nothing, extK=nothing, optimize_level=1,
    verbose=0
)
    kF = para.kF


    diagram = Diagram.diagram_parquet_response(:vertex4, para, partition,
        channels=channels, filter=filter, extK=extK, transferLoop=transferLoop, optimize_level=optimize_level)

    partition = diagram[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    println(partition)
    neighbor = UEG.neighbor(partition)


    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
    end

    ver4, result = Ver4.Spec_Jl(para, diagram;
        kamp=kamp, kamp2=kamp2, n=n, theta=theta, phi=phi,
        neval=neval, print=verbose,
        neighbor=neighbor,
        reweight_goal=reweight_goal
    )

    if isnothing(ver4) == false && verbose != -1
        # for (p, data) in ver4
        for p in partition
            data = ver4[p]
            printstyled("partition: $p\n", color=:yellow)
            for (ti, _theta) in enumerate(theta)
                printstyled("θ = $(_theta/π) π\n", color=:green)
                for (pidx, _phi) in enumerate(phi)
                    printstyled(" ϕ= $(_phi/π) π\n", color=:blue)
                    @printf("%12s    %16s    %16s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "di", "ex", "symmetric", "asymmetric")
                    for (ki, k) in enumerate(kamp)
                        factor = 1.0
                        d1, d2 = real(data[1, ti, pidx, ki]) * factor, real(data[2, ti, pidx, ki]) * factor
                        s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
                        di, ex = (s - a), (a) * 2.0
                        @printf("%12.6f    %16s    %16s    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$di", "$ex", "$s", "$a")
                    end
                end
            end
        end

        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kamp, n, theta, phi, ver4)
            end
        end
    end
    return ver4, result
end
