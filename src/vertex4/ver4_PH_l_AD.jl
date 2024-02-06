"""
Calculate vertex4 averged on the Fermi surface
"""

function integrandPH_AD(idx, var, config)
    para, kampgrid, kamp2grid, qgrid, lgrid, n = config.userdata[1:6]
    MaxLoopNum, extT_labels, spin_conventions, leafStat, leaf_maps, momLoopPool, root, funcGraphs! = config.userdata[7:end]
    # MaxLoopNum, extT_labels, leafStat, leaf_maps, momLoopPool, root, funcGraphs! = config.userdata[7:end]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    # error("$(varK.data[:, 1])")
    l = lgrid[var[4][1]]
    loopNum = config.dof[idx][1]
    extKidx = var[5][1]
    kamp = kampgrid[extKidx]
    kamp2 = kamp2grid[extKidx]
    qamp = qgrid[extKidx]
    varK.data[1, 1] = kamp
    varK.data[1, 2] = kamp
    varK.data[2, 2] = qamp
    #varK.data[1, 1], varK.data[1, 2] = kF, kF
    if dim == 3
        varK.data[1, 3] = kamp2 * x
        varK.data[2, 3] = kamp2 * sqrt(1 - x^2)
    else
        varK.data[1, 3] = kamp2 * cos(x)
        varK.data[2, 3] = kamp2 * sin(x)
    end

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])

    # println(varK.data[:, 1:MaxLoopNum])

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
            # idorder = diagid.order
            idorder = leafOrders[idx][i]
            # @assert idorder == leafOrders[idx][i]
            # leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder; idtype=Instant, tau_num=interactionTauNum(diagid.para))
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, idorder; idtype=Instant, tau_num=1)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    factor = legendfactor(x, l, dim)
    factor *= para.NF / (2π)^(dim * loopNum)

    graphfuncs! = funcGraphs![idx]
    graphfuncs!(root, leafval[idx])

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

    # for ri in 1:length(extT_labels[idx])
    #     wuu += root[ri] * phase(varT, extT_labels[idx][ri], n, β)
    # end
    # return wuu * factor
end

function measurePH_AD(idx, var, obs, weight, config)
    Lidx = var[4][1]
    Kidx = var[5][1]
    # obs[idx][Lidx, Kidx] += weight
    obs[idx][1, Lidx, Kidx] += weight.d
    obs[idx][2, Lidx, Kidx] += weight.e
end

function PH_AD(para::ParaMC, diagram;
    kamp=[para.kF,], #amplitude of k of the left legs
    kamp2=kamp, #amplitude of k of the right leg
    q=[0.0 for k in kamp],
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    root_dir=joinpath(@__DIR__, "source_codeAD/"),
    config=nothing,
    kwargs...
)
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram
    # partition, diagpara, FeynGraphs, extT_labels = diagram

    # UEG.MCinitialize!(para)
    if NoBubble in diagpara[1].filter
        UEG.MCinitialize!(para, false)
    else
        UEG.MCinitialize!(para, true)
    end

    for p in diagpara
        @assert diagpara[1].filter == p.filter "filter should be the same"
    end

    dim, β, kF, order = para.dim, para.β, para.kF, para.order
    Nl = length(l)
    Nk = length(kamp)

    # println(spin_conventions)

    @assert length(diagpara) == length(FeynGraphs) == length(extT_labels) == length(spin_conventions)
    # @assert length(diagpara) == length(FeynGraphs) == length(extT_labels)

    MaxLoopNum = maximum([key[1] for key in partition]) + 3
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()

    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key][1])
        push!(leaf_maps, leafmap)
    end
    leafStat, loopBasis = FeynmanDiagram.leafstates(leaf_maps, MaxLoopNum)
    # println(leafStat[3][1],partition[1])


    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)
    root = zeros(Float64, maximum(length.(extT_labels)))
    println("static compile has finished!")

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    K.data[:, 1] .= UEG.getK(kF, dim, 1)
    K.data[:, 2] .= UEG.getK(kF, dim, 1)
    K.data[:, 3] .= 0.0
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0

    if dim == 3
        X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    elseif dim == 2
        X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
    end
    L = MCIntegration.Discrete(1, Nl, alpha=alpha) # angular momentum
    AMP = MCIntegration.Discrete(1, Nk, alpha=alpha) # angular momentum

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nl, Nk) for p in diagpara]

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, L, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            # type=ComplexF64, # type of the integrand
            userdata=(para, kamp, kamp2, q, l, n, MaxLoopNum, extT_labels,
                spin_conventions, leafStat, leaf_maps, momLoopPool,
                # leafStat, leaf_maps, momLoopPool,
                root, funcGraphs!),
            kwargs...
        )
    end
    result = integrate(integrandPH_AD; measure=measurePH_AD, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        # if print >= 0
        #     report(result.config)
        #     report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
        #     report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        # end

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

function MC_PH_AD(para; kamp=[para.kF,], kamp2=kamp, q=[0.0 for k in kamp], n=[-1, 0, 0, -1], l=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],    # filter=[NoHartree, NoBubble, Proper],
    channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order), generate_type=:Parquet,
    verbose=0
)

    kF = para.kF
    _order = para.order

    # partition = UEG.partition(_order)


    if generate_type == :GV
        diagram = Ver4.diagramGV(para, partition; channels=channels, filter=filter)
    elseif generate_type == :Parquet
        diagram = Ver4.diagramParquet(para, partition; channels=channels, filter=filter)
    else
        error("generate_type $generate_type not implemented!")
    end

    partition = diagram[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    neighbor = UEG.neighbor(partition)

    println(partition)

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
        println(length(reweight_goal))
    end

    ver4, result = Ver4.PH_AD(para, diagram;
        kamp=kamp, kamp2=kamp2, q=q, n=n, l=l,
        neval=neval, print=verbose,
        neighbor=neighbor,
        reweight_goal=reweight_goal
    )

    if isnothing(ver4) == false
        # for (p, data) in ver4
        for p in partition
            data = ver4[p]
            printstyled("partition: $p\n", color=:yellow)
            for (li, _l) in enumerate(l)
                printstyled("l = $_l\n", color=:green)
                @printf("%12s    %16s    %16s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "di", "ex", "symmetric", "asymmetric")
                for (ki, k) in enumerate(kamp)
                    factor = 1.0
                    d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
                    s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
                    di, ex = (s - a), (a) * 2.0
                    @printf("%12.6f    %16s    %16s    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$di", "$ex", "$s", "$a")
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
                f[key] = (kamp, n, l, ver4)
            end
        end
    end
    return ver4, result
end
