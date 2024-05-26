function integrand_lavg_Clib(idx, var, config)
    para, chan, filter, kampgrid, kamp2grid, qgrid, lgrid, n = config.userdata[1:8]
    maxMomNum, extT_labels, spin_conventions, leafStat, leafval, momLoopPool, root, partition = config.userdata[9:end]

    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    l = lgrid[var[4][1]]
    loopNum = config.dof[idx][1]
    extKidx = var[5][1]
    kamp = kampgrid[extKidx]
    kamp2 = kamp2grid[extKidx]
    qamp = qgrid[extKidx]

    if chan == :PH
        varK.data[1, 1] = kamp
        varK.data[1, 2] = kamp
        varK.data[2, 2] = qamp
        if dim == 3
            varK.data[1, 3] = kamp2 * x
            varK.data[2, 3] = kamp2 * sqrt(1 - x^2)
        else
            varK.data[1, 3] = kamp2 * cos(x)
            varK.data[2, 3] = kamp2 * sin(x)
        end
    elseif chan == :PP
        varK.data[1, 1] = kamp
        varK.data[1, 3] = -kamp
        varK.data[2, 3] = qamp
        if dim == 3
            varK.data[1, 2] = kamp2 * x
            varK.data[2, 2] = kamp2 * sqrt(1 - x^2)
        else
            varK.data[1, 2] = kamp2 * cos(x)
            varK.data[2, 2] = kamp2 * sin(x)
        end
    else
        error("not implemented")
    end

    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])

    for (i, lfstat) in enumerate(leafStat[idx])
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx, tau_num = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx, lfstat.tau_num
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            order = lforders[1]
            leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, para, lforders; idtype=Instant, tau_num=tau_num)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    group = partition[idx]
    if para.isDynamic
        evalfuncParquetADDynamic_map[group](root, leafval[idx])
    elseif Proper in filter
        evalfunc_vertex4Proper_map[group](root, leafval[idx])
    else
        evalfunc_vertex4_map[group](root, leafval[idx])
    end

    factor = legendfactor(x, l, dim)
    factor *= para.NF / (2π)^(dim * loopNum)

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

function lavg_Clib(para::ParaMC, diagram;
    chan=:PH,  # or :PP
    kamp=[para.kF,], #amplitude of k of the left legs
    kamp2=kamp, #amplitude of k of the right leg
    q=[0.0 for k in kamp],
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    integrand::Function=integrand_lavg_Clib,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    kwargs...
)
    partition, diagpara, extT_labels, spin_conventions = diagram
    filter = diagpara[1].filter
    MaxOrder = 5
    

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
        @assert filter== p.filter "filter should be the same"
    end

    dim, β, kF = para.dim, para.β, para.kF
    Nl = length(l)
    Nk = length(kamp)

    @assert length(diagpara) == length(partition) == length(extT_labels) == length(spin_conventions)
    
    maxMomNum = maximum([key[1] for key in partition]) + 3


    df = CSV.read(root_dir * "loopBasis_vertex4_maxOrder$(MaxOrder).csv", DataFrame)
    loopBasis = [df[!, col][1:maxMomNum] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)

    leafstates = Vector{Vector{LeafStateADDynamic}}()
    leafvalues = Vector{Vector{Float64}}()
    for key in partition
        key_str = join(string.(key))
        df = CSV.read(root_dir * "leafinfo_vertex4_$key_str.csv", DataFrame)
        leafstates_par = Vector{LeafStateADDynamic}()
        for row in eachrow(df)
            push!(leafstates_par, LeafStateADDynamic(row[2], _StringtoIntVector(row[3]), row[4:end]..., 1))
        end
        push!(leafstates, leafstates_par)
        push!(leafvalues, df[!, names(df)[1]])
    end

    root = zeros(Float64, maximum(length.(extT_labels)))
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
    obs = [zeros(ComplexF64, 2, Nl, Nk) for _ in diagpara]

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, L, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            # type=ComplexF64, # type of the integrand
            userdata=(para, chan,filter, kamp, kamp2, q, l, n, maxMomNum, extT_labels,
                spin_conventions, leafstates, leafvalues, momLoopPool,
                root, partition),
            kwargs...
        )
    end
    result = integrate(integrand; measure=measure_lavg, config=config, solver=solver, neval=neval, print=print, kwargs...)

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

function MC_lavg_Clib(para; kamp=[para.kF,], kamp2=kamp, q=[0.0 for k in kamp], n=[-1, 0, 0, -1], l=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],    # filter=[NoHartree, NoBubble, Proper],
    channel=:PH,
    partition=UEG.partition(para.order), neighbor = UEG.neighbor(partition),
    transferLoop=nothing,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    verbose=0
)
    kF = para.kF

    if Proper in filter
        root_dir=joinpath(@__DIR__, "source_codeParquetAD_Proper/")
    end

    diaginfo = Ver4.diagram_loadinfo(para, partition,
        filter=filter, transferLoop=transferLoop, root_dir = root_dir)

    partition = diaginfo[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    println(partition)
    

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
    end

    ver4, result = Ver4.lavg_Clib(para, diaginfo; chan = channel,
        kamp=kamp, kamp2=kamp2, q=q, n=n, l=l,
        neval=neval, print=verbose,
        neighbor=neighbor,
        root_dir = root_dir, reweight_goal=reweight_goal
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
