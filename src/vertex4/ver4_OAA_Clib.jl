function integrand_OAA_Clib(idx, var, config)
    weight, factor = diagram_weight_OAA_Clib(idx, var, config)
    return weight
end

function diagram_weight_OAA_Clib(idx, var, config)
    paras, maxMomNum, extT_labels, spin_conventions, leafstates, leafval, momLoopPool, root, partition = config.userdata
    varK, varT, varX, varN = var[1], var[2], var[3], var[4]

    x = varX[1]
    n = varN[1]
    loopNum = config.dof[idx][1]

    para = paras[n]
    l = para.l
    param = para.para
    dim, β, me, λ, μ, e0, ϵ0 = param.dim, param.β, param.me, param.mass2, param.μ, param.e0, param.ϵ0

    k1, k2 = para.kamp
    if para.channel == :PH
        varK.data[1, 1], varK.data[1, 2] = k1, k1
        varK.data[1, 3] = k2 * x
        varK.data[2, 3] = k2 * sqrt(1 - x^2)
    elseif para.channel == :PP
        varK.data[1, 1], varK.data[1, 3] = k1, -k1
        varK.data[1, 2] = k2 * x
        varK.data[2, 2] = k2 * sqrt(1 - x^2)
    else
        error("not implemented")
    end

    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])
    # println(momLoopPool)
    for (i, lfstat) in enumerate(leafstates[idx])
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
            # println(kq, (k1, k2, x))
            # @assert dot(kq, kq) ≈ (k1^2 + k2^2 - 2k1 * k2 * x) "$(dot(kq, kq)) != $(k1^2 + k2^2 - 2k1 * k2 * x)"
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Instant, tau_num=tau_num)
        elseif lftype == 4 # dynamic bosonic
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            # println(kq, (k1, k2, x))
            # @assert dot(kq, kq) ≈ (k1^2 + k2^2 - 2k1 * k2 * x) "$(dot(kq, kq)) != $(k1^2 + k2^2 - 2k1 * k2 * x)"
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Dynamic, tau_num=tau_num)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    group = partition[idx]
    if param.isDynamic
        evalfuncParquetADDynamic_map[group](root, leafval[idx])
    else
        evalfunc_vertex4_map[group](root, leafval[idx])
    end

    factor = para.para.NF / (2π)^(dim * loopNum)
    factor *= legendfactor(x, l, dim)

    wuu = zero(ComplexF64)
    wud = zero(ComplexF64)
    for ri in 1:length(extT_labels[idx])
        if spin_conventions[idx][ri] == FeynmanDiagram.UpUp
            wuu += root[ri]
        elseif spin_conventions[idx][ri] == FeynmanDiagram.UpDown
            wud += root[ri]
        end
    end
    # println(wuu, wud)
    wuu, wud = factor * wuu, factor * wud
    return Weight{ComplexF64}(wuu, wud), factor
end

function measure_OAA_Clib(idx, var, obs, relative_weight, config)
    w, factor = diagram_weight_OAA_Clib(idx, var, config)
    inverse_probability = abs(relative_weight) / abs(w)
    paras, maxMomNum, extT_labels, spin_conventions, leafstates, leafvalues, momLoopPool, root, _ = config.userdata

    varT = var[2]
    n = var[4][1]
    para = paras[n]
    ωn = para.ωn #get the frequency grid to measure
    β = para.para.β
    for i in eachindex(ωn)
        wuu = zero(ComplexF64)
        wud = zero(ComplexF64)
        for ri in 1:length(extT_labels[idx])
            if spin_conventions[idx][ri] == FeynmanDiagram.UpUp
                wuu += root[ri] * phase(varT, extT_labels[idx][ri], ωn[i], β)
            elseif spin_conventions[idx][ri] == FeynmanDiagram.UpDown
                wud += root[ri] * phase(varT, extT_labels[idx][ri], ωn[i], β)
            end
        end
        obs[idx][1, i, n] += wuu * factor * inverse_probability
        obs[idx][2, i, n] += wud * factor * inverse_probability
    end
end

function one_angle_averaged_ParquetAD_Clib(paras::Vector{OneAngleAveraged}, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    measurefreq=5,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    # isLayered2D=false,
    integrand::Function=integrand_OAA_Clib,
    kwargs...
)
    dim, β, kF, order = paras[1].para.dim, paras[1].para.β, paras[1].para.kF, paras[1].para.order
    Nw = length(paras[1].ωn)
    partition, diagpara, extT_labels, spin_conventions = diagram
    for p in paras
        # p.para.isDynamic && UEG.MCinitialize!(p.para)
        if p.para.isDynamic
            if NoBubble in diagpara[1].filter
                UEG.MCinitialize!(p.para, false)
            else
                UEG.MCinitialize!(p.para, true)
            end
            root_dir = joinpath(root_dir, "dynamic/")
        end
        @assert length(p.ωn) == Nw "All parameters must have the same frequency list"
        @assert p.para.dim == dim "All parameters must have the same dimension"
        @assert p.para.β ≈ β "All parameters must have the same inverse temperature"
        @assert p.para.kF ≈ kF "All parameters must have the same Fermi momentum"
    end

    @assert length(diagpara) == length(extT_labels) == length(spin_conventions)

    maxMomNum = maximum([key[1] for key in partition]) + 2

    df = CSV.read(root_dir * "loopBasis_vertex4_maxOrder$(order).csv", DataFrame)
    loopBasis = [df[!, col] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)

    if paras[1].para.isDynamic
        f = jldopen(root_dir * "leafinfo_O$(order).jld2", "r")
        leafstates = f["leafstates"]
        leafvalues = f["values"]
    else
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
    end

    root = zeros(Float64, maximum(length.(extT_labels)))

    # all momentum will be sampled around the dimensionless Fermi momentum 1.0
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3) # the first three momenta are external
    K.data[:, 1] = UEG.getK(kF, dim, 1)
    K.data[:, 2] = UEG.getK(kF, dim, 1)
    K.data[:, 3] = UEG.getK(kF, dim, 1)
    # all time variables will be sampled within [0.0, 1.0]
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha) # the first one is external
    T.data[1] = 0.0
    if dim == 3
        X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    elseif dim == 2
        X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
    end
    N = MCIntegration.Discrete(1, length(paras), alpha=alpha) #index of paras

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nw, length(paras)) for p in diagpara]

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, N),
            dof=dof,
            obs=obs,
            type=Weight,
            measurefreq=measurefreq,
            userdata=(paras, maxMomNum, extT_labels,
                spin_conventions, leafstates, leafvalues, momLoopPool,
                root, partition),
            kwargs...
        )
    end
    result = integrate(integrand; measure=measure_OAA_Clib, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            # report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
            # report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            r = measurement.(real(avg), real(std))
            i = measurement.(imag(avg), imag(std))
            data = Complex.(r, i)
            datadict[key] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function MC_OAA_Clib(para, chan::Symbol;  # chan: :PH or :PP
    kamp=[para.kF,], kamp2=kamp, n=[[-1, 0, 0, -1]], l=0,
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree],
    # channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order),
    transferLoop=nothing,
    verbose=0
)
    diaginfo = Ver4.diagram_loadinfo(para, partition,
        filter=filter, transferLoop=transferLoop)

    partition = diaginfo[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
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

    paras = [Ver4.OneAngleAveraged(para, [kamp[1], kamp2[1]], n, chan, l),]
    ver4, result = Ver4.one_angle_averaged_ParquetAD_Clib(paras, diaginfo;
        neval=neval, print=verbose,
        neighbor=neighbor,
        reweight_goal=reweight_goal
    )
    if isnothing(ver4) == false
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

const dynamic_partition_map = Dict(
    # order of partitions changes, and matters for Clib 
    1 => [(1, 0, 0)],
    2 => [(2, 0, 0), (1, 0, 1), (1, 0, 0)],
    3 => [(1, 0, 2), (2, 0, 1), (2, 0, 0), (2, 1, 0), (1, 0, 1), (3, 0, 0), (1, 0, 0)]
)