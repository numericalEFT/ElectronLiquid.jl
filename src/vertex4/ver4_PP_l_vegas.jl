# Base.abs(w::Vector{ComplexF64}) = sqrt(sum(abs(i)^2 for i in w))
Base.abs(w::Vector{Weight}) = sqrt(sum(abs(i)^2 for i in w))

function integrand_PP_vegas(var, _weights, config)
    para, MaxLoopNum, extT_labels, spin_conventions, leafstates, leafval, momLoopPool, root, isLayered2D, partition, part_index, part_list = config.userdata
    for idx in 1:length(part_index)
        w, factor = diagram_weight_PP_vegas(idx, var, config)
        _weights[idx] = w
        # println(_weights)
    end
end

function diagram_weight_PP_vegas(pidx, var, config)
    para, MaxLoopNum, extT_labels, spin_conventions, leafstates, leafval, momLoopPool, root, isLayered2D, partition, part_index, part_list = config.userdata
    idx = part_index[pidx] # count only (pidx,0,0) diagrams
    varK, varT, varX = var[1], var[2], var[3]

    x = varX[1]
    loopNum = config.dof[pidx][1]

    l = para.l
    param = para.para
    dim, β, me, λ, μ, e0, ϵ0 = param.dim, param.β, param.me, param.mass2, param.μ, param.e0, param.ϵ0

    k1, k2 = para.kamp
    varK.data[1, 1], varK.data[1, 3] = k1, -k1
    varK.data[1, 2] = k2 * x
    varK.data[2, 2] = k2 * sqrt(1 - x^2)

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    # println(momLoopPool)
    for (i, lfstat) in enumerate(leafstates[idx])
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx, tau_num = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx, lfstat.tau_num
        if lftype == 0
            continue
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            order = lforders[1]
            leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Instant, tau_num=tau_num)
        elseif lftype == 4 # dynamic bosonic
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Dynamic, tau_num=tau_num)
        else
            error("this leaftype $lftype not implemented!")
        end
        # @assert !(isnan(leafval[idx][i])) "nan at $idx, $i"
    end

    group = (para.para.order, partition[idx]...)
    if param.isDynamic
        evalfuncParquetADDynamic_map[group](root, leafval[idx])
    else
        evalfuncParquetAD_map[group](root, leafval[idx])
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
    # wuu, wud = 1.0, 1.0
    # if l % 2 == 0
    #     ppweight = -wuu + 2wud
    # else
    #     ppweight = wuu
    # end
    return Weight{ComplexF64}(wuu, wud), factor
    # ppweight = abs(-wuu + 2wud) + abs(wuu)
    # return ppweight, factor
end

function measure_PP_vegas(var, obs, relative_weight, config)
    para, MaxLoopNum, extT_labels, spin_conventions, leafstates, leafval, momLoopPool, root, isLayered2D, partition, part_index, part_list = config.userdata
    for pidx in 1:length(part_index)
        w, factor = diagram_weight_PP_vegas(pidx, var, config)
        inverse_probability = abs(relative_weight) / abs(w)

        varK, varT, varX = var[1], var[2], var[3]
        x = varX[1]
        l = para.l
        param = para.para
        dim, β, me, λ, μ, e0, ϵ0 = param.dim, param.β, param.me, param.mass2, param.μ, param.e0, param.ϵ0
        k1, k2 = para.kamp
        varK.data[1, 1], varK.data[1, 3] = k1, -k1
        varK.data[1, 2] = k2 * x
        varK.data[2, 2] = k2 * sqrt(1 - x^2)

        for i in 1:length(part_list[pidx])
            idx = part_list[pidx][i]
            loopNum = config.dof[pidx][1]
            FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
            # println(momLoopPool)
            for (i, lfstat) in enumerate(leafstates[idx])
                lftype, lforders, leafτ_i, leafτ_o, leafMomIdx, tau_num = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx, lfstat.tau_num
                if lftype == 0
                    continue
                elseif lftype == 1 #fermionic 
                    τ = varT[leafτ_o] - varT[leafτ_i]
                    kq = FrontEnds.loop(momLoopPool, leafMomIdx)
                    ϵ = dot(kq, kq) / (2me) - μ
                    order = lforders[1]
                    leafval[idx][i] = Propagator.green_derive(τ, ϵ, β, order)
                elseif lftype == 2 #bosonic 
                    kq = FrontEnds.loop(momLoopPool, leafMomIdx)
                    τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
                    leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Instant, tau_num=tau_num)
                elseif lftype == 4 # dynamic bosonic
                    kq = FrontEnds.loop(momLoopPool, leafMomIdx)
                    τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
                    leafval[idx][i] = Propagator.interaction_derive(τ1, τ2, kq, param, lforders; idtype=Dynamic, tau_num=tau_num)
                else
                    error("this leaftype $lftype not implemented!")
                end
                # @assert !(isnan(leafval[idx][i])) "nan at $idx, $i"
            end
            # eval all (n, i, j) diagrams
            group = (para.para.order, partition[idx]...)
            param = para.para
            if param.isDynamic
                evalfuncParquetADDynamic_map[group](root, leafval[idx])
            else
                evalfuncParquetAD_map[group](root, leafval[idx])
            end
            ωn = para.ωn
            wuu = zero(ComplexF64)
            wud = zero(ComplexF64)
            for ri in 1:length(extT_labels[idx])
                if spin_conventions[idx][ri] == FeynmanDiagram.UpUp
                    wuu += root[ri] * phase(varT, extT_labels[idx][ri], ωn[1], β)
                    # p = phase(varT, extT_labels[idx][ri], ωn[1], β)
                    # @assert !(isnan(p)) "nan appear in p, T=$varT, ωn=$(ωn[1]), extT_label=$(extT_labels[idx][ri])"
                elseif spin_conventions[idx][ri] == FeynmanDiagram.UpDown
                    wud += root[ri] * phase(varT, extT_labels[idx][ri], ωn[1], β)
                end
            end
            # @assert !(isnan(wuu)) "nan appear in wuu"
            # @assert !(isnan(wud)) "nan appear in wud"
            # @assert !(isnan(inverse_probability)) "nan appear in inverse_probability"
            obs[pidx][1, i] += wuu * factor * inverse_probability
            obs[pidx][2, i] += wud * factor * inverse_probability
            println((obs[pidx][1, i], obs[pidx][2, i]))
        end
    end
end

function MC_PP_ParquetAD_vegas(para;
    partition=UEG.partition(para.order),
    filter=[NoHartree, NoBubble],
    # channel=[PHr, PHEr, PPr], # channel is determined at compiling
    kamp=para.kF, kamp2=kamp, n=[0, 1, -1], l=0,
    neval=1e6,
    measurefreq=10,
    alpha=3.0,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    print=0,
    config=nothing,
    isLayered2D=false,
    kwargs...
)
    ver4para = Ver4.OneAngleAveraged(para, [kamp, kamp2], [n,], :PP, l)
    diagram = Ver4.diagramParquet_load(para, partition; filter=filter)

    dim, β, kF, order = para.dim, para.β, para.kF, para.order
    partition, diagpara, extT_labels, spin_conventions = diagram

    if para.isDynamic
        if NoBubble in diagpara[1].filter
            UEG.MCinitialize!(para, false)
        else
            UEG.MCinitialize!(para, true)
        end
        root_dir = joinpath(@__DIR__, "source_codeParquetAD/dynamic/")
    end

    # if isnothing(reweight_goal)
    #     reweight_goal = Float64[]
    #     for o in 1:order
    #         # push!(reweight_goal, 8.0^(order + vOrder - 1))
    #         push!(reweight_goal, 1.0^(o - 1))
    #     end
    #     push!(reweight_goal, 1.0)
    #     println(length(reweight_goal))
    # end

    MaxLoopNum = maximum([key[1] for key in partition]) + 2
    df = CSV.read(root_dir * "loopBasis_ParquetADmaxOrder$(order).csv", DataFrame)
    loopBasis = [df[!, col] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)
    if para.isDynamic
        f = jldopen(root_dir * "leafinfo_O$(order).jld2", "r")
        leafstates = f["leafstates"]
        leafvalues = f["values"]
        for l in leafvalues
            l = l .* 0.0
        end
    else
        leafstates = Vector{Vector{Ver4.LeafStateADVer4Dynamic}}()
        leafvalues = Vector{Vector{Float64}}()
        for key in partition
            key_str = join(string.(key))
            df = CSV.read(root_dir * "leafinfo_O$(order)_Parquet$key_str.csv", DataFrame)
            leafstates_par = Vector{Ver4.LeafStateADVer4Dynamic}()
            for row in eachrow(df)
                push!(leafstates_par, Ver4.LeafStateADVer4Dynamic(row[2], Sigma._StringtoIntVector(row[3]), row[4:end]..., 1))
            end
            push!(leafstates, leafstates_par)
            push!(leafvalues, df[!, names(df)[1]])
        end
    end

    root = zeros(Float64, maximum(length.(extT_labels)))
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3) # the first three momenta are external
    K.data[:, 1] = UEG.getK(kF, dim, 1)
    K.data[:, 2] = UEG.getK(kF, dim, 1)
    K.data[:, 3] = UEG.getK(kF, dim, 1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha) # the first one is external
    T.data[1] = 0.0
    if dim == 3
        X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    elseif dim == 2
        # X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
        error("not implemented yet!")
    end

    part_index = [findall(x -> x == (n, 0, 0), partition)[1] for n in 1:order] # index of (n,0,0) partition
    max_part_num = maximum([length([p for p in partition if p[1] == partition[i][1]]) for i in part_index])
    part_list = [findall(x -> x[1] == n, partition) for n in 1:order]
    println(part_index)
    println(max_part_num)
    println(part_list)

    dof = [[diagpara[i].innerLoopNum, diagpara[i].totalTauNum - 1, 1] for i in part_index] # K, T, X
    obs = [zeros(ComplexF64, 2, max_part_num) for i in part_index] # uu and ud
    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X),
            dof=dof,
            obs=obs,
            type=Weight,
            # type=ComplexF64,
            # type=Float64,
            measurefreq=measurefreq,
            userdata=(ver4para, MaxLoopNum, extT_labels,
                spin_conventions, leafstates, leafvalues, momLoopPool,
                root, isLayered2D, partition, part_index, part_list),
            kwargs...
        )
    end
    # result = integrate(integrand_PP_vegas; measure=measure_PP_vegas, config=config, solver=:vegasmc, neval=neval, print=print, reweight_goal=reweight_goal, inplace=true, kwargs...)
    result = integrate(integrand_PP_vegas; measure=measure_PP_vegas, config=config, solver=:vegasmc, neval=neval, print=print, inplace=true, kwargs...)
    if isnothing(result) == false
        if print >= 0
            report(result.config)
            # report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
            # report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        end

        datadict = Dict{eltype(partition),Any}()
        println(result.mean)
        println(result.stdev)
        for o in 1:order
            println("collect data for order $o")
            println(part_list[o])
            for i in 1:length(part_list[o])
                p = partition[part_list[o][i]]
                avg, std = result.mean[o][:, i], result.stdev[o][:, i]
                rval = measurement.(real(avg), real(std))
                ival = measurement.(imag(avg), imag(std))
                data = Complex.(rval, ival)
                datadict[p] = data
            end
        end
        ver4 = datadict
    else
        ver4 = nothing
    end

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
