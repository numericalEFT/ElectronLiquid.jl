function integrand_ParquetAD_Clib(idx, var, config)
    weight, factor = diagram_weight_ParquetAD_Clib(idx, var, config)
    return weight
end

function diagram_weight_ParquetAD_Clib(idx, var, config)
    paras, maxMomNum, extT_labels, spin_conventions, leafstates, leafval, momLoopPool, root, isLayered2D, partition = config.userdata
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

    group = (para.para.order, partition[idx]...)
    if param.isDynamic
        evalfuncParquetADDynamic_map[group](root, leafval[idx])
    else
        evalfuncParquetAD_map[group](root, leafval[idx])
    end
    # get_eval_func(para.para.order, partition[idx])(root, leafval[idx])
    # graphfuncs! = funcGraphs![idx]
    # graphfuncs!(root, leafval[idx])
    # println(root)

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

function measure_ParquetAD_Clib(idx, var, obs, relative_weight, config)
    w, factor = diagram_weight_ParquetAD_Clib(idx, var, config)
    inverse_probability = abs(relative_weight) / abs(w)
    paras, maxMomNum, extT_labels, spin_conventions, leafstates, leafvalues, momLoopPool, root, graphfuncs! = config.userdata

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
    isLayered2D=false,
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
            root_dir = joinpath(@__DIR__, "source_codeParquetAD/dynamic/")
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
        leafstates = Vector{Vector{Ver4.LeafStateADVer4Dynamic}}()
        leafvalues = Vector{Vector{Float64}}()
        for key in partition
            key_str = join(string.(key))
            df = CSV.read(root_dir * "leafinfo_vertex4_$key_str.csv", DataFrame)
            leafstates_par = Vector{Ver4.LeafStateADVer4Dynamic}()
            for row in eachrow(df)
                push!(leafstates_par, Ver4.LeafStateADVer4Dynamic(row[2], _StringtoIntVector(row[3]), row[4:end]..., 1))
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
                root, isLayered2D, partition),
            kwargs...
        )
    end
    result = integrate(integrand_ParquetAD_Clib; measure=measure_ParquetAD_Clib, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

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

# function get_eval_func(order, partition)
#     # slow
#     key_str = join(string.(partition))
#     func_name = "eval_ver4O$(order)ParquetAD$(key_str)!"
#     return getfield(Ver4, Symbol(func_name))
# end

# include("source_codeParquetAD/Cwrapper_ver4O1ParquetAD.jl")
# include("source_codeParquetAD/Cwrapper_ver4O2ParquetAD.jl")
include("source_codeParquetAD/Cwrapper_ver4O3ParquetAD.jl")
include("source_codeParquetAD/Cwrapper_ver4O4ParquetAD.jl")
# include("source_codeParquetAD/Cwrapper_ver4O5ParquetAD.jl")
# include("source_codeParquetAD/Cwrapper_ver4O6ParquetAD.jl")

# include("source_codeParquetAD/dynamic/Cwrapper_ver4O1ParquetADDynamic.jl")
# include("source_codeParquetAD/dynamic/Cwrapper_ver4O2ParquetADDynamic.jl")
# include("source_codeParquetAD/dynamic/Cwrapper_ver4O3ParquetADDynamic.jl")

# provide dict of (order, partition...) => func
include("source_codeParquetAD/func_dict_ParquetAD.jl")
# include("source_codeParquetAD/dynamic/func_dict_ParquetADDynamic.jl")

const evalfuncParquetAD_vertex4_map = Dict(
    (4, 0, 0, 0) => eval_ver4O4ParquetAD000!,
    (4, 0, 0, 1) => eval_ver4O4ParquetAD001!,
    (4, 0, 0, 2) => eval_ver4O4ParquetAD002!,
    (4, 0, 0, 3) => eval_ver4O4ParquetAD003!,
    (4, 0, 0, 4) => eval_ver4O4ParquetAD004!,
    (4, 1, 0, 0) => eval_ver4O4ParquetAD100!,
    (4, 1, 0, 1) => eval_ver4O4ParquetAD101!,
    (4, 1, 0, 2) => eval_ver4O4ParquetAD102!,
    (4, 1, 0, 3) => eval_ver4O4ParquetAD103!,
    (4, 1, 1, 0) => eval_ver4O4ParquetAD110!,
    (4, 1, 1, 1) => eval_ver4O4ParquetAD111!,
    (4, 1, 1, 2) => eval_ver4O4ParquetAD112!,
    (4, 1, 2, 0) => eval_ver4O4ParquetAD120!,
    (4, 1, 2, 1) => eval_ver4O4ParquetAD121!,
    (4, 1, 3, 0) => eval_ver4O4ParquetAD130!,
    (4, 2, 0, 0) => eval_ver4O4ParquetAD200!,
    (4, 2, 0, 1) => eval_ver4O4ParquetAD201!,
    (4, 2, 0, 2) => eval_ver4O4ParquetAD202!,
    (4, 2, 1, 0) => eval_ver4O4ParquetAD210!,
    (4, 2, 1, 1) => eval_ver4O4ParquetAD211!,
    (4, 2, 2, 0) => eval_ver4O4ParquetAD220!,
    (4, 3, 0, 0) => eval_ver4O4ParquetAD300!,
    (4, 3, 0, 1) => eval_ver4O4ParquetAD301!,
    (4, 3, 1, 0) => eval_ver4O4ParquetAD310!,
    (4, 4, 0, 0) => eval_ver4O4ParquetAD400!
)

const dynamic_partition_map = Dict(
    # order of partitions changes, and matters for Clib 
    1 => [(1, 0, 0)],
    2 => [(2, 0, 0), (1, 0, 1), (1, 0, 0)],
    3 => [(1, 0, 2), (2, 0, 1), (2, 0, 0), (2, 1, 0), (1, 0, 1), (3, 0, 0), (1, 0, 0)]
)