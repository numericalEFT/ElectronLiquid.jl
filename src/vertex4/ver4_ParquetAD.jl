
function diagramParquet_load(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{FeynmanDiagram.Response}[]
    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0
    partition = Tuple{Int64,Int64,Int64}[]
    # partition = []
    jldopen(joinpath(@__DIR__, "source_codeParquetAD", "extT_spin_O$(paramc.order)_ParquetAD.jld2"), "r") do f
        for p in _partition
            key_str = join(string.(p))
            if key_str in keys(f)
                extT, spin = f[key_str]
                push!(partition, p)
                push!(diagpara, diagPara(paramc, p[1], filter, KinL - KoutL))
                push!(extT_labels, extT)
                push!(spin_conventions, spin)
            end
        end
    end
    return (partition, diagpara, extT_labels, spin_conventions)
end

function diagramParquet(paramc::ParaMC, _partition::Vector{T};
    channel=[PHr, PHEr, PPr],
    filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0

    FeynGraphs = FeynmanDiagram.GV.diagdict_parquet_ver4(:vertex4,
        _partition, channel=channel,
        filter=filter, isDynamic=paramc.isDynamic, koffset=3)
    extT_labels = Vector{Vector{Int}}[]
    spin_conventions = Vector{FeynmanDiagram.Response}[]
    partition = [k for k in keys(FeynGraphs)]
    for p in partition
        push!(diagpara, diagPara(paramc, p[1], filter, KinL - KoutL))
        push!(extT_labels, FeynGraphs[p][2])
        push!(spin_conventions, FeynGraphs[p][3])
    end
    return (partition, diagpara, FeynGraphs, extT_labels, spin_conventions)
end

function green_derive(τ, ϵ, β, order)
    if order == 0
        result = Propagator.green(τ, ϵ, β)
    elseif order == 1
        result = -Spectral.kernelFermiT_dω(τ, ϵ, β)
    elseif order == 2
        result = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
    elseif order == 3
        result = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
    elseif order == 4
        result = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
    elseif order == 5
        result = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
    else
        error("not implemented!")
        # result = Propagator.green(τ, ϵ, β) * 0.0
    end
    return result
end

function integrand_ParquetAD(idx, var, config)
    weight, factor = diagram_weight_ParquetAD(idx, var, config)
    return weight
end

function diagram_weight_ParquetAD(idx, var, config)
    paras, MaxLoopNum, extT_labels, spin_conventions, leafStat, momLoopPool, root, funcGraphs! = config.userdata
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
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

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    # println(momLoopPool)
    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ
            order = leafOrders[idx][i][1]
            leafval[idx][i] = green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            # println(kq, (k1, k2, x))
            # @assert dot(kq, kq) ≈ (k1^2 + k2^2 - 2k1 * k2 * x) "$(dot(kq, kq)) != $(k1^2 + k2^2 - 2k1 * k2 * x)"
            # order = lftype / 2 - 1
            order = leafOrders[idx][i][2]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[idx][i] = e0^2 / ϵ0 * invK * (λ * invK)^order
            elseif dim == 2
                if isLayered2D == false
                    invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                    leafval[idx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                else
                    if order == 0
                        q = sqrt(dot(kq, kq) + 1e-16)
                        invK = 1.0 / q
                        leafval[idx][i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[idx][i] = 0.0 # no high-order counterterms
                    end
                end
            else
                error("not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs! = funcGraphs![idx]
    graphfuncs!(root, leafval[idx])
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

function measure_ParquetAD(idx, var, obs, relative_weight, config)
    w, factor = diagram_weight_ParquetAD(idx, var, config)
    inverse_probability = abs(relative_weight) / abs(w)
    paras, MaxLoopNum, extT_labels, spin_conventions, leafStat, momLoopPool, root, graphfuncs! = config.userdata

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


function one_angle_averaged_ParquetAD(paras::Vector{OneAngleAveraged}, diagram;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    measurefreq=5,
    kwargs...
)

    dim, β, kF = paras[1].para.dim, paras[1].para.β, paras[1].para.kF
    Nw = length(paras[1].ωn)
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram
    for p in paras
        # p.para.isDynamic && UEG.MCinitialize!(p.para)
        if p.para.isDynamic
            if NoBubble in diagpara[1].filter
                UEG.MCinitialize!(p.para, false)
            else
                UEG.MCinitialize!(p.para, true)
            end
        end
        @assert length(p.ωn) == Nw "All parameters must have the same frequency list"
        @assert p.para.dim == dim "All parameters must have the same dimension"
        @assert p.para.β ≈ β "All parameters must have the same inverse temperature"
        @assert p.para.kF ≈ kF "All parameters must have the same Fermi momentum"
    end

    @assert length(diagpara) == length(FeynGraphs) == length(extT_labels) == length(spin_conventions)

    MaxLoopNum = maximum([key[1] for key in partition]) + 2
    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()

    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key][1])
        push!(leaf_maps, leafmap)
    end
    leafStat, loopbasis = FeynmanDiagram.leafstates_diagtree(leaf_maps, MaxLoopNum)
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopbasis)

    println("static compile has finished!")

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

    # if isnothing(neighbor)
    #     neighbor = UEG.neighbor(partition)
    # end
    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, N),
            dof=dof,
            obs=obs,
            type=Weight,
            measurefreq=measurefreq,
            userdata=(paras, MaxLoopNum, extT_labels,
                spin_conventions, leafStat, momLoopPool,
                root, funcGraphs!),
            kwargs...
        )
    end
    result = integrate(integrand_ParquetAD; measure=measure_ParquetAD, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

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

