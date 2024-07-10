@inline function phase_ver3(varT, extT, nqout, nkin, β)
    # println(extT)
    # tq, tkin, tkout = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    tkout, tkin, tq = varT[extT[1]], varT[extT[2]], varT[extT[3]]
    wqout, wkin = π * (2nqout + 1) / β, π * (2nkin + 1) / β
    wkout = wkin - wqout
    return exp(-1im * (tkin * wkin - tq * wqout - tkout * wkout))
end

# @inline function phase_ver3(varT, extT, n, β)
#     # println(extT)
#     return phase_ver3(varT, extT, n[1], n[2], β)
# end
@inline function interactionTauNum(type::AnalyticProperty)
    if type == Instant
        return 1
    else
        return 2
    end
end

"""
integrand of vertex3
"""
function integrand_ver3KW(idx, var, config)
    para, kin, nkin, qout, nqout = config.userdata[1:5]
    maxMomNum, extT_labels, spin_conventions, leafStat, leaf_maps, momLoopPool, root, funcGraphs! = config.userdata[6:end]

    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = leafStat
    varK, varT = var[1], var[2]
    loopNum = config.dof[idx][1]
    k1 = kin[var[3][1]]
    n1 = nkin[var[4][1]]
    q = qout[var[5][1]]
    nq = nqout[var[6][1]]

    varK.data[:, 1] .= q
    varK.data[1, 2] = k1

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

    # factor = para.NF / (2π)^(dim * (loopNum))
    factor = 1.0 / (2π)^(dim * (loopNum))
    graphfuncs! = funcGraphs![idx]
    graphfuncs!(root, leafval[idx])
    wuu = zero(ComplexF64)
    wud = zero(ComplexF64)
    for ri in 1:length(extT_labels[idx])
        if spin_conventions[idx][ri] == UpUp
            wuu += root[ri] * phase_ver3(varT, extT_labels[idx][ri], nq, n1, β)
        elseif spin_conventions[idx][ri] == UpDown
            wud += root[ri] * phase_ver3(varT, extT_labels[idx][ri], nq, n1, β)
        end
    end

    return Weight(wuu * factor, wud * factor)

    # return Weight(1.0, 1.0)
end

function measure_ver3KW(idx, var, obs, weight, config)
    KINidx = var[3][1]
    NKINidx = var[4][1]
    QOUTidx = var[5][1]
    NQOUTidx = var[6][1]
    obs[idx][1, KINidx, NKINidx, QOUTidx, NQOUTidx] += weight.d
    obs[idx][2, KINidx, NKINidx, QOUTidx, NQOUTidx] += weight.e
end

function KW(para::ParaMC, diagram;
    kin=[para.kF,], #amplitude of kin
    nkin=[0,], # matfreq of kin
    qout=[getK(0.0, para.dim, 1),],
    nqout=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    integrand::Function=integrand_ver3KW,
    kwargs...)

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    partition, diagpara, FeynGraphs, extT_labels, spin_conventions = diagram

    if NoBubble in diagpara[1].filter
        UEG.MCinitialize!(para, false)
    else
        UEG.MCinitialize!(para, true)
    end

    for p in diagpara
        @assert diagpara[1].filter == p.filter "filter should be the same"
    end

    @assert length(diagpara) == length(FeynGraphs) == length(extT_labels) == length(spin_conventions)

    Nkin = length(kin)
    Nnkin = length(nkin)
    Nqout = length(qout)
    Nnqout = length(nqout)

    maxMomNum = maximum([key[1] for key in partition]) + 2
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

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=2)
    K.data[:, 1] .= UEG.getK(kin[1], dim, 1)
    K.data[:, 2] .= UEG.getK(kin[1], dim, 1) .- qout[1]
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0

    KIN = MCIntegration.Discrete(1, Nkin, alpha=alpha)
    NKIN = MCIntegration.Discrete(1, Nnkin, alpha=alpha)
    QOUT = MCIntegration.Discrete(1, Nqout, alpha=alpha)
    NQOUT = MCIntegration.Discrete(1, Nnqout, alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nkin, Nnkin, Nqout, Nnqout) for p in diagpara]
    println("obs size:", size(obs[1]))

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, KIN, NKIN, QOUT, NQOUT),
            dof=dof,
            obs=obs,
            type=Weight,
            # type=ComplexF64, # type of the integrand
            userdata=(para, kin, nkin, qout, nqout, maxMomNum, extT_labels,
                spin_conventions, leafStat, leaf_maps, momLoopPool,
                root, funcGraphs!),
            kwargs...
        )
    end
    result = integrate(integrand; measure=measure_ver3KW, config=config, solver=solver, neval=neval, print=print, kwargs...)

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

function MC_KW_angle(para;
    kin=[para.kF,], nkin=[0,],
    qout=[get(0.0, para.dim, 1),], nqout=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree,],
    channels=[PHr, PHEr, PPr, Alli],
    partition=UEG.partition(para.order),
    transferLoop=nothing, extK=nothing, optimize_level=1,
    verbose=0)
    kF = para.kF
    diagram = Diagram.diagram_parquet_response(:vertex3, para, partition,
        channels=channels, filter=filter, extK=extK, transferLoop=transferLoop, optimize_level=optimize_level)
    partition = diagram[1]
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

    ver3, result = Ver3.KW(para, diagram;
        kin=kin, nkin=nkin,
        qout=qout, nqout=nqout,
        neval=neval, print=verbose,
        neighbor=neighbor, reweight_goal=reweight_goal)

    if isnothing(ver3) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (kin, nkin, qout, nqout, ver3)
            end
        end
    end

    return ver3, result
end