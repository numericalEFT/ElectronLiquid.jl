function integrand_parquetAD_Clib(idx, vars, config)
    varK, varT = vars
    para, maxMomNum = config.userdata[1:2]
    leafstates, leafval = config.userdata[3][idx], config.userdata[4][idx]
    momLoopPool, root, partition, isLayered2D = config.userdata[5:end]

    dim, β, me, μ = para.dim, para.β, para.me, para.μ
    loopNum = config.dof[idx][1]
    is_zero_order = partition[idx] == (0, 0, 0) ? true : false

    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])
    for (i, lfstat) in enumerate(leafstates)
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
        if lftype == 0
            continue
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            order = lforders[1]
            leafval[i] = Propagator.green_derive(τ, ϵ, β, order)
            if is_zero_order
                leafval[i] *= (ϵ + μ)
            end
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            τ2, τ1 = varT[leafτ_o], varT[leafτ_i]
            leafval[i] = Propagator.interaction_derive(τ1, τ2, kq, para, lforders;
                idtype=Instant, tau_num=1, isLayered=isLayered2D)
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    group = partition[idx]
    evalfunc_freeEnergy_map[group](root, leafval)

    factor = 1.0 / (2π)^(dim * loopNum)
    return root[1] * factor
end

function ParquetAD_Clib(para::ParaMC, diagram_info;
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    integrand::Function=integrand_parquetAD_Clib,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    name="freeEnergy",
    kwargs...
)
    MaxOrder = 5
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.GV"
    if para.isDynamic
        if NoBubble in diagpara[1].filter
            UEG.MCinitialize!(para, false)
        else
            UEG.MCinitialize!(para, true)
        end
    end
    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara = diagram_info

    maxMomNum = maximum([key[1] for key in partition]) + 1

    df = CSV.read(root_dir * "loopBasis_$(name)_maxOrder$(MaxOrder).csv", DataFrame)
    loopBasis = [df[!, col][1:maxMomNum] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)

    leafstates = Vector{Vector{LeafStateAD}}()
    leafvalues = Vector{Vector{Float64}}()
    for key in partition
        key_str = join(string.(key))
        df = CSV.read(root_dir * "leafinfo_$(name)_$key_str.csv", DataFrame)
        leafstates_par = Vector{LeafStateAD}()
        for row in eachrow(df)
            push!(leafstates_par, LeafStateAD(row[2], _StringtoIntVector(row[3]), row[4:end]...))
        end
        push!(leafstates, leafstates_par)
        push!(leafvalues, df[!, names(df)[1]])
    end

    root = zeros(Float64, 1)
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 100.0 * kF)
    T = Continuous(0.0, β; alpha=alpha, adapt=true)

    dof = [[p.innerLoopNum, p.totalTauNum] for p in diagpara] # K, T
    obs = [zero(Float64) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, maxMomNum, leafstates, leafvalues, momLoopPool, root, partition, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measure, print=print, neval=neval, solver=solver, kwargs...)

    if isnothing(result) == false
        if print >= 0
            report(result.config)
            println(report(result, pick=o -> first(o)))
            println(result)
        end

        datadict = Dict{eltype(partition),Any}()
        for (o, key) in enumerate(partition)
            avg, std = result.mean[o], result.stdev[o]
            if key[1] == 0
                datadict[key] = measurement(avg, std)
            else
                datadict[key] = measurement(avg, std) / β
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

@inline function _StringtoIntVector(str::AbstractString)
    pattern = r"[-+]?\d+"
    return [parse(Int, m.match) for m in eachmatch(pattern, str)]
end