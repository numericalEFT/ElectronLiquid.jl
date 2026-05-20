function integrandKT(idx, vars, config)
    varK, varT, ExtTidx, ExtKidx = vars
    para, kgrid, tgrid, maxMomNum, extT_labels = config.userdata[1:5]
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = config.userdata[6]
    momLoopPool, root = config.userdata[7:8]
    graphfuncs! = config.userdata[9][idx]
    isLayered2D = config.userdata[end]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    varT.data[2] = tgrid[ExtTidx[1]]
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
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
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

    graphfuncs!(root, leafval[idx])  # allocations due to run-time variable `idx`

    weight = sum(root[i] for i in eachindex(extT_labels[idx]))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function integrandKT_Clib(idx, vars, config)
    varK, varT, ExtTidx, ExtKidx = vars
    para, kgrid, tgrid, maxMomNum, extT_labels = config.userdata[1:5]
    leafstates, leafval = config.userdata[6][idx], config.userdata[7][idx]
    momLoopPool, root = config.userdata[8:9]
    isLayered2D = config.userdata[10]
    partition = config.userdata[11]

    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    varT.data[2] = tgrid[ExtTidx[1]]

    FrontEnds.update(momLoopPool, varK.data[:, 1:maxMomNum])
    for (i, lfstat) in enumerate(leafstates)
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            order = lforders[1]
            leafval[i] = Propagator.green_derive(τ, ϵ, β, order)
        elseif lftype == 2 #bosonic
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            order = lforders[2]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[i] = e0^2 / ϵ0 * invK * (λ * invK)^order
            elseif dim == 2
                if isLayered2D == false
                    invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                    leafval[i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                else
                    if order == 0
                        q = sqrt(dot(kq, kq) + 1e-16)
                        invK = 1.0 / q
                        leafval[i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[i] = 0.0 # no high-order counterterms
                    end
                end
            else
                error("not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    group = partition[idx]
    evalfuncParquetAD_sigma_map[group](root, leafval)

    weight = sum(root[i] for i in eachindex(extT_labels[idx]))
    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function measureKT(idx, vars, obs, weight, config)
    t = vars[3][1]  #imaginary time
    k = vars[4][1]  #K
    obs[idx][t, k] += weight
end

function KT(para::ParaMC, diagram;
    kgrid=[0.0,],
    tgrid=[para.β - 1e-10,], # must be (0, β)
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    integrand::Function=integrandKT,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Polarization.KT"
    para.isDynamic && UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, FeynGraphs, extT_labels = diagram
    maxMomNum = maximum([key[1] for key in partition]) + 1

    funcGraphs! = Dict{Int,Function}()
    leaf_maps = Vector{Dict{Int,Graph}}()
    for (i, key) in enumerate(partition)
        funcGraphs![i], leafmap = Compilers.compile(FeynGraphs[key])
        push!(leaf_maps, leafmap)
    end
    leafStat, loopbasis = FeynmanDiagram.leafstates(leaf_maps, maxMomNum)
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopbasis)

    println("static compile has finished!")

    root = zeros(Float64, maximum(length.(extT_labels)))
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=2, alpha=alpha)
    T.data[1] = 0.0
    T.data[2] = tgrid[1]
    ExtTidx = MCIntegration.Discrete(1, length(tgrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 2, 1, 1] for p in diagpara] # K, T, ExtTidx, ExtKidx
    # observable of polarization diagram of different permutations
    obs = [zeros(Float64, length(tgrid), length(kgrid)) for o in eachindex(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, ExtTidx, ExtKidx),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, tgrid, maxMomNum, extT_labels, leafStat, momLoopPool, root, funcGraphs!, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measureKT, print=print, neval=neval, solver=solver, kwargs...)

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
            data = measurement.(avg, std)
            datadict[key] = data
        end
        return datadict, result
    else
        return nothing, nothing
    end
end

function KT_Clib(para::ParaMC, diagram_info;
    kgrid=[para.kF,],
    tgrid=[para.β - 1e-10,], # must be (0, β)
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    integrand::Function=integrandKT_Clib,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    name="chargePolar",
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.ParquetAD_Clib"
    para.isDynamic && UEG.MCinitialize!(para)

    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, extT_labels = diagram_info
    maxMomNum = maximum([key[1] for key in partition]) + 1

    df = CSV.read(root_dir * "loopBasis_$(name)_maxOrder6.csv", DataFrame)
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

    root = zeros(Float64, maximum(length.(extT_labels)))
    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=2, alpha=alpha)
    T.data[1] = 0.0
    T.data[2] = tgrid[1]
    ExtTidx = MCIntegration.Discrete(1, length(tgrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 2, 1, 1] for p in diagpara] # K, T, ExtTidx, ExtKidx
    # observable of polarization diagram of different permutations
    obs = [zeros(Float64, length(tgrid), length(kgrid)) for _ in eachindex(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, ExtTidx, ExtKidx),
            dof=dof,
            type=Float64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, tgrid, maxMomNum, extT_labels, leafStat, momLoopPool, root, funcGraphs!, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measureKT, print=print, neval=neval, solver=solver, kwargs...)

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