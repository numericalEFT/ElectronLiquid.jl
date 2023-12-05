
struct LeafStateAD
    type::Int
    orders::Vector{Int}
    inTau_idx::Int
    outTau_idx::Int
    loop_idx::Int

    function LeafStateAD(type::Int, orders::Vector{Int}, inTau_idx::Int, outTau_idx::Int, loop_idx::Int)
        return new(type, orders, inTau_idx, outTau_idx, loop_idx)
    end
end

function ParquetAD(para::ParaMC, diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.GV"
    para.isDynamic && UEG.MCinitialize!(para)

    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, FeynGraphs, extT_labels = diagram
    MaxLoopNum = maximum([key[1] for key in partition]) + 1

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
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    # T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, X, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, ngrid, MaxLoopNum, extT_labels, leafStat, momLoopPool, root, funcGraphs!, isLayered2D),
            kwargs...
        )
    end

    result = integrate(integrandGV; config=config, measure=measureGV, print=print, neval=neval, solver=solver, kwargs...)

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

function integrandParquetAD_Clib(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafstates, leafval = config.userdata[6][idx], config.userdata[7][idx]
    momLoopPool, root = config.userdata[8:9]
    isLayered2D = config.userdata[10]
    partition = config.userdata[11]

    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lfstat) in enumerate(leafstates)
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
        if lftype == 0
            continue
            # elseif isodd(lftype) #fermionic 
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ
            # order = (lftype - 1) / 2
            order = lforders[1]
            if order == 0
                leafval[i] = Propagator.green(τ, ϵ, β)
            elseif order == 1
                leafval[i] = -Spectral.kernelFermiT_dω(τ, ϵ, β)
            elseif order == 2
                leafval[i] = Spectral.kernelFermiT_dω2(τ, ϵ, β) / 2.0
            elseif order == 3
                leafval[i] = -Spectral.kernelFermiT_dω3(τ, ϵ, β) / 6.0
            elseif order == 4
                leafval[i] = Spectral.kernelFermiT_dω4(τ, ϵ, β) / 24.0
            elseif order == 5
                leafval[i] = -Spectral.kernelFermiT_dω5(τ, ϵ, β) / 120.0
            else
                error("not implemented!")
            end
        elseif lftype == 2 #bosonic
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            # order = lftype / 2 - 1
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
    evalfuncParquetAD_map[group](root, leafval)

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))
    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function ParquetAD_Clib(para::ParaMC, diagram;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    root_dir=joinpath(@__DIR__, "source_codeGV/"),
    kwargs...
)
    @assert solver == :mcmc "Only :mcmc is supported for Sigma.ParquetAD_Clib"
    para.isDynamic && UEG.MCinitialize!(para)

    if isLayered2D
        @assert para.dim == 2 "Only 2D is supported for the tanh screened Coulomb interaction"
    end

    dim, β, kF = para.dim, para.β, para.kF
    partition, diagpara, extT_labels = diagram
    MaxLoopNum = maximum([key[1] for key in partition]) + 1

    df = CSV.read(root_dir * "loopBasis_ParquetADmaxOrder$(para.order).csv", DataFrame)
    loopBasis = [df[!, col] for col in names(df)]
    momLoopPool = FrontEnds.LoopPool(:K, dim, loopBasis)

    leafstates = Vector{Vector{LeafStateAD}}()
    leafvalues = Vector{Vector{Float64}}()
    for key in partition
        key_str = join(string.(key))
        df = CSV.read(root_dir * "leafinfo_Parquet$key_str.csv", DataFrame)
        leafstates_par = Vector{LeafStateAD}()
        for row in eachrow(df)
            push!(leafstates_par, LeafStateAD(row[2], _StringtoIntVector(row[3]), row[4:end]...))
        end
        push!(leafstates, leafstates_par)
        push!(leafvalues, df[!, names(df)[1]])
    end

    root = zeros(Float64, maximum(length.(extT_labels)))
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    # T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, X, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, length(ngrid), length(kgrid)) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, ngrid, MaxLoopNum, extT_labels, leafstates, leafvalues, momLoopPool, root, isLayered2D, partition),
            kwargs...
        )
    end

    result = integrate(integrandParquetAD_Clib; config=config, measure=measureGV, print=print, neval=neval, solver=solver, kwargs...)

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

@inline function _StringtoIntVector(str::AbstractString)
    pattern = r"[-+]?\d+"
    return [parse(Int, m.match) for m in eachmatch(pattern, str)]
end