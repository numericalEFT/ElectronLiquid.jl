function integrandKW_Reweight(idx,vars,config)
    weight = DiagramWeight(idx,vars,config)
    return weight
end

function DiagramWeight(pidx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafstates, leafval = config.userdata[6], config.userdata[7]
    momLoopPool, root = config.userdata[8:9]
    isLayered2D = config.userdata[10]
    partition = config.userdata[11]
    part_index = config.userdata[12]
    part_list = config.userdata[13]

    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]

    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    if para.isDynamic
        tau_num = 2
    else
        tau_num = 1
    end
    
    idx = part_index[pidx]
    for (i, lfstat) in enumerate(leafstates[idx])
        lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
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
            order = lforders[2]
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

    group = partition[idx]
    evalfuncParquetAD_sigma_map[group](root, leafval[idx])

    n = ngrid[varN[1]]
    wsigma = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))
    loopNum = config.dof[pidx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return wsigma * factor
end

function measureKW_Reweight(pidx, vars, obs, relative_weight, config) # for the mcmc algorithm
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafstates, leafval = config.userdata[6], config.userdata[7]
    momLoopPool, root = config.userdata[8:9]
    isLayered2D = config.userdata[10]
    partition = config.userdata[11]
    part_index = config.userdata[12]
    part_list = config.userdata[13]

    nidx = varN[1]  #matsubara frequency
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0
    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    if para.isDynamic
        tau_num = 2
    else
        tau_num = 1
    end

    weight = DiagramWeight(pidx, vars, config)
    loopNum = config.dof[pidx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    inverse_probability = abs(relative_weight) / abs(weight)
    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])

    for (j,iidx) in enumerate(part_list[pidx])
        for (i, lfstat) in enumerate(leafstates[iidx])
            lftype, lforders, leafτ_i, leafτ_o, leafMomIdx = lfstat.type, lfstat.orders, lfstat.inTau_idx, lfstat.outTau_idx, lfstat.loop_idx
            if lftype == 0
                continue
            # elseif isodd(lftype) #fermionic 
            elseif lftype == 1 #fermionic 
                τ = varT[leafτ_o] - varT[leafτ_i]
                kq = FrontEnds.loop(momLoopPool, leafMomIdx)
                ϵ = dot(kq, kq) / (2me) - μ
                order = lforders[1]
                println(lforders[1])
                leafval[iidx][i] = Propagator.green_derive(τ, ϵ, β, order)
            elseif lftype == 2 #bosonic
                kq = FrontEnds.loop(momLoopPool, leafMomIdx)
                order = lforders[2]
                if dim == 3
                    invK = 1.0 / (dot(kq, kq) + λ)
                    leafval[iidx][i] = e0^2 / ϵ0 * invK * (λ * invK)^order
                elseif dim == 2
                    if isLayered2D == false
                        invK = 1.0 / (sqrt(dot(kq, kq)) + λ)
                        leafval[iidx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                    else
                        if order == 0
                            q = sqrt(dot(kq, kq) + 1e-16)
                            invK = 1.0 / q
                            leafval[iidx][i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                        else
                            leafval[iidx][i] = 0.0 # no high-order counterterms
                        end
                    end
                else
                    error("not implemented!")
                end
            else
                error("this leaftype $lftype not implemented!")
            end
        end
        group = partition[iidx]
        evalfuncParquetAD_sigma_map[group](root, leafval[iidx])
        n = ngrid[nidx]
        wsigma = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[iidx]))
        obs[pidx][j, nidx, extidx] += wsigma * factor * inverse_probability
    end
end


function ParquetAD_Reweight(para::ParaMC, diagram_info;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    solver=:mcmc,
    isLayered2D::Bool=false,
    integrand::Function=integrandKW_Reweight,
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"),
    name="sigma",
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

    df = CSV.read(root_dir * "loopBasis_$(name)_maxOrder7.csv", DataFrame)
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
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kgrid[1]
    # T = MCIntegration.Continuous(0.0, β; grid=collect(LinRange(0.0, β, 1000)), offset=1, alpha=alpha)
    T = Continuous(0.0, β; alpha=alpha, adapt=true, offset=1)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    part_index = [findall(x -> x == (ni, 0, 0), partition)[1] for ni in 1:para.order]
    part_list = [findall(x -> x[1] == ni, partition) for ni in 1:para.order]
    max_part_num = maximum([length([p for p in partition if p[1] == partition[i][1]]) for i in part_index])


    dof = [[diagpara[i].innerLoopNum, diagpara[i].totalTauNum - 1, 1, 1] for i in part_index] # K, T, X, ExtKidx
    # observable of sigma diagram of different permutations
    obs = [zeros(ComplexF64, max_part_num, length(ngrid), length(kgrid)) for _ in 1:length(dof)]

    if isnothing(config)
        config = Configuration(;
            var=(K, T, X, ExtKidx),
            dof=dof,
            type=ComplexF64, # type of the integrand
            obs=obs,
            userdata=(para, kgrid, ngrid, maxMomNum, extT_labels, leafstates, leafvalues, momLoopPool, root, isLayered2D, partition, part_index, part_list),
            kwargs...
        )
    end

    result = integrate(integrand; config=config, measure=measureKW_Reweight, thermal_ratio=0.2, print=print, neval=neval, solver=solver, kwargs...)

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
        # for (o, key) in enumerate(partition)
        #     avg, std = result.mean[o], result.stdev[o]
        #     r = measurement.(real(avg), real(std))
        #     i = measurement.(imag(avg), imag(std))
        #     data = Complex.(r, i)
        #     datadict[key] = data
        # end
        for k in 1:length(dof)
            for (i,iidx) in enumerate(part_list[k])
                avg = result.mean[k][i, :, :]
                std = result.stdev[k][i, :, :]
                r = measurement.(real.(avg), real.(std))
                i = measurement.(imag.(avg), imag.(std))
                data = Complex.(r, i)
                datadict[partition[iidx]] = data
            end
        end
        return datadict, result
    else
        return nothing, nothing
    end
end


function MC_Reweight(para; kgrid=[para.kF,], ngrid=[0], neval=1e6, reweight_goal=nothing,
    # spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order),
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    root_dir=joinpath(@__DIR__, "source_codeParquetAD/"), verbose=-1
)
    @assert para.spin == 2 "Only spin-unpolarized case is supported for compiled C library"
    kF = para.kF
    neighbor = UEG.neighbor(partition)

    if isLayered2D
        @assert (para.dim == 2) "Only 2D systems supports the tanh screened Coulomb interaction"
    end

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            reweight_factor = 2.0^(2order + sOrder + vOrder - 2)
            if (order, sOrder, vOrder) == (1, 0, 0)
                reweight_factor = 4.0
            end
            push!(reweight_goal, reweight_factor)
        end
        push!(reweight_goal, 4.0)
    end

    diaginfo = Sigma.diagram_loadinfo(para, partition, root_dir=root_dir)
    sigma, result = Sigma.ParquetAD_Reweight(para, diaginfo;
        root_dir=root_dir, isLayered2D=isLayered2D,
        neighbor=neighbor, reweight_goal=reweight_goal,
        kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread, print=verbose)

    if isnothing(sigma) == false
        if isnothing(filename) == false
            jldopen(filename, "a+") do f
                key = "$(UEG.short(para))"
                if haskey(f, key)
                    @warn("replacing existing data for $key")
                    delete!(f, key)
                end
                f[key] = (ngrid, kgrid, sigma)
            end
        end
        for (ip, key) in enumerate(partition)
            println("Group ", key)
            @printf("%10s  %10s   %10s   %10s   %10s \n", "q/kF", "real(avg)", "err", "imag(avg)", "err")
            r, i = real(sigma[key]), imag(sigma[key])
            for (in, n) in enumerate(ngrid)
                println("n = $n")
                for (iq, q) in enumerate(kgrid)
                    @printf("%10.6f  %10.6f ± %10.6f   %10.6f ± %10.6f\n", q[1] / kF, r[in, iq].val, r[in, iq].err, i[in, iq].val, i[in, iq].err)
                end
            end
        end
    end
    return sigma, result
end