function diagramParquet_dk(paramc::ParaMC, _partition::Vector{T}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    extT_labels = Vector{Vector{Int}}[]

    jldopen(joinpath(@__DIR__, "source_codeGV", "extT_ParquetAD_dk.jld2"), "r") do f
        for p in _partition
            push!(diagpara, diagPara(paramc, p[1], filter))
            key_str = join(string.(p))
            push!(extT_labels, f[key_str])
        end
    end
    return (_partition, diagpara, extT_labels)
end

function diagramParquet(paramc::ParaMC, _partition::Vector{T}, spinPolarPara::Float64,
    extra_variables::Dict{String,Int}; filter=[FeynmanDiagram.NoHartree]) where {T}
    diagpara = Vector{DiagParaF64}()
    FeynGraphs = FeynmanDiagram.diagdict_parquet(:sigma, _partition, extra_variables, spinPolarPara=spinPolarPara, filter=filter, isDynamic=paramc.isDynamic)
    extT_labels = Vector{Vector{Int}}[]

    for p in _partition
        push!(diagpara, diagPara(paramc, p[1], filter))
        push!(extT_labels, FeynGraphs[p][2])
    end
    return (_partition, diagpara, FeynGraphs, extT_labels)
end

function integrand_dk(idx, vars, config)
    varK, varT, varN, ExtKidx = vars
    para, kgrid, ngrid, MaxLoopNum, extT_labels = config.userdata[1:5]
    leafval, leafType, leafOrders, leafτ_i, leafτ_o, leafMomIdx = config.userdata[6]
    momLoopPool, root = config.userdata[7:8]
    graphfuncs! = config.userdata[9][idx]
    isLayered2D = config.userdata[end]
    dim, β, me, λ, μ, e0, ϵ0 = para.dim, para.β, para.me, para.mass2, para.μ, para.e0, para.ϵ0

    extidx = ExtKidx[1]
    varK.data[1, 1] = kgrid[extidx]
    FrontEnds.update(momLoopPool, varK.data[:, 1:MaxLoopNum])
    for (i, lftype) in enumerate(leafType[idx])
        if lftype == 0
            continue
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o[idx][i]] - varT[leafτ_i[idx][i]]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            ϵ = dot(kq, kq) / (2me) - μ

            # println(momLoopPool[leafMomIdx[idx][i]])
            dmu_order, dk_order = leafOrders[idx][i][1], leafOrders[idx][i][3]
            order = dmu_order + dk_order
            if order == 0
                leafval[idx][i] = Propagator.green(τ, ϵ, β)
            elseif order == 1
                leafval[idx][i] = Spectral.kernelFermiT_dω(τ, ϵ, β) * (-1)^dmu_order
            elseif order == 2
                leafval[idx][i] = Spectral.kernelFermiT_dω2(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 3
                leafval[idx][i] = Spectral.kernelFermiT_dω3(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 4
                leafval[idx][i] = Spectral.kernelFermiT_dω4(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 5
                leafval[idx][i] = Spectral.kernelFermiT_dω5(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            else
                error("not implemented!")
            end
            if dk_order != 0
                leafval[idx][i] *= kq[1] / me * momLoopPool[leafMomIdx[idx][i]][1]
            end
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx[idx][i])
            order = leafOrders[idx][i][2]
            dk_order = leafOrders[idx][i][3]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[idx][i] = e0^2 / ϵ0 * invK * (λ * invK)^order
                if dk_order != 0
                    leafval[idx][i] *= -2(order + 1) * invK * kq[1] * momLoopPool[leafMomIdx[idx][i]][1]
                end
            elseif dim == 2
                abskq = sqrt(dot(kq, kq))
                if isLayered2D == false
                    invK = 1.0 / (abskq + λ)
                    leafval[idx][i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                    if dk_order != 0
                        leafval[idx][i] *= -(order + 1) * invK * kq[1] / abskq * momLoopPool[leafMomIdx[idx][i]][1]
                    end
                else
                    if order == 0
                        q = sqrt(abskq + 1e-16)
                        invK = 1.0 / q
                        leafval[idx][i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[idx][i] = 0.0 # no high-order counterterms
                    end
                    if dk_order != 0
                        error("not implemented!")
                    end
                end
            else
                error("not implemented!")
            end
        else
            error("this leaftype $lftype not implemented!")
        end
    end

    graphfuncs!(root, leafval[idx])

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function integrand_dk_Clib(idx, vars, config)
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
        elseif lftype == 1 #fermionic 
            τ = varT[leafτ_o] - varT[leafτ_i]
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            ϵ = dot(kq, kq) / (2me) - μ

            dmu_order, dk_order = leafOrders[1], leafOrders[3]
            order = dmu_order + dk_order
            if order == 0
                leafval[i] = Propagator.green(τ, ϵ, β)
            elseif order == 1
                leafval[i] = Spectral.kernelFermiT_dω(τ, ϵ, β) * (-1)^dmu_order
            elseif order == 2
                leafval[i] = Spectral.kernelFermiT_dω2(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 3
                leafval[i] = Spectral.kernelFermiT_dω3(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 4
                leafval[i] = Spectral.kernelFermiT_dω4(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            elseif order == 5
                leafval[i] = Spectral.kernelFermiT_dω5(τ, ϵ, β) * (-1)^dmu_order / factorial(dmu_order)
            else
                error("not implemented!")
            end
            if dk_order != 0
                leafval[i] *= kq[1] / me * momLoopPool[leafMomIdx][1]
            end
        elseif lftype == 2 #bosonic 
            kq = FrontEnds.loop(momLoopPool, leafMomIdx)
            order, dk_order = leafOrders[2], leafOrders[3]
            if dim == 3
                invK = 1.0 / (dot(kq, kq) + λ)
                leafval[i] = e0^2 / ϵ0 * invK * (λ * invK)^order
                if dk_order != 0
                    leafval[i] *= -2(order + 1) * invK * kq[1] * momLoopPool[leafMomIdx][1]
                end
            elseif dim == 2
                abskq = sqrt(dot(kq, kq))
                if isLayered2D == false
                    invK = 1.0 / (abskq + λ)
                    leafval[i] = e0^2 / 2ϵ0 * invK * (λ * invK)^order
                    if dk_order != 0
                        leafval[i] *= -(order + 1) * invK * kq[1] / abskq * momLoopPool[leafMomIdx][1]
                    end
                else
                    if order == 0
                        q = sqrt(abskq + 1e-16)
                        invK = 1.0 / q
                        leafval[i] = e0^2 / 2ϵ0 * invK * tanh(λ * q)
                    else
                        leafval[i] = 0.0 # no high-order counterterms
                    end
                    if dk_order != 0
                        error("not implemented!")
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
    # evalfuncParquetAD_dk_map[group](root, leafval)
    evalfuncParquetAD_map[group](root, leafval)

    n = ngrid[varN[1]]
    weight = sum(root[i] * phase(varT, extT, n, β) for (i, extT) in enumerate(extT_labels[idx]))

    loopNum = config.dof[idx][1]
    factor = 1.0 / (2π)^(dim * loopNum)
    return weight * factor
end

function MC_dk(para; kgrid=[para.kF,], ngrid=[0,], neval=1e6, reweight_goal=nothing,
    spinPolarPara::Float64=0.0, # spin-polarization parameter (n_up - n_down) / (n_up + n_down) ∈ [0,1]
    filename::Union{String,Nothing}=nothing, partition=UEG.partition(para.order), diagtype=:Parquet,
    isLayered2D=false, # whether to use the screened Coulomb interaction in 2D or not 
    isClib=true # whether to use compiled C library to calculate the Feynman diagram weight or not (only for spin-unpolarized case now)
)
    kF = para.kF
    neighbor = UEG.neighbor(partition)

    if isLayered2D
        @assert (para.dim == 2) && diagtype == :GV "Only 2D and GV diagrams supports the tanh screened Coulomb interaction"
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

    if diagtype == :Parquet
        if spinPolarPara == 0.0 && isClib
            diagram = Sigma.diagramParquet_dk(para, partition)
            sigma, result = Sigma.ParquetAD_Clib(para, diagram; isLayered2D=isLayered2D,
                integrand=integrand_dk_Clib, jobname="ParquetAD_dk",
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        else
            diagram = Sigma.diagramParquet(para, partition, spinPolarPara, Dict("k" => 1))
            sigma, result = Sigma.ParquetAD(para, diagram; isLayered2D=isLayered2D,
                integrand=integrand_dk,
                neighbor=neighbor, reweight_goal=reweight_goal,
                kgrid=kgrid, ngrid=ngrid, neval=neval, parallel=:nothread)
        end
    else
        error("unknown diagrams' generated type")
    end

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
            # r, i = real(sigma[(key..., 1)]), imag(sigma[(key..., 1)])
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