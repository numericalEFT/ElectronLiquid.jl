"""
Calculate vertex4 averged on the Fermi surface
"""
function integrandPH(idx, var, config)
    para, diag, root, extT, kampgrid, kamp2grid, qgrid, lgrid, n = config.userdata

    dim, β = para.dim, para.β
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    # error("$(varK.data[:, 1])")
    l = lgrid[var[4][1]]
    loopNum = config.dof[idx][1]
    extKidx = var[5][1]
    kamp = kampgrid[extKidx]
    kamp2 = kamp2grid[extKidx]
    qamp = qgrid[extKidx]
    varK.data[1, 1] = kamp
    varK.data[1, 2] = kamp
    varK.data[2, 2] = qamp
    #varK.data[1, 1], varK.data[1, 2] = kF, kF
    if dim == 3
        varK.data[1:2, 3] = [kamp2 * x, kamp2 * sqrt(1 - x^2)]
    else
        varK.data[1:2, 3] = [kamp2 * cos(x), kamp2 * sin(x)]
    end

    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)

    factor = legendfactor(x, l, dim)

    if !isempty(rootuu)
        wuu = factor * sum(weight[root] * phase(varT, extTuu[ri], n, β) for (ri, root) in enumerate(rootuu))
    else
        wuu = zero(ComplexF64)
    end
    if !isempty(rootud)
        wud = factor * sum(weight[root] * phase(varT, extTud[ri], n, β) for (ri, root) in enumerate(rootud))
    else
        wud = zero(ComplexF64)
    end

    factor = para.NF / (2π)^(dim * loopNum)
    return Weight(wuu * factor, wud * factor)
end

function measurePH(idx, var, obs, weight, config)
    Lidx = var[4][1]
    Kidx = var[5][1]

    obs[idx][1, Lidx, Kidx] += weight.d
    obs[idx][2, Lidx, Kidx] += weight.e
end

function PH(para::ParaMC, diagram;
    kamp=[para.kF,], #amplitude of k of the left legs
    kamp2=kamp, #amplitude of k of the right leg
    q=[0.0 for k in kamp],
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    print=0,
    alpha=3.0, #learning ratio
    config=nothing,
    kwargs...
)
    partition, diagpara, diag, root, extT = diagram

    # UEG.MCinitialize!(para)
    if NoBubble in diagpara[1].filter
        UEG.MCinitialize!(para, false)
    else
        UEG.MCinitialize!(para, true)
    end

    for p in diagpara
        @assert diagpara[1].filter == p.filter "filter should be the same"
    end

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    Nl = length(l)
    Nk = length(kamp)

    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    K = MCIntegration.FermiK(dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    K.data[:, 1] .= UEG.getK(kF, dim, 1)
    K.data[:, 2] .= UEG.getK(kF, dim, 1)
    K.data[:, 3] .= 0.0
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    if dim == 3
        X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    elseif dim == 2
        X = MCIntegration.Continuous(0.0, 2π, alpha=alpha) #x=θ
    end
    L = MCIntegration.Discrete(1, Nl, alpha=alpha) # angular momentum
    AMP = MCIntegration.Discrete(1, Nk, alpha=alpha) # angular momentum

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = [zeros(ComplexF64, 2, Nl, Nk) for p in diagpara]

    if isnothing(config)
        config = MCIntegration.Configuration(;
            var=(K, T, X, L, AMP),
            dof=dof,
            obs=obs,
            type=Weight,
            userdata=(para, diag, root, extT, kamp, kamp2, q, l, n),
            kwargs...
        )
    end
    result = integrate(integrandPH; measure=measurePH, config=config, solver=:mcmc, neval=neval, print=print, kwargs...)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        # if print >= 0
        #     report(result.config)
        #     report(result; pick=o -> (real(o[1, 1, 1])), name="uu")
        #     report(result; pick=o -> (real(o[2, 1, 1])), name="ud")
        # end

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

function MC_PH(para; kamp=[para.kF,], kamp2=kamp, q=[0.0 for k in kamp], n=[-1, 0, 0, -1], l=[0,],
    neval=1e6, filename::Union{String,Nothing}=nothing, reweight_goal=nothing,
    filter=[NoHartree, NoBubble, Proper],
    channel=[PHr, PHEr, PPr],
    partition=UEG.partition(para.order),
    verbose=0
)

    kF = para.kF
    _order = para.order

    # partition = UEG.partition(_order)


    diagram = Ver4.diagram(para, partition; channel=channel, filter=filter)

    partition = diagram[1] # diagram like (1, 1, 0) is absent, so the partition will be modified
    neighbor = UEG.neighbor(partition)

    if isnothing(reweight_goal)
        reweight_goal = Float64[]
        for (order, sOrder, vOrder) in partition
            # push!(reweight_goal, 8.0^(order + vOrder - 1))
            push!(reweight_goal, 8.0^(order - 1))
        end
        push!(reweight_goal, 1.0)
        println(length(reweight_goal))
    end

    ver4, result = Ver4.PH(para, diagram;
        kamp=kamp, kamp2=kamp2, q=q, n=n, l=l,
        neval=neval, print=verbose,
        neighbor=neighbor,
        reweight_goal=reweight_goal
    )

    if isnothing(ver4) == false
        for (p, data) in ver4
            printstyled("permutation: $p\n", color=:yellow)
            for (li, _l) in enumerate(l)
                printstyled("l = $_l\n", color=:green)
                @printf("%12s    %16s    %16s    %16s    %16s    %16s    %16s\n", "k/kF", "uu", "ud", "di", "ex", "symmetric", "asymmetric")
                for (ki, k) in enumerate(kamp)
                    factor = 1.0
                    d1, d2 = real(data[1, li, ki]) * factor, real(data[2, li, ki]) * factor
                    s, a = (d1 + d2) / 2.0, (d1 - d2) / 2.0
                    di, ex = (s - a), (a) * 2.0
                    @printf("%12.6f    %16s    %16s    %16s    %16s    %16s    %16s\n", k / kF, "$d1", "$d2", "$di", "$ex", "$s", "$a")
                end
            end
        end

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
