"""
Calculate vertex4 averged on the Fermi surface
"""
function integrandPH(config)
    para, diag, root, extT, kamp, lgrid, n = config.para

    kF, β = para.kF, para.β
    idx = config.curr
    var = config.var
    varK, varT = var[1], var[2]
    x = config.var[3][1]
    varK.data[:, 3] = [kamp * x, kamp * sqrt(1 - x^2), 0.0]
    # error("$(varK.data[:, 1])")
    varL = config.var[4][1]
    l = lgrid[varL]

    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)

    factor = 1.0
    if l == 0
        factor = 1.0 / 2
    elseif l == 1
        factor = x / 2.0
    else
        error("not implemented")
    end

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

    return Weight(wuu * para.NF, wud * para.NF)
end

function measurePH(config)
    factor = 1.0 / config.reweight[config.curr]

    varL = config.var[4][1]

    o = config.curr
    weight = integrandPH(config)
    config.observable[o, 1, varL] += weight.d / abs(weight) * factor
    config.observable[o, 2, varL] += weight.e / abs(weight) * factor
end

function PH(para::ParaMC, diagram;
    kamp=para.kF,
    n=[0, 0, 0],
    l=[0,],
    neval=1e6, #number of evaluations
    niter=10, block=16, print=0,
    alpha=3.0, #learning ratio
    reweight_goal=nothing,
    # reweight_goal=[1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0],
    config=nothing,
    neighbor=nothing,
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF, NF = para.dim, para.β, para.kF, para.NF
    Nl = length(l)

    partition, diagpara, diag, root, extT = diagram
    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    K.data[:, 1] .= UEG.getK(kamp, para.dim, 1)
    K.data[:, 2] .= UEG.getK(kamp, para.dim, 1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Continuous(-1.0, 1.0, alpha=alpha) #x=cos(θ)
    L = MCIntegration.Discrete(1, Nl, alpha=alpha) # angular momentum

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), 2, Nl)

    if isnothing(neighbor)
        neighbor = UEG.neighbor(UEG.partition(para.order))
    end
    if isnothing(config)
        config = MCIntegration.Configuration((K, T, X, L), dof, obs;
            para=(para, diag, root, extT, kamp, l, n),
            neighbor=neighbor, reweight_goal=reweight_goal
        )
    end
    result = MCIntegration.sample(config, integrandPH, measurePH; neval=neval, niter=niter, print=print, block=block)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        if print >= 0
            # println(MCIntegration.summary(result, [o -> real(o[i, 1, 1, 1, 1, 1, 1, 1]) for i in 1:length(dof)]))
            MCIntegration.summary(result.config)
            MCIntegration.summary(result, [o -> (real(o[i, 1, 1])) for i in 1:length(dof)], ["uu$i" for i in 1:length(dof)])
            MCIntegration.summary(result, [o -> (real(o[i, 2, 1])) for i in 1:length(dof)], ["ud$i" for i in 1:length(dof)])
        end

        # println("UpUp ver4: ")
        # for li in 1:N
        #     @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 1], std[li, 1])
        # end
        # println("UpDown ver4: ")
        # for li in 1:N
        #     @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], avg[li, 2], std[li, 2])
        # end

        # println("S ver4: ")
        # for li in 1:N
        #     @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] + avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        # end
        # println("A ver4: ")
        # for li in 1:N
        #     @printf("%8.4f   %8.4f ±%8.4f\n", grid[li], (avg[li, 1] - avg[li, 2]) / 2, (std[li, 1] + std[li, 2]) / 2)
        # end

        # jldopen("dataF.jld2", "a+") do f
        #     key = "$(UEG.short(para))"
        #     if haskey(f, key)
        #         @warn("replacing existing data for $key")
        #         delete!(f, key)
        #     end
        #     f[key] = (para, avg, std)
        # end
        avg, std = result.mean, result.stdev
        r = measurement.(real(avg), real(std))
        i = measurement.(imag(avg), imag(std))
        data = Complex.(r, i)
        datadict = Dict{eltype(partition),typeof(data[1, :, :])}()
        for i in 1:length(dof)
            datadict[partition[i]] = data[i, :, :]
        end
        return datadict, result
    end

end

# p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=Order, mass2=0.001)
# MC(p)