"""
Calculate exchange vertex4 averged on the Fermi surface
"""
# const lgrid = [1, 2]
# const Nl = length(lgrid)

@inline function phase(varT, extT, ninL, noutL, ninR, β)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]], varT[extT[OUTR]]
    winL, woutL, winR = π * (2ninL + 1) / β, π * (2noutL + 1) / β, π * (2ninR + 1) / β
    woutR = winL + winR - woutL
    return exp(-1im * (tInL * winL - tOutL * woutL + tInR * winR - tOutR * woutR))
end

function integrandKW(config)
    para, diag, root, extT, kinL, koutL, kinR, ninL, noutL, ninR = config.para

    kF, β = para.kF, para.β
    idx = config.curr
    var = config.var
    varK, varT = var[1], var[2]
    _kinL, _koutL, _kinR = kinL[var[3][1]], koutL[var[4][1]], kinR[var[5][1]]
    _ninL, _noutL, _ninR = ninL[var[6][1]], noutL[var[7][1]], ninR[var[8][1]]

    # println(kinL, ", ", koutL, ", ", kinR)
    varK.data[:, 1] .= _kinL
    varK.data[:, 2] .= _koutL
    varK.data[:, 3] .= _kinR

    diagram = diag[idx]
    weight = diagram.node.current
    rootuu, rootud = root[1][idx], root[2][idx]
    extTuu, extTud = extT[1][idx], extT[2][idx]

    # varK.data[:, 2] = [kF * x, kF * sqrt(1 - x^2), 0.0]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    if !isempty(rootuu)
        wuu = sum(weight[root] * phase(varT, extTuu[ri], _ninL, _noutL, _ninR, β) for (ri, root) in enumerate(rootuu))
    else
        wuu = zero(ComplexF64)
    end
    if !isempty(rootud)
        wud = sum(weight[root] * phase(varT, extTud[ri], _ninL, _noutL, _ninR, β) for (ri, root) in enumerate(rootud))
    else
        wud = zero(ComplexF64)
    end
    # println(_kinL, ", ", _koutL, ", ", _kinR)
    # println(weight[3], ", ", weight[4])
    # println(idx)
    # println(wuu, ",  ", wud)
    # error("stop")
    # return Weight(wuu * para.NF, wud * para.NF)
    return Weight(zero(ComplexF64), wud * para.NF)
end

function measureKW(config)
    factor = 1.0 / config.reweight[config.curr]
    var = config.var
    kinL, koutL, kinR = var[3][1], var[4][1], var[5][1]
    ninL, noutL, ninR = var[6][1], var[7][1], var[8][1]
    # println(config.observable[1][1])
    # if config.curr == 1
    o = config.curr
    weight = integrandKW(config)
    config.observable[o, 1, kinL, koutL, kinR, ninL, noutL, ninR] += weight.d / abs(weight) * factor
    config.observable[o, 2, kinL, koutL, kinR, ninL, noutL, ninR] += weight.e / abs(weight) * factor
    # else
    #     return
    # end
end

function KW(para::ParaMC, diagram;
    kinL=[[para.kF, 0.0, 0.0][1:para.dim],],
    ninL=[0,],
    koutL=[[para.kF, 0.0, 0.0][1:para.dim],],
    noutL=[0,],
    kinR=[[para.kF, 0.0, 0.0][1:para.dim],],
    ninR=[0,],
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

    partition, diagpara, diag, root, extT = diagram
    @assert length(diagpara) == length(diag) == length(root[1]) == length(extT[1])
    @assert length(root[1]) == length(root[2])
    @assert length(extT[1]) == length(extT[2])

    K = MCIntegration.FermiK(para.dim, kF, 0.2 * kF, 10.0 * kF, offset=3)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    # X = MCIntegration.Continuous(-1.0, 1.0) #x=cos(θ)

    NkinL, NwinL, NkoutL, NwoutL, NkinR, NwinR = length(kinL), length(ninL), length(koutL), length(noutL), length(kinR), length(ninR)

    vKinL = MCIntegration.Discrete(1, NkinL, alpha=alpha)
    vWinL = MCIntegration.Discrete(1, NwinL, alpha=alpha)
    vKoutL = MCIntegration.Discrete(1, NkoutL, alpha=alpha)
    vWoutL = MCIntegration.Discrete(1, NwoutL, alpha=alpha)
    vKinR = MCIntegration.Discrete(1, NkinR, alpha=alpha)
    vWinR = MCIntegration.Discrete(1, NwinR, alpha=alpha)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1, 1, 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), 2, NkinL, NkoutL, NkinR, NwinL, NwoutL, NwinR) # observable for the Fock diagram 

    if isnothing(neighbor)
        neighbor = UEG.neighbor(UEG.partition(para.order))
    end
    if isnothing(config)
        config = MCIntegration.Configuration((K, T, vKinL, vKoutL, vKinR, vWinL, vWoutL, vWinR), dof, obs;
            para=(para, diag, root, extT, kinL, koutL, kinR, ninL, noutL, ninR),
            neighbor=neighbor, reweight_goal=reweight_goal
        )
    end
    result = MCIntegration.sample(config, integrandKW, measureKW; neval=neval, niter=niter, print=print, block=block)

    # function info(idx, di)
    #     return @sprintf("   %8.4f ±%8.4f", avg[idx, di], std[idx, di])
    # end

    if isnothing(result) == false
        # avg, std = result.mean, result.stdev
        # avg *= para.NFstar
        # std *= para.NFstar
        # N = size(avg)[1]
        # grid = lgrid
        if print >= 0
            # println(MCIntegration.summary(result, [o -> real(o[i, 1, 1, 1, 1, 1, 1, 1]) for i in 1:length(dof)]))
            MCIntegration.summary(result.config)
            MCIntegration.summary(result, [o -> (real(o[i, 1, 1, 1, 1, 1, 1, 1])) for i in 1:length(dof)], ["uu$i" for i in 1:length(dof)])
            MCIntegration.summary(result, [o -> (real(o[i, 2, 1, 1, 1, 1, 1, 1])) for i in 1:length(dof)], ["ud$i" for i in 1:length(dof)])
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
        return result
    end

end

# p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=Order, mass2=0.001)
# MC(p)