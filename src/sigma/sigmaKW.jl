function integrandKW(config)
    para, diag = config.para
    diagram = diag[config.curr]
    weight = diagram.node.current
    object = diagram.node.object
    l = config.var[3][1]
    varK, varT = config.var[1], config.var[2]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, object[r].para.extT, l, para.β) for r in diagram.root)
    return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measureKW(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[3][1]
    k = config.var[4][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrandKW(config)
    config.observable[o, l+1, k] += weight / abs(weight) * factor
end

function sigmaKW(para::ParaMC;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    niter=10, block=16, print=0,
    reweight_goal=[1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0]
)
    UEG.MCinitialize!(para)
    dim, β, kF = para.dim, para.β, para.kF
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=3.0)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(lgrid[1], lgrid[end], alpha=3.0)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid))

    # println("building diagram ...")
    @time diagpara, diag, root, extT = sigmaDiag(para.order)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), Nl, length(kgrid)) # observable for the Fock diagram 

    ngb = UEG.neighbor(UEG.partition(para.order))
    config = MCIntegration.Configuration((K, T, X, ExtKidx), dof, obs;
        neighbor=ngb, para=(para, diag),
        reweight_goal=reweight_goal[1:length(dof)+1]
    )

    result = MCIntegration.sample(config, integrandKW, measureKW; neval=neval, niter=niter, print=print, block=block)

    if isnothing(result) == false

        # jldsave("data.jld2", order=Order, partition=UEG.partition(Order), avg=avg, std=std)
        avg, std = result.mean, result.stdev
        if print >= 0
            println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))
        end

        # jldopen("dataCT.jld2", "a+") do f
        #     key = "$(UEG.short(para))"
        #     if haskey(f, key)
        #         @warn("replacing existing data for $key")
        #         delete!(f, key)
        #     end
        #     f[key] = (para, avg, std)
        # end
        return avg, std
    end

end

#p = ParaMC(rs=5.0, beta=25.0, Fs=-0.585, order=Order, mass2=1e-5)
#MC(p)
#p = ParaMC(rs=5.0, beta=25.0, Fs=-1.0, order=Order, mass2=1e-5)
#MC(p)
# p = ParaMC(rs=5.0, beta=25.0, Fs=-0.0, order=Order, mass2=0.001)
# MC(p)
# p = ParaMC(rs=5.0, beta=25.0, Fs=-1.0, order=Order, mass2=0.001)
# MC(p)
# p = ParaMC(rs=5.0, beta=25.0, Fs=-1.0, order=Order, mass2=0.01)
# MC(p)
