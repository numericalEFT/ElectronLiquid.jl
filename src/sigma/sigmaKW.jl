diagPara(order::Int) = GenericPara(diagType=SigmaDiag, innerLoopNum=order, hasTau=true, loopDim=UEG.Dim, spin=UEG.Spin, firstLoopIdx=2,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, UEG.IsDynamic ? [Instant, Dynamic] : [Instant,]),],  #instant charge-charge interaction
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    # NoFock,
    ]
)

function sigmaDiag(order)
    println("Build the order $order diagrams into an experssion tree ...")
    _partition = UEG.partition(order)
    println("Diagram set: ", _partition)

    sigma = Dict()
    for p in _partition
        d = Parquet.sigma(diagPara(p[1])).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if UEG.IsFock == false
            sigma[p] = d
        else # remove the Fock subdiagrams
            if p == (1, 0, 0) # the Fock diagram itself should not be removed
                sigma[p] = d
            else
                sigma[p] = DiagTree.removeHatreeFock!(d)
            end
        end
    end

    sigma = [sigma[p] for p in _partition]
    diagpara = [diags[1].id.para for diags in sigma]
    diag = [ExprTree.build(diags) for diags in sigma] # DiagTree to ExprTree
    root = [d.root for d in diag] #get the list of root nodes
    #assign the external Tau to the corresponding diagrams
    extT = [[diag[ri].node.object[idx].para.extT for idx in r] for (ri, r) in enumerate(root)]
    # extT2 = [[diag[ri].node.object[idx].siteidx for idx in r] for (ri, r) in enumerate(root)]
    # println(extT, " verus ", extT2)
    return _partition, diagpara, diag, root, extT
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end

function integrandKW(config)
    para, diag, kgrid, ngrid = config.para
    diagram = diag[config.curr]
    weight = diagram.node.current
    object = diagram.node.object
    l = config.var[3][1]
    k = config.var[4][1]
    varK, varT = config.var[1], config.var[2]
    varK.data[1, 1] = kgrid[k]
    wn = ngrid[l]

    ExprTree.evalKT!(diagram, varK.data, varT.data, para)
    w = sum(weight[r] * phase(varT, object[r].para.extT, wn, para.β) for r in diagram.root)
    return w #the current implementation of sigma has an additional minus sign compared to the standard defintion
end

function measureKW(config)
    factor = 1.0 / config.reweight[config.curr]
    l = config.var[3][1]
    k = config.var[4][1]
    # println(config.observable[1][1])
    o = config.curr
    weight = integrandKW(config)
    config.observable[o, l, k] += weight / abs(weight) * factor
end

function sigmaKW(para::ParaMC;
    kgrid=[para.kF,],
    ngrid=[0,],
    neval=1e6, #number of evaluations
    niter=10, block=16, print=0,
    alpha=3.0, #learning ratio
    reweight_goal=[1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 4.0, 2.0],
    kwargs...
)
    UEG.MCinitialize!(para)

    dim, β, kF = para.dim, para.β, para.kF
    K = MCIntegration.FermiK(dim, kF, 0.5 * kF, 10.0 * kF, offset=1)
    K.data[:, 1] .= 0.0
    K.data[1, 1] = kF
    # T = MCIntegration.Tau(β, β / 2.0, offset=1)
    T = MCIntegration.Continuous(0.0, β, offset=1, alpha=alpha)
    T.data[1] = 0.0
    X = MCIntegration.Discrete(1, length(ngrid), alpha=alpha)
    ExtKidx = MCIntegration.Discrete(1, length(kgrid), alpha=alpha)

    # println("building diagram ...")
    @time partition, diagpara, diag, root, extT = sigmaDiag(para.order)

    dof = [[p.innerLoopNum, p.totalTauNum - 1, 1, 1] for p in diagpara] # K, T, ExtKidx
    obs = zeros(ComplexF64, length(dof), length(ngrid), length(kgrid)) # observable for the Fock diagram 

    ngb = UEG.neighbor(UEG.partition(para.order))
    config = MCIntegration.Configuration((K, T, X, ExtKidx), dof, obs;
        neighbor=ngb, para=(para, diag, kgrid, ngrid),
        reweight_goal=reweight_goal[1:length(dof)+1]
    )

    result = MCIntegration.sample(config, integrandKW, measureKW; neval=neval, niter=niter, print=print, block=block)

    if isnothing(result) == false

        avg, std = result.mean, result.stdev
        if print >= 0
            MCIntegration.summary(result.config)
            println(MCIntegration.summary(result, [o -> real(o[i]) for i in 1:length(dof)]))
        end

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


# function muZ(para::ParaMC; kwargs...)
#     function zfactor(data, β)
#         return @. (imag(data[:, 2, 1]) - imag(data[:, 1, 1])) / (2π / β)
#     end

#     function mu(data)
#         return real(data[:, 1, 1])
#     end

#     data = sigmaKW(para; kgrid=[para.kF,], ngrid=[0, 1], kwargs...)

#     return mu(data), zfactor(data, para.β)
# end

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
