
diagPara(p::ParaMC) = GenericPara(diagType=SigmaDiag, innerLoopNum=p.order, hasTau=true, loopDim=p.dim, spin=p.spin, firstLoopIdx=2,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, p.isDynamic ? [Instant, Dynamic] : [Instant,]),],  #instant charge-charge interaction
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    # NoFock,
    ]
)

function sigmaDiag(order)
    println("Build the diagrams into an experssion tree ...")

    _partition = partition(order)
    println("Diagram set: ", _partition)

    sigma = Dict()
    for p in _partition
        d = Parquet.sigma(diagPara(p[1])).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if isFock == false
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
    return diagpara, diag, root, extT
end

@inline function phase(varT, extT, l, β)
    tin, tout = varT[extT[1]], varT[extT[2]]
    return exp(1im * π * (2l + 1) / β * (tout - tin))
end
