println("Build the diagrams into an experssion tree ...")
# include("../common/eval.jl")

# const Order = 1

diagPara(order) = GenericPara(diagType=Ver4Diag, innerLoopNum=order, hasTau=true, loopDim=UEG.Dim, spin=UEG.Spin, firstLoopIdx=3,
    interaction=[FeynmanDiagram.Interaction(ChargeCharge, [
        Instant,
        Dynamic
    ]),],  #instant charge-charge interaction
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    ],
    transferLoop=KinL - KoutL
)

function ver4Diag(order)
    println("Build the vertex4 diagrams into an experssion tree ...")
    _partition = UEG.partition(order)
    println("Diagram set: ", _partition)

    KinL = KoutL = [1.0, 0, 0]
    KinR = KoutR = [0, 1.0, 0]
    legK = [KinL, KoutL, KinR, KoutR]
    chan = [PHr, PHEr, PPr]

    ver4 = Dict()
    for p in _partition
        d = Parquet.ver4(diagPara(p[1]), legK, chan).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if UEG.IsFock == false
            ver4[p] = d
        else # remove the Fock subdiagrams
            ver4[p] = DiagTree.removeHatreeFock!(d)
        end
    end
end

# const diagpara = [diagPara(o) for o in 1:Order]
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PHr, PHEr, PPr]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PHEr,]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PPr,]) for i in 1:Order]   #diagram of different orders
# ver4 = [Parquet.vertex4(diagpara[i], legK, [PHr,]) for i in 1:Order]   #diagram of different orders
#different order has different set of K, T variables, thus must have different exprtrees
# println(ver4)

# const diag = [ExprTree.build(ver4[o].diagram) for o in 1:Order]    #experssion tree representation of diagrams 
# const rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
# const rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
# #assign the external Tau to the corresponding diagrams
# const extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
# const extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]

# exit(0)


@inline function phase(varT, extT, β, isF)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]],
    varT[extT[OUTR]]
    if (isF)
        return cos(π / β * ((tInL + tOutL) - (tInR + tOutR)))
    else
        return cos(π / β * ((tInL - tOutL) + (tInR - tOutR)))
    end
end