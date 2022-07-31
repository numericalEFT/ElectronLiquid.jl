module Ver4

using Printf, LinearAlgebra
using ..Parameters
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
using ..Measurements

# push!(LOAD_PATH, "../common/")
using ..UEG

function diagPara(order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, UEG.IsDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return GenericPara(
        diagType=Ver4Diag,
        innerLoopNum=order,
        hasTau=true,
        loopDim=UEG.Dim,
        spin=UEG.Spin,
        firstLoopIdx=3,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

function diagram(_partition;
    channel=[PHr, PHEr, PPr],
    filter=[
    # Girreducible,
    # Proper,   #one interaction irreduble diagrams or not
    # NoBubble, #allow the bubble diagram or not
    ]
)
    println("Build the vertex4 diagrams into an experssion tree ...")
    # _partition = UEG.partition(order)
    println("Diagram set: ", _partition)

    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0
    legK = [KinL, KoutL, KinR]

    ver4 = Dict()
    diagpara = Vector{GenericPara}()
    for p in _partition
        para = diagPara(p[1], filter, KinL - KoutL)
        push!(diagpara, para)
        legK = [DiagTree.getK(para.totalLoopNum, 1), DiagTree.getK(para.totalLoopNum, 2), DiagTree.getK(para.totalLoopNum, 3)]
        d = Parquet.vertex4(para, legK, channel).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if UEG.IsFock == false
            ver4[p] = d
        else # remove the Fock subdiagrams
            ver4[p] = DiagTree.removeHatreeFock!(d)
        end
    end

    diag = [ExprTree.build(ver4[p]) for p in _partition]    #experssion tree representation of diagrams 
    rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
    rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
    #assign the external Tau to the corresponding diagrams
    extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
    extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
    return (diagpara, diag, [rootuu, rootud], [extTuu, extTud])
end


# include("common.jl")

# include("ver4_avg.jl")

end