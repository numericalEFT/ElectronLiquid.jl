module Ver4

using Printf, LinearAlgebra
using ..StaticArrays
using ..Parameters
using ..CompositeGrids
using ..ElectronGas
using ..MCIntegration
using ..Lehmann

using ..FeynmanDiagram
using ..Measurements

# push!(LOAD_PATH, "../common/")
using ..UEG
using ..Propagator

function diagPara(order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, UEG.IsDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return GenericPara(
        diagType=Ver4Diag,
        innerLoopNum=order - 1,
        hasTau=true,
        loopDim=UEG.Dim,
        spin=UEG.Spin,
        firstLoopIdx=4,
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

    ver4 = []
    diagpara = Vector{GenericPara}()
    partition = []
    for p in _partition
        para = diagPara(p[1], filter, KinL - KoutL)
        legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
        d = Parquet.vertex4(para, legK, channel).diagram
        # println("1 $p ", d)
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        # println("2 $p ", d)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        # println("3 $p ", d)
        if isempty(d) == false
            if UEG.IsFock # remove the Fock subdiagrams
                d = DiagTree.removeHatreeFock!(d)
            end
            push!(diagpara, para)
            push!(ver4, d)
            push!(partition, p)
        else
            @warn("partition $p doesn't have any diagram. It will be ignored.")
        end
        # rootuu = [idx for idx in diag.root if diag.node.object[idx].para.response == UpUp] #select the diagram with upup
        # rootud = [idx for idx in diag.root if diag.node.object[idx].para.response == UpDown] #select the diagram with updown
        # #assign the external Tau to the corresponding diagrams
        # extTuu = [diag.node.object[idx].para.extT for idx in rootuu]
        # extTud = [diag.node.object[idx].para.extT for idx in rootud]
    end



    diag = [ExprTree.build(d) for d in ver4]    #experssion tree representation of diagrams 
    rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
    rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
    #assign the external Tau to the corresponding diagrams
    extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
    extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
    return (partition, diagpara, diag, [rootuu, rootud], [extTuu, extTud])
end

mutable struct Weight{T} <: FieldVector{2,T}
    d::T
    e::T
    Weight{T}() where {T} = new{T}(0.0, 0.0)
    Weight(d::T, e::T) where {T} = new{T}(d, e)
end

const Base.zero(::Type{Weight}) = Weight(0.0, 0.0)
const Base.abs(w::Weight) = abs(w.d) + abs(w.e) # define abs(Weight)



# include("common.jl")

# include("ver4_avg.jl")
include("ver4KW.jl")

end