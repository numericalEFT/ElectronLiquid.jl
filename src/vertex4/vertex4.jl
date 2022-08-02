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

function diagPara(para::ParaMC, order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return GenericPara(
        diagType=Ver4Diag,
        innerLoopNum=order - 1,
        hasTau=true,
        loopDim=para.dim,
        spin=para.spin,
        firstLoopIdx=4,
        interaction=inter,
        filter=filter,
        transferLoop=transferLoop
    )
end

function diagram(paramc::ParaMC, _partition::AbstractVector;
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
        para = diagPara(paramc, p[1], filter, KinL - KoutL)
        legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
        d = Parquet.vertex4(para, legK, channel).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)
        if isempty(d) == false
            if paramc.isFock # remove the Fock subdiagrams
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


@inline function phase(varT, extT, ninL, noutL, ninR, β)
    # println(extT)
    tInL, tOutL, tInR, tOutR = varT[extT[INL]], varT[extT[OUTL]], varT[extT[INR]], varT[extT[OUTR]]
    winL, woutL, winR = π * (2ninL + 1) / β, π * (2noutL + 1) / β, π * (2ninR + 1) / β
    woutR = winL + winR - woutL
    return exp(-1im * (tInL * winL - tOutL * woutL + tInR * winR - tOutR * woutR))
end

@inline function phase(varT, extT, n, β)
    # println(extT)
    return phase(varT, extT, n[1], n[2], n[3], β)
end


# include("common.jl")

# include("ver4_avg.jl")
include("ver4KW.jl")
include("ver4_PH_l.jl")
include("exchange_interaction.jl")

end