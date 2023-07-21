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
import ..ExprTreeF64

import ..Weight

function diagPara(para::ParaMC, order, filter, transferLoop)
    inter = [FeynmanDiagram.Interaction(ChargeCharge, para.isDynamic ? [Instant, Dynamic] : [Instant,]),]  #instant charge-charge interaction
    return DiagParaF64(
        type=Ver4Diag,
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

function diagram(paramc::ParaMC, _partition::Vector{T};
    channel=[PHr, PHEr, PPr],
    filter=[
        NoHartree,
        # Girreducible,
        # Proper,   #one interaction irreduble diagrams or not
        # NoBubble, #allow the bubble diagram or not
    ],
    dR=false # whether to use the derivative of the renormalized interaction, for the RG purpose
) where {T}
    # println("Build the vertex4 diagrams into an expression tree ...")
    # _partition = UEG.partition(order)
    # println("Diagram set: ", _partition)

    KinL, KoutL, KinR = zeros(16), zeros(16), zeros(16)
    KinL[1], KoutL[2], KinR[3] = 1.0, 1.0, 1.0
    legK = [KinL, KoutL, KinR]

    diag = Vector{ExprTreeF64}()
    diagpara = Vector{DiagParaF64}()
    partition = Vector{T}()
    for p in _partition
        para = diagPara(paramc, p[1], filter, KinL - KoutL)
        legK = [DiagTree.getK(para.totalLoopNum + 3, 1), DiagTree.getK(para.totalLoopNum + 3, 2), DiagTree.getK(para.totalLoopNum + 3, 3)]
        d::Vector{Diagram{Float64}} = Parquet.vertex4(para, legK, channel).diagram
        d = DiagTree.derivative(d, BareGreenId, p[2], index=1)
        d = DiagTree.derivative(d, BareInteractionId, p[3], index=2)

        # the Taylor expansion should be d^n f(x) / dx^n / n!, so there is a factor of 1/n! for each derivative
        for _d in d
            _d.factor *= 1 / factorial(p[2]) / factorial(p[3])
        end

        if dR
            d = DiagTree.derivative(d, BareInteractionId, 1, index=3)
        end
        if isempty(d) == false
            if paramc.isFock # remove the Fock subdiagrams
                DiagTree.removeHartreeFock!(d)
            end
            push!(diagpara, para)
            push!(diag, ExprTree.build(d))
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

    # diag = [ExprTree.build(d) for d in ver4]    #expression tree representation of diagrams 
    rootuu = [[idx for idx in d.root if d.node.object[idx].para.response == UpUp] for d in diag] #select the diagram with upup
    rootud = [[idx for idx in d.root if d.node.object[idx].para.response == UpDown] for d in diag] #select the diagram with updown
    #assign the external Tau to the corresponding diagrams
    extTuu = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootuu)]
    extTud = [[diag[ri].node.object[idx].para.extT for idx in root] for (ri, root) in enumerate(rootud)]
    return (partition, diagpara, diag, [rootuu, rootud], [extTuu, extTud])
end


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

@inline ud2sa(Wuu, Wud) = @. (Wuu + Wud) / 2, (Wuu - Wud) / 2
@inline sa2ud(Ws, Wa) = @. Ws + Wa, Ws - Wa

# include("common.jl")

# include("ver4_avg.jl")
include("ver4KW.jl")
include("ver4_PH_l.jl")
include("ver4_PH_l_df.jl")
include("exchange_interaction.jl")

end